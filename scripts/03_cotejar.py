#!/usr/bin/env python3
"""
A3-A4 · Cotejo: construir la base SQLite y resolver nombres contra secuencias.

Qué hace
--------
1. Crea datos/venomir_v2.sqlite con el esquema de sql/esquema.sql
2. Carga los TSV de datos/normalizado/ en la tabla `entrada`
3. Deriva `molecula` (una fila por secuencia única) y `alias`
4. Detecta y anota los casos de cotejo en la tabla `marca`
5. Escribe datos/casos_cotejo.tsv con el recuento de cada caso

Idea que hay que tener clara
----------------------------
`entrada` es sagrada: es lo que dijo cada fuente. Todo lo demás
(molecula, alias, marca) se BORRA y se vuelve a calcular en cada
ejecución. Así, cambiar una regla de cotejo no obliga a descargar nada
otra vez, y el resultado siempre es consecuencia de reglas escritas,
no de arreglos a mano.

Uso
---
    python scripts/03_cotejar.py
    python scripts/03_cotejar.py --bd otra_ruta.sqlite
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from datetime import date
from pathlib import Path

import yaml

RAIZ = Path(__file__).resolve().parent.parent
CONFIG = RAIZ / "config" / "fuentes.yml"
ESQUEMA = RAIZ / "sql" / "esquema.sql"
DIR_NORM = RAIZ / "datos" / "normalizado"
BD_POR_DEFECTO = RAIZ / "datos" / "venomir_v2.sqlite"
MANIFEST = RAIZ / "datos" / "crudo" / "MANIFEST.tsv"


# ---------------------------------------------------------------------
def abrir_bd(ruta: Path) -> sqlite3.Connection:
    ruta.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(ruta)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")    # SQLite las ignora si no se pide
    con.executescript(ESQUEMA.read_text(encoding="utf-8"))
    return con


def cargar_fuentes(con: sqlite3.Connection, cfg: dict) -> dict[tuple[str, str], int]:
    """Inserta en `fuente` lo que dice la configuración, más el sha256 del
    MANIFEST cuando lo hay. Devuelve {(nombre, version): id}."""
    huellas: dict[tuple[str, str], str] = {}
    if MANIFEST.exists():
        with open(MANIFEST, encoding="utf-8", newline="") as f:
            for fila in csv.DictReader(f, delimiter="\t"):
                # una fuente tiene varios archivos; basta con guardar uno
                huellas.setdefault((fila["fuente"], fila["version"]), fila["sha256"])

    ids = {}
    hoy = date.today().isoformat()
    for f in cfg["fuentes"]:
        if f.get("estado") == "pendiente":
            continue
        clave = (f["nombre"], str(f["version"]))
        con.execute(
            """INSERT OR IGNORE INTO fuente
               (nombre, version, fecha_fuente, fecha_descarga, url, sha256, prioridad, cita)
               VALUES (?,?,?,?,?,?,?,?)""",
            (clave[0], clave[1], str(f.get("fecha", "")), hoy,
             f.get("url_base", "") or f.get("wayback", {}).get("url_patron", ""),
             huellas.get(clave, ""), f["prioridad"], f.get("cita", "")),
        )
        ids[clave] = con.execute(
            "SELECT id FROM fuente WHERE nombre=? AND version=?", clave
        ).fetchone()["id"]
    return ids


def cargar_entradas(con: sqlite3.Connection, ids: dict, universo: dict) -> None:
    """Carga los TSV normalizados, quedándose solo con el universo definido
    en config/fuentes.yml (humano + virus). Los archivos de datos/normalizado/
    siguen completos: el filtro se aplica aquí, así que cambiar la regla y
    volver a ejecutar es inmediato."""
    exactas = set(universo.get("especies_exactas", []))
    patron = re.compile(universo.get("especies_patron", "."), re.IGNORECASE)

    def en_universo(especie: str) -> bool:
        # VIRmiRNA es una base de miRNAs virales: todo lo suyo entra aunque
        # escriba el nombre de especie de otra manera ('HSV1', 'EBV'...).
        return especie in exactas or bool(patron.search(especie or ""))

    con.execute("DELETE FROM entrada")
    especies_aceptadas: set[str] = set()
    for (nombre, version), fuente_id in ids.items():
        ruta = DIR_NORM / f"{nombre}_{version}.tsv"
        if not ruta.exists():
            print(f"   ! falta {ruta.name}")
            continue
        es_base_viral = nombre in ("virmirna", "antivirmir")
        filas, descartadas = [], 0
        with open(ruta, encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if not (es_base_viral or en_universo(r["especie_texto"])):
                    descartadas += 1
                    continue
                especies_aceptadas.add(r["especie_texto"])
                filas.append((fuente_id, r["id_fuente"], r["nombre"], r["nombre_norm"],
                              r["tipo"], r["secuencia"] or None, r["especie_texto"],
                              r["extra"]))
        con.executemany(
            """INSERT OR IGNORE INTO entrada
               (fuente_id, id_fuente, nombre, nombre_norm, tipo, secuencia, especie_texto, extra)
               VALUES (?,?,?,?,?,?,?,?)""", filas)
        print(f"   {nombre} {version}: {len(filas)} entradas "
              f"({descartadas} fuera del universo)")

    destino = RAIZ / "datos" / "universo_especies.tsv"
    destino.write_text(
        "especie\n" + "\n".join(sorted(e for e in especies_aceptadas if e)) + "\n",
        encoding="utf-8")
    print(f"   especies en el universo: {len(especies_aceptadas)}"
          f"  ->  {destino.relative_to(RAIZ)}")


# ---------------------------------------------------------------------
# EL COTEJO
# ---------------------------------------------------------------------
def derivar(con: sqlite3.Connection) -> None:
    """molecula y alias se calculan enteros desde `entrada`."""
    for tabla in ("marca", "alias", "molecula_taxon", "madura_en_horquilla", "molecula"):
        con.execute(f"DELETE FROM {tabla}")

    # --- molecula: una fila por secuencia única -----------------------
    # Aquí ocurre el cotejo: dos entradas con la misma secuencia acaban
    # en la MISMA molécula, aunque se llamen distinto y vengan de bases
    # distintas. El tipo (maduro/horquilla) forma parte de la identidad.
    con.execute("""
        INSERT INTO molecula (secuencia, longitud, tipo)
        SELECT secuencia, LENGTH(secuencia), tipo
        FROM entrada
        WHERE secuencia IS NOT NULL
        GROUP BY secuencia, tipo
    """)

    # --- alias: todos los nombres de cada molécula --------------------
    con.execute("""
        INSERT OR IGNORE INTO alias (molecula_id, fuente_id, nombre, nombre_norm, id_fuente)
        SELECT m.id, e.fuente_id, e.nombre, e.nombre_norm, e.id_fuente
        FROM entrada e
        JOIN molecula m ON m.secuencia = e.secuencia AND m.tipo = e.tipo
    """)

    # --- nombre principal: gana la fuente de menor `prioridad` --------
    # Criterio escrito en config/fuentes.yml, no decidido aquí.
    # En caso de empate dentro de la misma fuente, el nombre alfabético
    # menor, solo para que el resultado sea determinista.
    con.execute("""
        UPDATE alias SET es_principal = 1
        WHERE id IN (
            SELECT a.id FROM alias a
            JOIN fuente f ON f.id = a.fuente_id
            WHERE a.id = (
                SELECT a2.id FROM alias a2
                JOIN fuente f2 ON f2.id = a2.fuente_id
                WHERE a2.molecula_id = a.molecula_id
                ORDER BY f2.prioridad, a2.nombre
                LIMIT 1)
        )
    """)
    con.commit()


def marcar_casos(con: sqlite3.Connection) -> None:
    """Cada caso del diseño, con su consulta. Todas escriben en `marca`."""

    # CASO 1 · misma secuencia, nombres distintos (ignorando el formato).
    # Se cuenta por nombre_norm: si solo cambian mayúsculas o guiones es
    # el caso 6, no el 1.
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle, valor)
        SELECT a.molecula_id, 'caso1_multinombre',
               GROUP_CONCAT(DISTINCT a.nombre),
               COUNT(DISTINCT a.nombre_norm)
        FROM alias a
        GROUP BY a.molecula_id
        HAVING COUNT(DISTINCT a.nombre_norm) > 1
    """)

    # CASO 6 · el nombre solo cambia en formato (mayúsculas, guion bajo).
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle, valor)
        SELECT a.molecula_id, 'caso6_formato_nombre',
               GROUP_CONCAT(DISTINCT a.nombre), COUNT(DISTINCT a.nombre)
        FROM alias a
        GROUP BY a.molecula_id, a.nombre_norm
        HAVING COUNT(DISTINCT a.nombre) > 1
    """)

    # CASO 4 · mismo nombre, secuencias distintas.
    # Es el inverso del caso 1 y hay que mirarlo: puede ser un error de
    # la fuente, una isoforma mal nombrada o un cambio entre versiones.
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle, valor)
        SELECT a.molecula_id, 'caso4_conflicto_nombre', a.nombre_norm, n.n
        FROM alias a
        JOIN (SELECT nombre_norm, COUNT(DISTINCT molecula_id) AS n
              FROM alias GROUP BY nombre_norm HAVING n > 1) n
          ON n.nombre_norm = a.nombre_norm
    """)

    # CASO 5 · renombrado entre versiones de la MISMA base de datos.
    # De aquí sale la tabla de equivalencias que traduce la matriz vieja.
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle)
        SELECT v.molecula_id, 'caso5_renombrado',
               v.antiguo || ' -> ' || v.nuevo
        FROM (
            SELECT a1.molecula_id,
                   a1.nombre AS antiguo, a2.nombre AS nuevo
            FROM alias a1
            JOIN fuente f1 ON f1.id = a1.fuente_id
            JOIN alias a2 ON a2.molecula_id = a1.molecula_id
            JOIN fuente f2 ON f2.id = a2.fuente_id
            WHERE f1.nombre = f2.nombre          -- misma base de datos
              AND f1.version <> f2.version       -- distinta versión
              AND f1.version < f2.version        -- antigua -> nueva
              AND a1.nombre <> a2.nombre
        ) v
    """)

    # CASO 12 · el mismo miRNA repetido DENTRO de una misma fuente, con
    # dos identificadores distintos. No estaba previsto en el diseño:
    # apareció al validar contra la Tabla S1 del artículo (VIRmiRNA tiene
    # 5, probablemente el mismo miRNA descrito en dos publicaciones).
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle, valor)
        SELECT a.molecula_id, 'caso12_duplicado_fuente',
               a.nombre || ' (' || GROUP_CONCAT(DISTINCT e.id_fuente) || ')',
               COUNT(DISTINCT e.id_fuente)
        FROM alias a
        JOIN entrada e ON e.fuente_id = a.fuente_id AND e.nombre = a.nombre
        JOIN molecula m ON m.id = a.molecula_id
                       AND m.secuencia = e.secuencia AND m.tipo = e.tipo
        GROUP BY a.molecula_id, a.fuente_id, a.nombre
        HAVING COUNT(DISTINCT e.id_fuente) > 1
    """)

    # CASO 2 (provisional) · la misma secuencia aparece atribuida a más
    # de una especie. Provisional porque aquí se compara el TEXTO de la
    # especie tal como lo escribe cada fuente; se resolverá de verdad en
    # el paso A6, con los taxid de NCBI.
    con.execute("""
        INSERT INTO marca (molecula_id, tipo, detalle, valor)
        SELECT m.id, 'caso2_multiespecie',
               GROUP_CONCAT(DISTINCT e.especie_texto), COUNT(DISTINCT e.especie_texto)
        FROM molecula m
        JOIN entrada e ON e.secuencia = m.secuencia AND e.tipo = m.tipo
        WHERE e.especie_texto <> ''
        GROUP BY m.id
        HAVING COUNT(DISTINCT e.especie_texto) > 1
    """)

    con.commit()


def marcar_isoformas(con: sqlite3.Connection) -> None:
    """CASO 7 · una secuencia contenida en otra (isoformas).

    No se hace en SQL porque comparar "una cadena dentro de otra" con LIKE
    obliga a recorrer todas las parejas y es lentísimo. En Python se agrupa
    por longitud y se comprueba solo lo necesario.
    """
    filas = con.execute(
        "SELECT id, secuencia FROM molecula WHERE tipo='maduro'").fetchall()
    secuencias = [(r["id"], r["secuencia"]) for r in filas]
    # se ordena de corta a larga: solo puede estar contenida la más corta
    secuencias.sort(key=lambda x: len(x[1]))
    largas = [(i, s) for i, s in secuencias]
    marcas = []
    for pos, (id_corta, corta) in enumerate(secuencias):
        for id_larga, larga in largas[pos + 1:]:
            if len(larga) == len(corta):
                continue
            if corta in larga:
                marcas.append((id_corta, "caso7_isoforma",
                               f"contenida en la molécula {id_larga}", id_larga))
    con.executemany(
        """INSERT INTO marca (molecula_id, tipo, detalle, ref_molecula_id)
           VALUES (?,?,?,?)""", marcas)
    con.commit()
    print(f"   caso 7 (isoformas): {len(marcas)} relaciones")


# ---------------------------------------------------------------------
def informe(con: sqlite3.Connection) -> None:
    destino = RAIZ / "datos" / "casos_cotejo.tsv"
    filas = con.execute("SELECT * FROM v_casos").fetchall()
    with open(destino, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["caso", "n_moleculas"])
        for r in filas:
            w.writerow([r["tipo"], r["n_moleculas"]])

    print("\n--- resumen ---")
    for tabla in ("entrada", "molecula", "alias", "marca"):
        n = con.execute(f"SELECT COUNT(*) AS n FROM {tabla}").fetchone()["n"]
        print(f"   {tabla:<10} {n:>7}")
    print("\n--- casos de cotejo ---")
    for r in filas:
        print(f"   {r['tipo']:<24} {r['n_moleculas']:>6}")
    print(f"\nEscrito: {destino.relative_to(RAIZ)}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Construye la base VenomiR v2")
    ap.add_argument("--bd", type=Path, default=BD_POR_DEFECTO)
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    con = abrir_bd(args.bd)
    print(f"Base: {args.bd}")

    print("\n1. fuentes")
    ids = cargar_fuentes(con, cfg)
    print("\n2. entradas")
    cargar_entradas(con, ids, cfg.get('universo', {}))
    con.commit()
    print("\n3. moléculas y alias")
    derivar(con)
    print("\n4. casos de cotejo")
    marcar_casos(con)
    marcar_isoformas(con)

    con.execute(
        "INSERT OR REPLACE INTO version_base (version, fecha, notas) VALUES (?,?,?)",
        (cfg["release"], date.today().isoformat(), "generada por scripts/03_cotejar.py"))
    con.commit()
    informe(con)
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
