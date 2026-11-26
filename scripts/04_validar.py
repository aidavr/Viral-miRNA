#!/usr/bin/env python3
"""
A8 · Validación: comprobar que la base está completa y es coherente.

Por qué este script existe
--------------------------
Un pipeline que "funciona" sin comprobaciones no sirve de nada: si una
descarga se queda a medias, el resultado sale igualmente, solo que con
menos datos y sin avisar. Aquí se comprueba contra una referencia
EXTERNA e independiente: la Tabla S1 del artículo de VIRmiRNA, que dice
cuántos miRNAs tiene cada virus. Si recuperamos 1.308 de 1.308, la
descarga está completa.

Comprobaciones
--------------
  1. VIRmiRNA: nº de miRNAs por virus frente a la Tabla S1 del artículo
  2. miRBase: los virales de la 23.1 deben ser 530 maduros / 320 horquillas
  3. Integridad: sin secuencias vacías, sin moléculas sin alias,
     un único nombre principal por molécula
  4. Coherencia: toda entrada con secuencia tiene su molécula

Devuelve código 0 si todo pasa, 1 si algo falla.

Uso
---
    python scripts/04_validar.py
"""

from __future__ import annotations

import csv
import sqlite3
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
BD = RAIZ / "datos" / "venomir_v2.sqlite"
TABLA_S1 = RAIZ / "datos" / "referencia" / "virmirna_tabla_s1.tsv"

# miRBase 23.1, contado al descargarla el 25/11/2026. Si cambia, es que
# la fuente ha cambiado y hay que mirarlo, no ajustar el número.
MIRBASE_VIRAL_ESPERADO = {"maduro": 530, "horquilla": 320}

fallos: list[str] = []
avisos: list[str] = []


def comprobar(condicion: bool, mensaje: str, critico: bool = True) -> None:
    if condicion:
        print(f"   OK   {mensaje}")
    elif critico:
        fallos.append(mensaje)
        print(f"   FALLA {mensaje}")
    else:
        avisos.append(mensaje)
        print(f"   aviso {mensaje}")


def main() -> int:
    if not BD.exists():
        print(f"No existe {BD}. Ejecuta antes scripts/03_cotejar.py")
        return 1
    con = sqlite3.connect(BD)
    con.row_factory = sqlite3.Row

    # --- 1. VIRmiRNA frente al artículo -------------------------------
    print("\n1. VIRmiRNA frente a la Tabla S1 del artículo")
    esperado = {}
    with open(TABLA_S1, encoding="utf-8", newline="") as f:
        for fila in csv.DictReader(f, delimiter="\t"):
            esperado[fila["codigo"]] = int(fila["n_mirnas"])

    obtenido = {
        r["especie_texto"]: r["n"]
        for r in con.execute("""
            SELECT e.especie_texto, COUNT(DISTINCT e.id_fuente) AS n
            FROM entrada e JOIN fuente f ON f.id = e.fuente_id
            WHERE f.nombre = 'virmirna'
            GROUP BY e.especie_texto""")
    }
    total_obt, total_esp = sum(obtenido.values()), sum(esperado.values())
    comprobar(total_obt == total_esp,
              f"total de miRNAs de VIRmiRNA: {total_obt} (esperado {total_esp})")
    for codigo, n in sorted(esperado.items()):
        if obtenido.get(codigo, 0) != n:
            comprobar(False, f"   {codigo}: {obtenido.get(codigo, 0)} (esperado {n})")
    comprobar(len(obtenido) == len(esperado),
              f"nº de virus: {len(obtenido)} (esperado {len(esperado)})")

    # --- 2. miRBase viral ---------------------------------------------
    print("\n2. miRBase 23.1, parte viral")
    for tipo, n_esp in MIRBASE_VIRAL_ESPERADO.items():
        n = con.execute("""
            SELECT COUNT(*) AS n FROM entrada e JOIN fuente f ON f.id = e.fuente_id
            WHERE f.nombre='mirbase' AND f.version='23.1' AND e.tipo=?
              AND e.especie_texto <> 'Homo sapiens'""", (tipo,)).fetchone()["n"]
        comprobar(n == n_esp, f"{tipo}s virales: {n} (esperado {n_esp})")

    # --- 3. integridad interna ----------------------------------------
    print("\n3. Integridad de la base")
    n = con.execute(
        "SELECT COUNT(*) AS n FROM molecula WHERE secuencia IS NULL OR secuencia=''"
    ).fetchone()["n"]
    comprobar(n == 0, f"moléculas sin secuencia: {n}")

    n = con.execute("""SELECT COUNT(*) AS n FROM molecula m
                       WHERE NOT EXISTS (SELECT 1 FROM alias a WHERE a.molecula_id = m.id)"""
                    ).fetchone()["n"]
    comprobar(n == 0, f"moléculas sin ningún nombre: {n}")

    n = con.execute("""SELECT COUNT(*) AS n FROM (
                         SELECT molecula_id FROM alias WHERE es_principal = 1
                         GROUP BY molecula_id HAVING COUNT(*) <> 1)"""
                    ).fetchone()["n"]
    comprobar(n == 0, f"moléculas sin un único nombre principal: {n}")

    n = con.execute("""SELECT COUNT(*) AS n FROM molecula
                       WHERE secuencia GLOB '*[^ACGUN]*'""").fetchone()["n"]
    comprobar(n == 0, f"secuencias con caracteres no válidos: {n}")

    # --- 4. coherencia entrada <-> molecula ---------------------------
    print("\n4. Coherencia")
    n = con.execute("""SELECT COUNT(*) AS n FROM entrada e
                       WHERE e.secuencia IS NOT NULL
                         AND NOT EXISTS (SELECT 1 FROM molecula m
                                         WHERE m.secuencia = e.secuencia AND m.tipo = e.tipo)"""
                    ).fetchone()["n"]
    comprobar(n == 0, f"entradas con secuencia pero sin molécula: {n}")

    # Hay un alias por entrada, salvo por los registros duplicados dentro
    # de una misma fuente (caso 12): dos entradas con el mismo nombre y la
    # misma secuencia comparten un único alias. La diferencia tiene que
    # cuadrar exactamente con el número de duplicados.
    n_ent = con.execute("SELECT COUNT(*) AS n FROM entrada").fetchone()["n"]
    n_ali = con.execute("SELECT COUNT(*) AS n FROM alias").fetchone()["n"]
    n_dup = con.execute("""SELECT COALESCE(SUM(valor - 1), 0) AS n FROM marca
                           WHERE tipo = 'caso12_duplicado_fuente'""").fetchone()["n"]
    comprobar(n_ent - n_ali == n_dup,
              f"entradas ({n_ent}) - alias ({n_ali}) = duplicados en la fuente ({int(n_dup)})")

    # --- resumen ------------------------------------------------------
    print("\n" + "=" * 60)
    if fallos:
        print(f"{len(fallos)} comprobaciones FALLAN:")
        for f_ in fallos:
            print(f"   - {f_}")
    else:
        print("Todas las comprobaciones críticas pasan.")
    if avisos:
        print(f"\n{len(avisos)} avisos (no bloquean):")
        for a in avisos:
            print(f"   - {a}")
    con.close()
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
