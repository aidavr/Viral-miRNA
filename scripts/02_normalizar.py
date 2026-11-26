#!/usr/bin/env python3
"""
A2 · Normalización: cada fuente, a un formato común.

Qué hace
--------
Lee lo descargado en datos/crudo/ y escribe datos/normalizado/<fuente>_<version>.tsv
con SIEMPRE las mismas columnas:

    id_fuente  nombre  nombre_norm  tipo  secuencia  especie_texto  extra

Por qué este paso existe
------------------------
Cada base publica en un formato distinto (FASTA, tabla HTML, Excel...).
Si el cotejo tuviera que entender todos los formatos a la vez sería
imposible de mantener. Separando "leer cada formato" (aquí) de "decidir
qué es lo mismo" (03_cotejar.py), añadir una fuente nueva solo obliga a
escribir una función de lectura.

Dos normalizaciones importantes
-------------------------------
1. SECUENCIA: mayúsculas y en RNA (T -> U). Sin esto, 'ugagguag' y
   'TGAGGTAG' parecerían moléculas distintas siendo la misma.
2. nombre_norm: minúsculas y '_' -> '-'. Solo sirve para DETECTAR
   nombres que son el mismo salvo formato (caso 6, p. ej.
   rlcv-miR-rL1-21-3p frente a rlcv-miR-rl1-21-3p). El nombre original
   NUNCA se pierde: se guarda en `nombre`.

Uso
---
    python scripts/02_normalizar.py
"""

from __future__ import annotations

import csv
import html
import json
import re
import sys
from pathlib import Path

import yaml

RAIZ = Path(__file__).resolve().parent.parent
CONFIG = RAIZ / "config" / "fuentes.yml"
DIR_CRUDO = RAIZ / "datos" / "crudo"
DIR_NORM = RAIZ / "datos" / "normalizado"

COLUMNAS = ["id_fuente", "nombre", "nombre_norm", "tipo", "secuencia",
            "especie_texto", "extra"]


# ---------------------------------------------------------------------
# normalizadores
# ---------------------------------------------------------------------
def norm_secuencia(s: str) -> tuple[str, str]:
    """Normaliza una secuencia a RNA en mayúsculas.

    Devuelve (secuencia, nota). La nota queda vacía si no hubo nada raro.

    Regla de oro: NUNCA descartar en silencio. En VIRmiRNA hay 5 secuencias
    escritas con notación especial y, si se tiran sin avisar, el recuento ya
    no cuadra con el artículo y nadie sabe por qué:
      · 'uuaguguu...gagugg(c)'  el paréntesis marca una base final opcional
        -> se quita y se anota (la forma con esa base es una isoforma, que
           se detecta sola en el caso 7)
      · 'cccggagc...guucgb'     'b' es un código IUPAC ambiguo (C, G o U)
        -> pasa a N, que es lo que entienden bowtie y miRDeep2
    """
    s = re.sub(r"\s+", "", s or "").upper().replace("T", "U")
    if not s:
        return "", ""
    notas = []
    if "(" in s:
        s = re.sub(r"\([^)]*\)", "", s)
        notas.append("base final opcional entre paréntesis, eliminada")
    ambiguas = set(re.findall(r"[^ACGUN]", s))
    if ambiguas:
        if not ambiguas <= set("RYSWKMBDHV"):      # códigos IUPAC válidos
            return "", ""                           # esto no es una secuencia
        s = re.sub(r"[^ACGUN]", "N", s)
        notas.append(f"códigos IUPAC ambiguos ({','.join(sorted(ambiguas))}) -> N")
    return s, "; ".join(notas)


def norm_nombre(n: str) -> str:
    """Clave para comparar nombres que solo difieren en formato."""
    n = (n or "").strip().lstrip(">").lower().replace("_", "-")
    return re.sub(r"-{2,}", "-", n).strip("-")


def quitar_html(texto: str) -> str:
    """miRBase sirve los .fa envueltos en HTML (comprobado el 25/11/2026).
    Deshace <br>, quita etiquetas y traduce &amp; y compañía."""
    texto = re.sub(r"(?is)<script.*?</script>|<style.*?</style>", "", texto)
    texto = re.sub(r"(?i)<br\s*/?>", "\n", texto)
    texto = re.sub(r"<[^>]+>", "", texto)
    return html.unescape(texto)


# ---------------------------------------------------------------------
# lectores, uno por formato
# ---------------------------------------------------------------------
def especie_de_cabecera(campos: list[str]) -> str:
    """Saca el nombre de especie de una cabecera de miRBase.

    El formato cambia entre los dos archivos, y es un fallo fácil de cometer:
        mature.fa   >ebv-miR-BHRF1-1-5p MIMAT0000995 Epstein Barr virus miR-BHRF1-1-5p
        hairpin.fa  >ebv-mir-BHRF1-1    MI0001064    Epstein Barr virus miR-BHRF1-1 stem-loop
                                                     \\_____ especie ____/ \\__ sobra __/
    El nombre de especie ocupa un número variable de palabras ('Homo sapiens'
    son 2, 'Epstein Barr virus' son 3), así que no se puede cortar por
    posición fija: hay que quitar por el final 'stem-loop' (si está) y el
    nombre corto del miRNA, y quedarse con lo de en medio.
    """
    resto = campos[2:]
    if resto and resto[-1] == "stem-loop":
        resto = resto[:-1]
    return " ".join(resto[:-1]) if len(resto) > 1 else ""


def leer_fasta_mirbase(ruta: Path, tipo: str) -> list[dict]:
    """FASTA de miRBase (ver especie_de_cabecera para el formato)."""
    texto = quitar_html(ruta.read_text(encoding="utf-8", errors="replace"))
    registros, nombre, acc, especie, piezas = [], None, "", "", []

    def cerrar():
        if nombre:
            seq, nota = norm_secuencia("".join(piezas))
            registros.append({
                "id_fuente": acc, "nombre": nombre, "nombre_norm": norm_nombre(nombre),
                "tipo": tipo, "secuencia": seq, "especie_texto": especie,
                "extra": json.dumps({"nota_secuencia": nota}, ensure_ascii=False) if nota else "",
            })

    for linea in texto.splitlines():
        linea = linea.strip()
        if not linea:
            continue
        if linea.startswith(">"):
            cerrar()
            campos = linea[1:].split()
            nombre = campos[0] if campos else ""
            acc = campos[1] if len(campos) > 1 and campos[1].startswith("MI") else ""
            especie = especie_de_cabecera(campos)
            piezas = []
        else:
            piezas.append(linea)
    cerrar()
    return [r for r in registros if r["nombre"] and r["secuencia"]]


def unificar_especie_por_prefijo(registros: list[dict]) -> list[dict]:
    """Corrige las especies que salen cortadas de la cabecera.

    Unas pocas cabeceras de miRBase no llevan el nombre corto del miRNA, así
    que al quitar el último campo se pierde una palabra de la especie y
    aparecen variantes como 'Mouse gammaherpesvirus' junto a
    'Mouse gammaherpesvirus 68'.

    Todos los miRNAs de un mismo organismo comparten el prefijo del nombre
    ('mghv-', 'hsa-'), que miRBase asigna de forma estable. Así que para cada
    prefijo se toma la forma de especie MÁS FRECUENTE y se aplica a todos.
    Es una solución provisional: en el paso A6 la especie vendrá del taxid
    de NCBI, que es el criterio correcto.
    """
    from collections import Counter, defaultdict
    por_prefijo: dict[str, Counter] = defaultdict(Counter)
    for r in registros:
        prefijo = r["nombre"].split("-")[0].lower()
        if r["especie_texto"]:
            por_prefijo[prefijo][r["especie_texto"]] += 1
    mayoritaria = {p: c.most_common(1)[0][0] for p, c in por_prefijo.items()}
    for r in registros:
        prefijo = r["nombre"].split("-")[0].lower()
        if prefijo in mayoritaria:
            r["especie_texto"] = mayoritaria[prefijo]
    return registros


def leer_tabla_virmirna(ruta: Path) -> list[dict]:
    """Tabla HTML de VIRmiRNA (snapshot archivado).
    Columnas: VMR_ID | miRNA | Virus | miRNA_Seq | Length | GC |
              seed miRBase | seed VIRmiRNA | Método | Diana | PMID
    Se parte por '<tr' porque las filas del servidor no cierran </tr>:
    un parser estricto no las vería.
    """
    bruto = ruta.read_text(encoding="utf-8", errors="replace")
    bruto = re.sub(r"(?is)<script.*?</script>|<style.*?</style>", "", bruto)
    registros = []
    for trozo in re.split(r"(?i)<tr[^>]*>", bruto)[1:]:
        celdas = [
            " ".join(html.unescape(re.sub(r"<[^>]+>", " ", c)).split())
            for c in re.findall(r"(?is)<td[^>]*>(.*?)(?=</td>|<td|</tr|$)", trozo)
        ]
        if len(celdas) < 4 or not celdas[0].startswith("VMR_"):
            continue
        seq, nota = norm_secuencia(celdas[3])
        if not seq:
            print(f"   ! {celdas[0]} ({celdas[1]}): secuencia ilegible {celdas[3]!r}")
            continue
        extra = {"virus_fuente": celdas[2]}
        if nota:
            extra["nota_secuencia"] = nota
            extra["secuencia_original"] = celdas[3]
        if len(celdas) >= 11:
            extra["metodo"] = celdas[8]
            extra["diana"] = celdas[9]
            extra["pmid"] = celdas[10]
        registros.append({
            "id_fuente": celdas[0], "nombre": celdas[1],
            "nombre_norm": norm_nombre(celdas[1]), "tipo": "maduro",
            "secuencia": seq, "especie_texto": celdas[2],
            "extra": json.dumps(extra, ensure_ascii=False),
        })
    return registros


# ---------------------------------------------------------------------
def escribir(registros: list[dict], destino: Path) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    with open(destino, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS, delimiter="\t")
        w.writeheader()
        for r in registros:
            w.writerow({c: r.get(c, "") for c in COLUMNAS})


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    resumen = []

    for fuente in cfg["fuentes"]:
        nombre, version = fuente["nombre"], str(fuente["version"])
        if fuente.get("estado") == "pendiente":
            continue
        carpeta = DIR_CRUDO / f"{nombre}_{version}"
        if not carpeta.exists():
            print(f"!! falta {carpeta}: ejecuta antes 01_descargar.py")
            continue

        registros: list[dict] = []
        if nombre == "mirbase":
            for archivo, tipo in (("mature.fa", "maduro"), ("hairpin.fa", "horquilla")):
                ruta = carpeta / archivo
                if ruta.exists():
                    registros += leer_fasta_mirbase(ruta, tipo)
            # se unifica con los DOS archivos juntos: las cabeceras de
            # mature.fa son fiables y corrigen las de hairpin.fa
            registros = unificar_especie_por_prefijo(registros)
        elif nombre == "virmirna":
            for ruta in sorted(carpeta.glob("browse_*.html")):
                registros += leer_tabla_virmirna(ruta)
        else:
            print(f"!! sin lector para {nombre}")
            continue

        destino = DIR_NORM / f"{nombre}_{version}.tsv"
        escribir(registros, destino)
        n_sec = len({r["secuencia"] for r in registros})
        print(f"{nombre} {version}: {len(registros)} registros, {n_sec} secuencias únicas"
              f"  ->  {destino.relative_to(RAIZ)}")
        resumen.append((nombre, version, len(registros), n_sec))

    if not resumen:
        print("No se normalizó nada.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
