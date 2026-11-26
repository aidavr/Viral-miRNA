#!/usr/bin/env python3
"""
A1 · Descarga versionada de las fuentes.

Qué hace
--------
Lee config/fuentes.yml, descarga cada archivo en datos/crudo/<fuente>_<version>/
y escribe datos/crudo/MANIFEST.tsv con la procedencia de cada archivo:
URL, fecha de descarga, tamaño y sha256.

Por qué un MANIFEST
-------------------
Es lo que convierte una descarga en un dato citable. El sha256 es una
huella del contenido: si mañana la fuente cambia el archivo sin avisar,
la huella no coincide y nos enteramos. Sin esto, "descargado de miRBase"
no significa nada reproducible.

El script es idempotente: si el archivo ya está y su huella coincide con
la del MANIFEST, no vuelve a descargarlo. Con --forzar se baja todo otra vez.

Uso
---
    python scripts/01_descargar.py
    python scripts/01_descargar.py --solo mirbase
    python scripts/01_descargar.py --forzar
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import yaml

RAIZ = Path(__file__).resolve().parent.parent
CONFIG = RAIZ / "config" / "fuentes.yml"
DIR_CRUDO = RAIZ / "datos" / "crudo"
MANIFEST = DIR_CRUDO / "MANIFEST.tsv"

# Identificarse es buena práctica al descargar de servidores públicos.
AGENTE = "venomir-v2 (pipeline de referencia de vmiRNAs; contacto: ver README)"
ESPERA = 2.5        # segundos entre descargas, para no saturar el servidor
REINTENTOS = 5


# ---------------------------------------------------------------------
# utilidades
# ---------------------------------------------------------------------
def sha256_de(ruta: Path) -> str:
    """Huella del contenido del archivo (por bloques, para no cargarlo entero)."""
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def descargar(url: str, destino: Path) -> None:
    """Descarga con reintentos. Escribe primero en .parcial para no dejar
    archivos a medias si algo falla."""
    destino.parent.mkdir(parents=True, exist_ok=True)
    tmp = destino.with_suffix(destino.suffix + ".parcial")
    ultimo_error: Exception | None = None
    for intento in range(1, REINTENTOS + 1):
        try:
            peticion = urllib.request.Request(url, headers={"User-Agent": AGENTE})
            with urllib.request.urlopen(peticion, timeout=120) as r, open(tmp, "wb") as f:
                while True:
                    bloque = r.read(1 << 16)
                    if not bloque:
                        break
                    f.write(bloque)
            tmp.replace(destino)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            ultimo_error = e
            print(f"      intento {intento}/{REINTENTOS} falló: {e}")
            time.sleep(6 * intento)   # espera creciente
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f"no se pudo descargar {url}: {ultimo_error}")


def leer_manifest() -> dict[str, dict]:
    if not MANIFEST.exists():
        return {}
    with open(MANIFEST, encoding="utf-8", newline="") as f:
        return {fila["ruta"]: fila for fila in csv.DictReader(f, delimiter="\t")}


def escribir_manifest(filas: dict[str, dict]) -> None:
    campos = ["ruta", "fuente", "version", "url", "fecha_descarga", "bytes", "sha256"]
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos, delimiter="\t")
        w.writeheader()
        for ruta in sorted(filas):
            w.writerow({c: filas[ruta].get(c, "") for c in campos})


# ---------------------------------------------------------------------
# Wayback Machine
# ---------------------------------------------------------------------
# VIRmiRNA ya no está en línea. Wayback guarda copias con fecha; cada
# copia tiene una URL permanente del tipo
#     http://web.archive.org/web/<timestamp>id_/<url original>
# El sufijo "id_" pide el contenido ORIGINAL, sin la barra de Wayback.
# Pedimos el snapshot más cercano a una fecha y anotamos cuál nos dio,
# de modo que la descarga se puede repetir exactamente igual.
CDX = "https://web.archive.org/cdx/search/cdx"


def snapshot_mas_cercano(url: str, fecha_objetivo: str) -> str | None:
    """Devuelve el timestamp del snapshot con código 200 más cercano."""
    consulta = urllib.parse.urlencode(
        {"url": url, "output": "json", "fl": "timestamp,statuscode", "limit": "200"}
    )
    peticion = urllib.request.Request(f"{CDX}?{consulta}", headers={"User-Agent": AGENTE})
    with urllib.request.urlopen(peticion, timeout=120) as r:
        datos = json.loads(r.read().decode("utf-8") or "[]")
    if len(datos) < 2:
        return None
    validos = [fila[0] for fila in datos[1:] if fila[1] in ("200", "-")]
    if not validos:
        return None
    return min(validos, key=lambda ts: abs(int(ts[:8]) - int(fecha_objetivo)))


def descargar_virmirna(fuente: dict, filas: dict, forzar: bool) -> list[str]:
    """Descarga una página archivada por virus.

    Wayback corta la conexión si se le piden muchas cosas seguidas, así que
    un fallo NO aborta la descarga entera: se anota el virus y se sigue.
    Como el script es idempotente, basta con volver a ejecutarlo para
    completar lo que falte. Devuelve los códigos que fallaron.
    """
    cfg = fuente["wayback"]
    carpeta = DIR_CRUDO / f"virmirna_{fuente['version']}"
    hoy = date.today().isoformat()
    fallidos: list[str] = []
    for codigo in cfg["codigos"]:
        url = cfg["url_patron"].format(codigo=codigo)
        destino = carpeta / f"browse_{codigo}.html"
        clave = str(destino.relative_to(RAIZ)).replace("\\", "/")
        if destino.exists() and not forzar and clave in filas:
            print(f"   = {codigo}")
            continue
        try:
            ts = snapshot_mas_cercano(url, cfg["fecha_objetivo"])
        except Exception as e:
            print(f"   ! {codigo}: no se pudo consultar el archivo ({e})")
            fallidos.append(codigo)
            time.sleep(ESPERA * 3)
            continue
        if ts is None:
            print(f"   ! {codigo}: sin copia archivada")
            continue
        url_wb = f"https://web.archive.org/web/{ts}id_/{url}"
        print(f"   + {codigo}  (snapshot {ts})")
        try:
            descargar(url_wb, destino)
        except RuntimeError as e:
            print(f"   ! {codigo}: {e}")
            fallidos.append(codigo)
            time.sleep(ESPERA * 3)
            continue
        filas[clave] = {
            "ruta": clave,
            "fuente": "virmirna",
            "version": fuente["version"],
            "url": url_wb,          # se guarda la URL CON el timestamp: es permanente
            "fecha_descarga": hoy,
            "bytes": destino.stat().st_size,
            "sha256": sha256_de(destino),
        }
        time.sleep(ESPERA)
    return fallidos


# ---------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Descarga las fuentes de VenomiR v2")
    ap.add_argument("--solo", help="descargar solo esta fuente (p. ej. mirbase)")
    ap.add_argument("--forzar", action="store_true", help="volver a descargar todo")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    filas = leer_manifest()
    hoy = date.today().isoformat()
    pendientes: list[str] = []
    incompletos: list[str] = []

    for fuente in cfg["fuentes"]:
        nombre, version = fuente["nombre"], str(fuente["version"])
        etiqueta = f"{nombre} {version}"
        if args.solo and args.solo != nombre:
            continue
        if fuente.get("estado") == "pendiente":
            pendientes.append(f"{etiqueta}: {fuente['descripcion']}")
            continue

        print(f"\n== {etiqueta} ==")
        if "wayback" in fuente:
            fallidos = descargar_virmirna(fuente, filas, args.forzar)
            escribir_manifest(filas)          # guardar lo ya conseguido
            incompletos.extend(f"{etiqueta}: {c}" for c in fallidos)
            continue

        carpeta = DIR_CRUDO / f"{nombre}_{version}"
        for archivo in fuente.get("archivos", []):
            destino = carpeta / archivo["nombre"]
            clave = str(destino.relative_to(RAIZ)).replace("\\", "/")
            if destino.exists() and not args.forzar and clave in filas:
                if sha256_de(destino) == filas[clave]["sha256"]:
                    print(f"   = {archivo['nombre']}")
                    continue
                print(f"   ! {archivo['nombre']}: la huella no coincide, se baja de nuevo")
            print(f"   + {archivo['nombre']}")
            descargar(archivo["url"], destino)
            filas[clave] = {
                "ruta": clave,
                "fuente": nombre,
                "version": version,
                "url": archivo["url"],
                "fecha_descarga": hoy,
                "bytes": destino.stat().st_size,
                "sha256": sha256_de(destino),
            }
            time.sleep(ESPERA)

    escribir_manifest(filas)
    print(f"\nMANIFEST: {MANIFEST}  ({len(filas)} archivos)")
    if pendientes:
        print("\nFuentes pendientes de acceso (no descargadas):")
        for p in pendientes:
            print(f"   - {p}")
    if incompletos:
        print(f"\n{len(incompletos)} descargas fallaron. El archivo corta la conexión "
              f"cuando se le pide mucho seguido: vuelve a ejecutar el script y solo "
              f"bajará lo que falta.")
        for c in incompletos:
            print(f"   - {c}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
