# VenomiR v2

> El repositorio se llama **Viral-miRNA**; *VenomiR v2* es el nombre de la base de datos
> que construye, por continuidad con el pipeline VenomiR anterior.

Base de datos de referencia de microRNAs (miRNAs) **humanos y virales**, construida
cotejando varias fuentes **por secuencia**, no por nombre.

Sirve para cuantificar miRNAs virales en datos de small RNA-seq sin confundirlos con
miRNAs humanos parecidos: el problema que motiva todo esto.

> **Estado:** en construcción. Funcionan los pasos A1–A4 y A8 (descarga, normalización,
> cotejo y validación). Pendientes: horquillas desde genoma, taxonomía NCBI, marca de
> similitud con el genoma humano y exportación a FASTA/GFF3.

---

## El problema

En el análisis previo, un miRNA atribuido al herpesvirus del pavo (`hvt-miR-H14-3p`,
488.233 lecturas en el 100 % de las muestras) resultó ser **hsa-miR-221-3p**, un miRNA
humano. Y otro atribuido al virus de la mancha blanca del camarón era una repetición
`(GAA)n` del genoma humano.

Pasó por dos motivos:

1. Las lecturas solo se alinearon contra genomas **virales**. Si una lectura humana se
   parece a un miRNA viral, no tiene otro sitio al que ir: se cuenta como viral.
2. La referencia se construyó uniendo bases de datos **por nombre**. Cuando dos bases
   llaman distinto a la misma secuencia, esa molécula se cuenta dos veces; cuando dos
   virus comparten secuencia, no queda registro de la ambigüedad.

VenomiR v2 ataca el segundo problema (la referencia). El primero se resuelve en el
pipeline de muestras, mapeando humano y virus **a la vez**.

---

## La idea central: la unidad es la secuencia

```
        miRBase 23.1              VIRmiRNA
   ebv-miR-BART1-3p          rlcv-miR-rl1-6-3p
   UAGCACCGCUAUCCACUAUGUC    UAGCACCGCUAUCCACUAUGUC
                  \             /
                   \           /
              UNA SOLA MOLÉCULA
              + dos alias (uno por fuente)
              + marca: "caso 1, misma secuencia con nombres distintos"
```

Nada se descarta y nada se arregla a mano: lo dudoso se **anota** en la tabla `marca`,
y de ahí sale el informe de casos que acompaña a cada versión.

---

## Instalación

Hace falta Python 3.10 o superior y una única dependencia:

```bash
pip install pyyaml
```

o, con conda:

```bash
conda env create -f environment.yml
conda activate venomir-v2
```

## Uso

```bash
python scripts/01_descargar.py     # descarga las fuentes y escribe el MANIFEST
python scripts/02_normalizar.py    # cada fuente, al mismo formato
python scripts/03_cotejar.py       # construye la base SQLite y cruza las fuentes
python scripts/04_validar.py       # comprueba que está completa y es coherente
```

`01_descargar.py` es idempotente: si se corta a medias (el archivo de Internet Archive
limita las peticiones seguidas), basta con volver a ejecutarlo y solo baja lo que falta.

## Qué sale

| Archivo | Qué es |
|---|---|
| `datos/crudo/MANIFEST.tsv` | Procedencia: URL, fecha, tamaño y `sha256` de cada archivo descargado |
| `datos/normalizado/*.tsv` | Cada fuente en el mismo formato de columnas |
| `datos/venomir_v2.sqlite` | **La base**: moléculas, alias, marcas |
| `datos/casos_cotejo.tsv` | Recuento de cada caso de cotejo |
| `datos/universo_especies.tsv` | Las especies aceptadas, para revisarlas |

## Resultado de la ejecución actual

Universo: humano + virus (miRBase trae 323 especies; el resto se deja fuera porque
los miRNAs están muy conservados entre animales y eso tapa lo que interesa ver).

| | |
|---|---|
| Entradas (registros de las fuentes) | 11.874 |
| Moléculas (secuencias únicas) | 6.841 |
| Alias (nombres) | 11.869 |
| Especies | 87 |

Casos de cotejo detectados:

| Caso | Moléculas | Qué es |
|---|---|---|
| 1 · misma secuencia, distinto nombre | 1.071 | `ebv-miR-BART1-3p` = `rlcv-miR-rl1-6-3p` |
| 2 · misma secuencia, varias especies | 575 | provisional, por texto de especie |
| 4 · mismo nombre, distinta secuencia | 3.427 | casi todo, cambios de secuencia entre miRBase 22 y 23.1 |
| 5 · renombrado entre versiones | 851 | `hcmv-miR-UL148D` → `hcmv-miR-UL148D-3p` |
| 6 · el nombre solo cambia de formato | 269 | `rlcv-miR-rL1-21-3p` / `rlcv-miR-rl1-21-3p` |
| 7 · isoforma (una secuencia dentro de otra) | 838 | |
| 12 · duplicado dentro de una fuente | 5 | **no estaba previsto**: lo encontró la validación |

## Validación

`04_validar.py` comprueba la base contra referencias **externas**, no contra sí misma:

- **VIRmiRNA**: 1.308 miRNAs de 44 virus, los mismos que la Tabla S1 del artículo original.
- **miRBase 23.1 viral**: 530 maduros y 320 horquillas.
- Integridad: ninguna molécula sin secuencia ni sin nombre, un único nombre principal
  por molécula, y la diferencia entre entradas y alias cuadra con los duplicados.

Un pipeline sin estas comprobaciones produce un resultado aunque la descarga se quede a
medias, y nadie se entera.

---

## Fuentes

| Fuente | Versión | Papel |
|---|---|---|
| [miRBase](https://mirbase.org) | 23.1 (sep 2026) | Referencia principal, humana y viral |
| miRBase | 22 (mar 2018) | Solo para la tabla de equivalencias de nombres |
| [VIRmiRNA](https://doi.org/10.1093/database/bau103) | 2014 | miRNAs virales validados experimentalmente |
| AntiVIRmiR / VIRmiRNA2 | 2022 | **Pendiente**: servidor caído |
| NCBI RefSeq viral, GRCh38, NCBI Taxonomy | | Pendientes (pasos A5–A7) |

Las versiones exactas están en [`config/fuentes.yml`](config/fuentes.yml), que es lo
único que hay que tocar para generar otra versión de la base.

**Sobre VIRmiRNA:** su servidor original lleva caído desde al menos noviembre de 2026.
No se parte de ninguna copia privada: se recupera de las páginas archivadas en Internet
Archive (que son URLs permanentes y citables) y se valida contra la tabla publicada en
el artículo. Está explicado en [`docs/02_fuentes.md`](docs/02_fuentes.md).

## Documentación

- [`docs/01_esquema_sql.md`](docs/01_esquema_sql.md) — cómo está montada la base y por qué
- [`docs/02_fuentes.md`](docs/02_fuentes.md) — de dónde sale cada dato y qué problemas tiene
- [`sql/esquema.sql`](sql/esquema.sql) — el esquema, comentado

## Principios FAIR

| | Cómo |
|---|---|
| **Localizable** | Repositorio público, cada versión con su etiqueta y (previsto) DOI en Zenodo |
| **Accesible** | Todo son archivos de texto y un SQLite; sin servidor, sin registro |
| **Interoperable** | FASTA, GFF3 y TSV; identificadores de miRBase, VMR_ID y taxid de NCBI |
| **Reutilizable** | Procedencia por archivo (`MANIFEST.tsv`), licencia MIT, reglas escritas en código y no aplicadas a mano |

## Licencia y citación

Código bajo licencia MIT (ver [`LICENSE`](LICENSE)). Los datos pertenecen a sus fuentes
originales y hay que citarlas: están en `config/fuentes.yml`, campo `cita`.

Para citar este repositorio, ver [`CITATION.cff`](CITATION.cff).
