# Registro de cambios

Formato: [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/).

## [Sin publicar]

### Añadido
- Esquema de la base (`sql/esquema.sql`): moléculas, alias, taxonomía, marcas.
- A1 descarga versionada con `MANIFEST.tsv` (URL, fecha, sha256).
- A2 normalización de miRBase (FASTA) y VIRmiRNA (tabla HTML archivada).
- A3-A4 cotejo por secuencia y detección de los casos 1, 2, 4, 5, 6, 7 y 12.
- A8 validación contra la Tabla S1 del artículo de VIRmiRNA y contra miRBase.

### Encontrado al construirlo
- VIRmiRNA tiene **5 miRNAs duplicados** con dos identificadores cada uno
  (`prv-miR-10-1` = VMR_0607 y VMR_0608, y otros cuatro). Es el caso 12, que no
  estaba previsto en el diseño: lo detectó la validación al no cuadrar con el
  artículo (1.303 frente a 1.308).
- miRBase sirve los archivos `.fa` envueltos en HTML; hay que limpiarlos.
- Las cabeceras de `hairpin.fa` acaban en `stem-loop`, las de `mature.fa` no:
  si se recorta igual en los dos, el nombre del miRNA acaba dentro de la especie.
- miRBase **22** trae los nombres de especie truncados (`Herpes Simplex` en vez de
  `Herpes Simplex Virus 1`); la 23.1 ya los corrige. Es un motivo más para
  resolver la especie por taxid de NCBI y no por el texto de la cabecera.
- 5 secuencias de VIRmiRNA usan notación especial: una base final opcional entre
  paréntesis y códigos IUPAC ambiguos (`b`). Se normalizan y se deja anotado.

### Pendiente
- A5 horquillas desde los genomas virales (RNAfold).
- A6 taxonomía con NCBI Taxonomy (sustituye al filtro por texto de especie).
- A7 marca de similitud con miRNAs y genoma humanos.
- A9 exportación a FASTA, GFF3 y genoma combinado.
- Acceso a AntiVIRmiR / VIRmiRNA2 (servidor caído).
