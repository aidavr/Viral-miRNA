# De dónde sale cada dato

Documenta la procedencia de cada fuente, los problemas que tiene y cómo se han
resuelto. Es la parte "R" de FAIR (reutilizable): sin esto, los datos no se pueden
reinterpretar ni corregir.

---

## miRBase 23.1 (septiembre de 2026)

Referencia principal, tanto humana como viral. Se descarga de
`https://mirbase.org/download/CURRENT/`.

**Hallazgo importante:** la versión 23.1 **no añade ninguna secuencia viral nueva**
respecto a la 22. Comprobado descargando y comparando las dos: los mismos 530 maduros
virales de 34 especies. Lo que sí cambia son **84 nombres**, que pasan a llevar sufijo
de brazo:

```
ebv-miR-BART12   ->  ebv-miR-BART12-3p
hcmv-miR-UL148D  ->  hcmv-miR-UL148D-3p
```

Por eso actualizar la base de datos, por sí solo, no cambia los resultados del análisis
viral. Lo que los cambia es el mapeo competitivo contra el genoma humano.

### Tres cosas del formato que hay que saber

1. **El servidor sirve los `.fa` envueltos en HTML.** No son FASTA limpios: llevan
   `<br>` y entidades HTML. Hay que limpiarlos antes de parsearlos.

2. **Las cabeceras no tienen el mismo formato en los dos archivos:**

   ```
   mature.fa    >ebv-miR-BHRF1-1-5p MIMAT0000995 Epstein Barr virus miR-BHRF1-1-5p
   hairpin.fa   >ebv-mir-BHRF1-1    MI0001064    Epstein Barr virus miR-BHRF1-1 stem-loop
   ```

   El nombre de especie ocupa un número variable de palabras (`Homo sapiens` son 2,
   `Epstein Barr virus` son 3), así que no se puede cortar por posición: hay que quitar
   `stem-loop` si está, y después el nombre corto del miRNA.

3. **miRBase 22 trae los nombres de especie truncados.** Esto es de la fuente, no del
   parser:

   ```
   miRBase 22    >hsv1-miR-H1-3p MIMAT0015220 Herpes Simplex miR-H1-3p
   miRBase 23.1  >hsv1-miR-H1-3p MIMAT0015220 Herpes Simplex Virus 1 miR-H1-3p
   ```

   La base guarda lo que dice cada fuente, así que ambas formas conviven. Es un motivo
   más para resolver la especie por **taxid de NCBI** (paso A6) y no por el texto.

### Por qué también se descarga la 22

No aporta secuencias. Se usa solo para construir la tabla de equivalencias de nombres:
la matriz de miRDeep2 que ya tenemos lleva los nombres de la 22, y hay que poder
traducirla.

---

## VIRmiRNA (2014)

miRNAs virales validados experimentalmente: 1.308 miRNAs de 44 virus.

### El problema

El servidor original (`crdd.osdd.net`) **lleva caído desde al menos el 25/11/2026**
(conexión rechazada, no es un fallo temporal). Tampoco responde el servidor del grupo
que lo mantiene (`bioinfo.imtech.res.in`).

Existía una copia de estos datos en el pipeline anterior del grupo, **pero no se usa**:
no está publicada en ningún sitio, así que nadie de fuera podría reproducir el
resultado partiendo de ella. Un pipeline que depende de un archivo privado no es
reproducible, por muy bien escrito que esté.

### La solución: dos fuentes citables

**1. Las páginas archivadas en Internet Archive.** Guardan las tablas completas de cada
virus, con todas las columnas que hacen falta:

```
VMR_ID | miRNA | Virus | miRNA_Seq | Length | GC | ... | Método | Diana | PMID
```

Una URL de Wayback con marca de tiempo es **permanente e inmutable**:

```
https://web.archive.org/web/20170721045627id_/http://crdd.osdd.net/servers/virmirna/browse.php?by=HSV1&TYPE=Virus
                            └── fecha ──┘└┘
                                      el sufijo "id_" pide el contenido
                                      original, sin la barra de Wayback
```

El script pide el snapshot más cercano a una fecha y **anota en el MANIFEST el que
obtuvo**, así que la descarga se repite exactamente igual.

**2. La Tabla S1 del artículo**, que da el número de miRNAs de cada virus. Como el
artículo tiene DOI y el material suplementario está archivado en PMC, es una referencia
estable. Sirve para **comprobar que la recuperación está completa**:

> Qureshi A, Thakur N, Monga I, Thakur A, Kumar M. VIRmiRNA: a comprehensive resource
> for experimentally validated viral miRNAs and their targets. *Database* 2014;2014:bau103.
> [doi:10.1093/database/bau103](https://doi.org/10.1093/database/bau103)

La tabla está en `datos/referencia/virmirna_tabla_s1.tsv`, con los PMID de los trabajos
originales de cada virus. El resultado: **1.308 de 1.308, en los 44 virus**.

### Problemas de los datos

**Cinco secuencias con notación especial.** Se normalizan y se deja anotado de dónde
vienen; no se descartan:

| Ejemplo | Qué es | Qué se hace |
|---|---|---|
| `uuaguguuguccccccgagugg(c)` | base final opcional entre paréntesis | se quita (la forma larga ya aparece como isoforma, caso 7) |
| `cccggagcucgaaacccgguucgb` | `b` es un código IUPAC ambiguo (C, G o U) | pasa a `N`, que es lo que entienden bowtie y miRDeep2 |

**Cinco miRNAs duplicados.** El mismo nombre y la misma secuencia con **dos VMR_ID**:

```
prv-miR-10-1        VMR_0607  y  VMR_0608
prv-miR-8-3p        VMR_0789  y  VMR_0790
prv-miR-1-5p        VMR_0661  y  VMR_0662
rlcv-miR-rl1-34-3p  VMR_1037  y  VMR_1038
rlcv-miR-rl1-35-5p  VMR_1050  y  VMR_1051
```

Probablemente el mismo miRNA descrito en dos publicaciones. **Este caso no estaba
previsto en el diseño:** apareció porque la validación contra el artículo no cuadraba
(1.303 en vez de 1.308). Es exactamente para lo que sirve validar contra una referencia
externa.

**Nota sobre dos miRNAs de PRV.** La presentación del pipeline anterior decía que
`prv-miR-9-3p` y `prv-miR-5-p` se habían descartado a mano por estar mal anotados, pero
seguían en su tabla. Aquí no se descarta nada a mano: los dos entran, y la base registra
que tienen varias secuencias con el mismo nombre (caso 4). Quien los use, decide.

---

## AntiVIRmiR / VIRmiRNA2 (2022) — pendiente

Es la actualización de VIRmiRNA, con cuatro partes:

| Parte | Qué es | Papel aquí |
|---|---|---|
| **VIRmiRNA2** | ~3.340 miRNAs virales validados, de 50–56 virus | **secuencias** para la referencia |
| **AntiVmiR** | ~1.642 miRNAs del huésped con actividad antiviral o proviral, 34 virus | se cotejan como una fuente más, y aportan la marca antiviral/proviral |
| **DEmiRVIR** | expresión diferencial de miRNAs durante infecciones | anotación |
| **VIRmiRTar2** | dianas de los miRNAs virales | análisis de dianas |

**Estado: sin acceso.** El servidor no responde y el material suplementario del
artículo es un PDF sin los datos. Opciones pendientes: reintentar más adelante, buscar
copias archivadas (solo consta la portada), o escribir a los autores.

> Thakur A, Kumar M. AntiVIRmiR: a repository of host antiviral miRNAs and their
> expression along with experimentally validated viral miRNAs and their targets.
> *Front Genet* 2022;13:971852.
> [doi:10.3389/fgene.2022.971852](https://doi.org/10.3389/fgene.2022.971852)

**Importante:** cuando se consiga el acceso, sus miRNAs humanos se cotejan igual que
todo lo demás. No se da por hecho que ya estén en miRBase: pueden venir con nombres
antiguos, con otro formato o no estar.

---

## El universo: humano + virus

miRBase trae 323 especies. Cargarlas todas hace que los casos de cotejo dejen de
significar nada: los miRNAs están muy conservados entre animales (`let-7` es idéntico
en decenas de especies), así que "misma secuencia en varias especies" le pasa a casi
todo y tapa lo que sí interesa ver, que es la confusión **entre virus** y **entre virus
y humano**.

Con todas las especies salían 8.993 casos de "misma secuencia, varias especies"; con el
universo acotado, 575, y además son los interesantes.

El humano entra porque es el **competidor**: sin él no se puede saber si una lectura
atribuida a un virus es en realidad humana.

El filtro actual es provisional (por el texto del nombre de especie) y se sustituirá en
el paso A6 por el linaje de NCBI Taxonomy. Mientras tanto, `datos/universo_especies.tsv`
lista todas las especies aceptadas para poder revisarlas.

---

## Pendientes

| Fuente | Para qué | Paso |
|---|---|---|
| NCBI RefSeq viral 238 | recuperar horquillas, comprobar la anotación, predecir *novel* | A5, A6 |
| GRCh38.p14 / GENCODE 50 | competencia de mapeo, marca de similitud con lo humano | A7 |
| NCBI Taxonomy | especie, subfamilia y hospedador correctos | A6 |
| AntiVIRmiR | miRNAs virales validados y marca antiviral | — |
