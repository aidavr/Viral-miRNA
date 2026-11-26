# Cómo está montada la base, y por qué

Este documento explica el diseño de la base de datos. Está escrito para
entenderlo desde cero, no solo para consultarlo.

---

## 1. ¿Por qué una base de datos y no un Excel?

Podríamos tener una tabla con una fila por miRNA y columnas `nombre`,
`secuencia`, `fuente`, `especie`. El problema aparece en cuanto los datos son
reales:

- Un mismo miRNA está en **tres fuentes**, con **tres nombres**. ¿Una fila o tres?
  Si pones tres, cuentas el mismo miRNA tres veces. Si pones una, pierdes dos nombres.
- Una secuencia pertenece a **dos virus** (HSV-1 y HSV-2 comparten miRNAs). ¿Pones
  `"HSV-1, HSV-2"` en la celda? Entonces ya no puedes filtrar por virus sin partir texto.
- Un maduro está dentro de **varias horquillas**.

Todo esto son relaciones **de uno a muchos** o **de muchos a muchos**. Una tabla plana
las representa mal, y acabas metiendo listas dentro de celdas. El modelo relacional
existe justo para esto: **una tabla por tipo de cosa, y tablas que las relacionan**.

---

## 2. El concepto clave: clave primaria y clave foránea

```sql
CREATE TABLE molecula (
    id         INTEGER PRIMARY KEY,   -- clave primaria: identifica la fila
    secuencia  TEXT NOT NULL UNIQUE,
    ...
);

CREATE TABLE alias (
    id           INTEGER PRIMARY KEY,
    molecula_id  INTEGER NOT NULL REFERENCES molecula(id),   -- clave foránea
    nombre       TEXT NOT NULL,
    ...
);
```

- **Clave primaria** (`PRIMARY KEY`): identifica una fila sin ambigüedad.
- **Clave foránea** (`REFERENCES`): dice que `alias.molecula_id` apunta a una fila de
  `molecula`. La base impide guardar un alias que apunte a una molécula inexistente.

Así, una molécula puede tener muchos alias, y cada alias sabe de quién es.

> **Aviso práctico de SQLite:** las claves foráneas están **desactivadas por defecto**.
> Hay que pedirlo en cada conexión:
> ```sql
> PRAGMA foreign_keys = ON;
> ```
> Si no, la base acepta datos incoherentes sin protestar. Es el fallo más típico al
> empezar con SQLite.

---

## 3. Las tablas, por capas

```
  CAPA 1 · lo que dijo cada fuente (no se toca)
  ┌────────────┐        ┌─────────────┐
  │  fuente    │◄───────│  entrada    │   una fila por registro descargado
  └────────────┘        └─────────────┘
                              │  (se agrupa por secuencia)
  CAPA 2 · lo que decidimos   ▼
  ┌────────────┐        ┌─────────────┐        ┌───────────┐
  │   alias    │───────►│  molecula   │◄───────│   marca   │
  └────────────┘        └─────────────┘        └───────────┘
                           │        │
                           ▼        ▼
               molecula_taxon    madura_en_horquilla
                     │
                     ▼
                  taxon
```

**Capa 1 (`fuente`, `entrada`) es sagrada:** es el registro fiel de lo que publicó cada
base de datos. No se corrige ni se limpia.

**Capa 2 se borra y se recalcula entera** cada vez que se ejecuta `03_cotejar.py`. Esto
tiene una consecuencia muy buena: cambiar una regla de cotejo cuesta 2 segundos y no
obliga a descargar nada. Y garantiza que el resultado siempre es consecuencia de reglas
escritas en el código, nunca de un arreglo manual que nadie recuerda.

| Tabla | Una fila es... |
|---|---|
| `fuente` | una versión concreta de una base de datos (miRBase 23.1) |
| `entrada` | un registro tal cual lo publicó esa fuente |
| `molecula` | **una secuencia única** |
| `alias` | un nombre que una fuente le da a una molécula |
| `taxon` | una especie (con su taxid de NCBI) |
| `molecula_taxon` | "esta secuencia pertenece a esta especie" |
| `madura_en_horquilla` | "este maduro está en esta posición de esta horquilla" |
| `marca` | algo digno de anotar sobre una molécula |

---

## 4. El cotejo, en tres consultas

Esto es el corazón del pipeline. Está en `scripts/03_cotejar.py`.

### 4.1 Crear una molécula por secuencia única

```sql
INSERT INTO molecula (secuencia, longitud, tipo)
SELECT secuencia, LENGTH(secuencia), tipo
FROM entrada
WHERE secuencia IS NOT NULL
GROUP BY secuencia, tipo;
```

`GROUP BY secuencia, tipo` es **todo el cotejo**. Agrupa las entradas que comparten
secuencia, vengan de donde vengan y se llamen como se llamen. Si `ebv-miR-BART1-3p` de
miRBase y `rlcv-miR-rl1-6-3p` de VIRmiRNA tienen la misma secuencia, caen en el mismo
grupo y producen **una sola** molécula.

Se agrupa también por `tipo` porque un maduro y una horquilla son cosas distintas
aunque (en teoría) coincidiera la secuencia.

### 4.2 Enganchar cada nombre a su molécula

```sql
INSERT INTO alias (molecula_id, fuente_id, nombre, nombre_norm, id_fuente)
SELECT m.id, e.fuente_id, e.nombre, e.nombre_norm, e.id_fuente
FROM entrada e
JOIN molecula m ON m.secuencia = e.secuencia AND m.tipo = e.tipo;
```

`JOIN ... ON` empareja filas de dos tablas según una condición. Aquí: para cada entrada,
busca la molécula con su misma secuencia y tipo. Así cada nombre queda colgado de su
molécula.

### 4.3 Elegir el nombre principal

```sql
UPDATE alias SET es_principal = 1
WHERE id IN (
    SELECT a.id FROM alias a
    WHERE a.id = (
        SELECT a2.id FROM alias a2
        JOIN fuente f2 ON f2.id = a2.fuente_id
        WHERE a2.molecula_id = a.molecula_id
        ORDER BY f2.prioridad, a2.nombre
        LIMIT 1)
);
```

Para cada molécula se ordenan sus alias por la **prioridad de la fuente** (definida en
`config/fuentes.yml`: miRBase manda) y se marca el primero. `LIMIT 1` se queda con uno.

El segundo criterio (`a2.nombre`, orden alfabético) no tiene sentido biológico: está
para que, si hay empate, el resultado sea **siempre el mismo**. Un pipeline que da un
resultado distinto en cada ejecución no es reproducible.

---

## 5. Detectar los casos: `GROUP BY` + `HAVING`

Casi todos los casos de cotejo son la misma consulta con otra condición:

```sql
SELECT a.molecula_id
FROM alias a
GROUP BY a.molecula_id
HAVING COUNT(DISTINCT a.nombre_norm) > 1;
```

- `GROUP BY a.molecula_id` → un grupo por molécula.
- `COUNT(DISTINCT a.nombre_norm)` → cuántos nombres diferentes tiene.
- `HAVING ... > 1` → quedarse solo con las que tienen más de uno.

> **`WHERE` o `HAVING`?** `WHERE` filtra **filas antes** de agrupar; `HAVING` filtra
> **grupos después**. Para "moléculas con más de un nombre" hace falta contar primero,
> así que es `HAVING`.

El caso contrario (mismo nombre, secuencias distintas) es la misma idea al revés:

```sql
SELECT nombre_norm, COUNT(DISTINCT molecula_id) AS n
FROM alias
GROUP BY nombre_norm
HAVING n > 1;
```

### Lo que SQL no hace bien

El **caso 7** (una secuencia contenida en otra, es decir isoformas) se podría escribir
con `LIKE '%' || corta || '%'`, pero obligaría a comparar todas las parejas de
secuencias, sin poder usar índices. Con ~7.000 moléculas son ~24 millones de
comparaciones de texto. Se hace en Python, ordenando por longitud. **Regla general: SQL
para agrupar y relacionar; un lenguaje de programación para comparar contenidos.**

---

## 6. Vistas: consultas con nombre

Una vista no guarda datos; es una consulta guardada que se usa como si fuera una tabla:

```sql
CREATE VIEW v_equivalencias AS
SELECT a.nombre AS nombre_alias,
       p.nombre AS nombre_principal
FROM alias a
JOIN alias p ON p.molecula_id = a.molecula_id AND p.es_principal = 1
WHERE a.es_principal = 0;
```

Fíjate en que `alias` aparece **dos veces** con alias distintos (`a` y `p`): es un
*self join*, unir una tabla consigo misma. Aquí sirve para poner en la misma fila un
nombre secundario y el principal de su molécula.

Y ya se puede consultar como una tabla:

```sql
SELECT * FROM v_equivalencias WHERE nombre_alias = 'hcmv-miR-UL148D';
-- hcmv-miR-UL148D  ->  hcmv-miR-UL148D-3p
```

Esta vista es la que traducirá la matriz de miRDeep2 antigua a los nombres nuevos.

---

## 7. Cómo mirar la base por tu cuenta

Desde la línea de órdenes (si tienes instalado `sqlite3`):

```bash
sqlite3 datos/venomir_v2.sqlite
```

```sql
.tables                     -- qué tablas hay
.schema molecula            -- cómo es una tabla
.headers on
.mode column                -- salida legible

-- las moléculas con más nombres distintos
SELECT nombre_principal, n_alias, especies
FROM v_molecula
ORDER BY n_alias DESC
LIMIT 10;

-- todos los nombres de un miRNA concreto
SELECT f.nombre AS fuente, f.version, a.nombre, a.es_principal
FROM alias a
JOIN fuente f ON f.id = a.fuente_id
WHERE a.molecula_id = (SELECT molecula_id FROM alias
                       WHERE nombre = 'ebv-miR-BART1-3p' LIMIT 1);

-- los duplicados que encontró la validación
SELECT detalle FROM marca WHERE tipo = 'caso12_duplicado_fuente';
```

Desde Python:

```python
import sqlite3
con = sqlite3.connect("datos/venomir_v2.sqlite")
con.row_factory = sqlite3.Row          # para acceder por nombre de columna
for r in con.execute("SELECT * FROM v_casos"):
    print(r["tipo"], r["n_moleculas"])
```

Desde R:

```r
library(RSQLite)
con <- dbConnect(SQLite(), "datos/venomir_v2.sqlite")
dbGetQuery(con, "SELECT * FROM v_casos")
```

Para verla con ratón: [DB Browser for SQLite](https://sqlitebrowser.org/), gratuito.

---

## 8. Por qué SQLite y no PostgreSQL

| | SQLite | PostgreSQL |
|---|---|---|
| Qué es | **Un archivo** | Un servidor en marcha |
| Instalación | Ninguna (viene con Python) | Instalar, configurar, usuarios |
| Llevarlo al HPC | `scp base.sqlite` | No se puede "copiar" |
| Varios escribiendo a la vez | Mal | Bien |
| Tamaño razonable | Hasta decenas de GB | Lo que haga falta |

Aquí solo escribe un script, la base ocupa pocos MB y tiene que viajar a un clúster de
cálculo. SQLite es la opción correcta; PostgreSQL sería complicarse sin ganar nada.

**Y un matiz importante:** la base SQLite es la **fuente**, no el producto. Los pipelines
de nf-core no leen bases de datos: reciben archivos (`--mature`, `--hairpin`,
`--mirna_gtf`). Por eso el paso A9 exportará FASTA, GFF3 y TSV desde el SQLite. La base
es donde se razona; los archivos son lo que se entrega.
