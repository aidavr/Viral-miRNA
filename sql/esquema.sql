-- =====================================================================
-- VenomiR v2 — esquema de la base de datos
-- =====================================================================
-- Idea central: la UNIDAD de la base es la SECUENCIA, no el nombre.
-- Un mismo miRNA aparece en varias bases de datos con nombres distintos
-- (ebv-miR-BART1-3p = rlcv-miR-rl1-6-3p). Si la unidad fuese el nombre,
-- lo contaríamos dos veces. Por eso:
--     molecula  = una secuencia única   (la verdad biológica)
--     alias     = los nombres que le da cada fuente  (la etiqueta)
--
-- Convenciones del esquema:
--   * nombres de tabla y columna en minúscula y en español
--   * toda tabla derivada se puede reconstruir desde `entrada`
--   * nada se borra: lo dudoso se registra en `marca`
--
-- SQLite: se activan las claves foráneas en cada conexión con
--     PRAGMA foreign_keys = ON;
-- (SQLite las ignora por defecto, es un fallo clásico).
-- =====================================================================

PRAGMA foreign_keys = ON;


-- ---------------------------------------------------------------------
-- 1. fuente — de dónde viene cada dato (PROCEDENCIA)
-- ---------------------------------------------------------------------
-- Sin esta tabla la base no es reproducible: hay que poder decir
-- "esta secuencia viene de miRBase 23.1, descargada el 26/11/2026,
-- de esta URL, y el archivo tenía este sha256".
-- `sha256` es la huella del archivo: si alguien repite la descarga y
-- obtiene otra huella, el contenido ha cambiado.
CREATE TABLE IF NOT EXISTS fuente (
    id            INTEGER PRIMARY KEY,
    nombre        TEXT NOT NULL,            -- 'mirbase', 'virmirna', ...
    version       TEXT NOT NULL,            -- '23.1', '2014', 'wayback-20170721'
    fecha_fuente  TEXT,                     -- fecha de publicación de la fuente
    fecha_descarga TEXT NOT NULL,           -- cuándo la bajamos nosotros (ISO 8601)
    url           TEXT NOT NULL,
    sha256        TEXT,                     -- huella del archivo descargado
    prioridad     INTEGER NOT NULL,         -- 1 = manda su nombre. Ver `alias.es_principal`
    cita          TEXT,                     -- DOI o PMID del artículo que la describe
    UNIQUE (nombre, version)
);


-- ---------------------------------------------------------------------
-- 2. entrada — los datos TAL CUAL vienen, sin tocar
-- ---------------------------------------------------------------------
-- Es la capa "cruda normalizada": una fila por registro de cada fuente,
-- con las columnas homogeneizadas pero SIN decidir todavía nada.
-- Todo lo demás se deriva de aquí, así que si cambiamos una regla de
-- cotejo basta con volver a derivar, sin descargar otra vez.
CREATE TABLE IF NOT EXISTS entrada (
    id            INTEGER PRIMARY KEY,
    fuente_id     INTEGER NOT NULL REFERENCES fuente(id),
    id_fuente     TEXT,                     -- identificador propio de la fuente (MIMAT..., VMR_0289)
    nombre        TEXT NOT NULL,            -- nombre tal cual ('ebv-miR-BART1-3p')
    nombre_norm   TEXT NOT NULL,            -- nombre normalizado para comparar (ver docs)
    tipo          TEXT NOT NULL CHECK (tipo IN ('maduro','horquilla')),
    secuencia     TEXT,                     -- RNA en mayúsculas (A,C,G,U). NULL si la fuente no la da
    especie_texto TEXT,                     -- nombre de especie tal cual lo escribe la fuente
    extra         TEXT,                     -- JSON con lo que sea propio de esa fuente (PMID, método...)
    -- El identificador de la fuente forma parte de la clave. Si no se
    -- incluye, dos registros distintos de la misma base con el mismo
    -- nombre y secuencia se pisan y se pierde uno SIN AVISAR. Pasa de
    -- verdad: VIRmiRNA tiene 5 miRNAs repetidos con dos VMR_ID cada uno
    -- (p. ej. prv-miR-10-1 = VMR_0607 y VMR_0608). Lo detectó
    -- scripts/04_validar.py al no cuadrar con la Tabla S1 del artículo.
    UNIQUE (fuente_id, id_fuente, nombre, tipo, secuencia)
);

CREATE INDEX IF NOT EXISTS ix_entrada_sec    ON entrada(secuencia);
CREATE INDEX IF NOT EXISTS ix_entrada_nomnor ON entrada(nombre_norm);


-- ---------------------------------------------------------------------
-- 3. molecula — una fila por SECUENCIA única
-- ---------------------------------------------------------------------
-- Aquí es donde se resuelve el caso 1 (misma secuencia, distinto nombre):
-- todas las entradas con la misma secuencia colapsan en una molécula.
--   nivel: 'canonica'   maduro con horquilla conocida en la fuente
--          'derivada'   su horquilla la hemos reconstruido del genoma
--          'solo_maduro' no tiene horquilla: se contará por mapeo directo
CREATE TABLE IF NOT EXISTS molecula (
    id            INTEGER PRIMARY KEY,
    secuencia     TEXT NOT NULL UNIQUE,     -- la clave real de la base
    longitud      INTEGER NOT NULL,
    tipo          TEXT NOT NULL CHECK (tipo IN ('maduro','horquilla')),
    nivel         TEXT CHECK (nivel IN ('canonica','derivada','solo_maduro'))
);


-- ---------------------------------------------------------------------
-- 4. alias — todos los nombres de cada molécula
-- ---------------------------------------------------------------------
-- Una molécula tiene N alias (uno por fuente y versión que la nombra).
-- `es_principal = 1` marca el nombre que usaremos en los FASTA; se elige
-- por `fuente.prioridad` (miRBase manda).
-- Esta tabla es la que permite traducir la matriz antigua a la nueva:
-- es la tabla de equivalencias.
CREATE TABLE IF NOT EXISTS alias (
    id            INTEGER PRIMARY KEY,
    molecula_id   INTEGER NOT NULL REFERENCES molecula(id),
    fuente_id     INTEGER NOT NULL REFERENCES fuente(id),
    nombre        TEXT NOT NULL,
    nombre_norm   TEXT NOT NULL,
    id_fuente     TEXT,
    es_principal  INTEGER NOT NULL DEFAULT 0 CHECK (es_principal IN (0,1)),
    UNIQUE (molecula_id, fuente_id, nombre)
);

CREATE INDEX IF NOT EXISTS ix_alias_mol    ON alias(molecula_id);
CREATE INDEX IF NOT EXISTS ix_alias_nomnor ON alias(nombre_norm);


-- ---------------------------------------------------------------------
-- 5. taxon y molecula_taxon — a qué virus pertenece
-- ---------------------------------------------------------------------
-- Separadas porque una secuencia puede estar en VARIAS especies
-- (caso 2: jcv-miR-J1-3p = bkv-miR-B1-3p). Si la especie fuese una
-- columna de `molecula` no podríamos representarlo.
CREATE TABLE IF NOT EXISTS taxon (
    taxid         INTEGER PRIMARY KEY,      -- taxid de NCBI
    especie       TEXT NOT NULL,
    genero        TEXT,
    subfamilia    TEXT,                     -- Alpha/Beta/Gammaherpesvirinae...
    familia       TEXT,
    hospedador    TEXT CHECK (hospedador IN
                   ('humano','otro_vertebrado','invertebrado','planta','desconocido'))
);

CREATE TABLE IF NOT EXISTS molecula_taxon (
    molecula_id   INTEGER NOT NULL REFERENCES molecula(id),
    taxid         INTEGER NOT NULL REFERENCES taxon(taxid),
    asignacion    TEXT NOT NULL CHECK (asignacion IN ('unica','principal','ambigua')),
    PRIMARY KEY (molecula_id, taxid)
);


-- ---------------------------------------------------------------------
-- 6. madura_en_horquilla — dónde cae el maduro dentro del precursor
-- ---------------------------------------------------------------------
-- miRDeep2 y mirtop cuentan un maduro A TRAVÉS de su horquilla, así que
-- necesitan saber en qué posición empieza y acaba. De esta tabla sale
-- el GFF3 que se le pasa a mirtop (--mirna_gtf).
--   origen: 'fuente'   la relación la da la base de datos original
--           'derivada' la hemos calculado alineando contra el genoma
CREATE TABLE IF NOT EXISTS madura_en_horquilla (
    maduro_id     INTEGER NOT NULL REFERENCES molecula(id),
    horquilla_id  INTEGER NOT NULL REFERENCES molecula(id),
    inicio        INTEGER,                  -- 1-based, dentro de la horquilla
    fin           INTEGER,
    brazo         TEXT CHECK (brazo IN ('5p','3p','desconocido')),
    origen        TEXT NOT NULL CHECK (origen IN ('fuente','derivada')),
    PRIMARY KEY (maduro_id, horquilla_id)
);


-- ---------------------------------------------------------------------
-- 7. marca — todo lo dudoso, lo raro y lo que hay que vigilar
-- ---------------------------------------------------------------------
-- En vez de borrar o "arreglar" en silencio, se anota. Esta tabla es la
-- que genera el descriptivo de los casos de cotejo, y la que permite
-- filtrar después ("dame solo lo que no está marcado como humano").
-- Tipos previstos:
--   caso1_multinombre      misma secuencia con nombres distintos
--   caso2_multiespecie     misma secuencia en varias especies
--   caso3_multilocus       misma secuencia en varios loci del mismo virus
--   caso4_conflicto_nombre mismo nombre con secuencias distintas
--   caso5_renombrado       cambio de nombre entre versiones de una fuente
--   caso6_formato_nombre   solo cambia mayúsculas/guiones
--   caso7_isoforma         una secuencia contenida en otra
--   caso8_similar_humano   se parece a un miRNA humano (valor = E)
--   caso9_sin_horquilla    maduro sin precursor conocido
--   caso10_hospedador      virus de hospedador no humano
--   antiviral / proviral   anotación funcional (AntiVIRmiR)
CREATE TABLE IF NOT EXISTS marca (
    id            INTEGER PRIMARY KEY,
    molecula_id   INTEGER NOT NULL REFERENCES molecula(id),
    tipo          TEXT NOT NULL,
    detalle       TEXT,                     -- texto legible
    valor         REAL,                     -- número si lo hay (E-valor, nº de diferencias)
    ref_molecula_id INTEGER REFERENCES molecula(id)   -- la otra molécula implicada, si aplica
);

CREATE INDEX IF NOT EXISTS ix_marca_mol  ON marca(molecula_id);
CREATE INDEX IF NOT EXISTS ix_marca_tipo ON marca(tipo);


-- ---------------------------------------------------------------------
-- 8. version_base — qué versión de la base es esta
-- ---------------------------------------------------------------------
-- (No se llama `release` porque RELEASE es una palabra reservada de SQL
--  y obligaría a entrecomillarla en cada consulta.)
CREATE TABLE IF NOT EXISTS version_base (
    version       TEXT PRIMARY KEY,         -- '2026.12'
    fecha         TEXT NOT NULL,
    commit_git    TEXT,                     -- commit exacto del código que la generó
    notas         TEXT
);


-- =====================================================================
-- VISTAS — consultas frecuentes, guardadas con nombre
-- =====================================================================
-- Una vista no guarda datos: es una consulta con nombre. Útil para no
-- repetir los mismos JOIN y para que los scripts lean algo legible.

-- v_molecula: la ficha completa de cada molécula (nombre principal,
-- cuántos alias tiene, su especie y si está marcada como parecida a humano)
DROP VIEW IF EXISTS v_molecula;
CREATE VIEW v_molecula AS
SELECT
    m.id,
    m.secuencia,
    m.longitud,
    m.tipo,
    m.nivel,
    (SELECT a.nombre FROM alias a
      WHERE a.molecula_id = m.id AND a.es_principal = 1 LIMIT 1)        AS nombre_principal,
    (SELECT COUNT(*) FROM alias a WHERE a.molecula_id = m.id)           AS n_alias,
    (SELECT COUNT(DISTINCT a.fuente_id) FROM alias a
      WHERE a.molecula_id = m.id)                                       AS n_fuentes,
    (SELECT GROUP_CONCAT(t.especie, ' | ') FROM molecula_taxon mt
       JOIN taxon t ON t.taxid = mt.taxid
      WHERE mt.molecula_id = m.id)                                      AS especies,
    EXISTS (SELECT 1 FROM marca k
             WHERE k.molecula_id = m.id AND k.tipo = 'caso8_similar_humano')
                                                                        AS parecida_humano
FROM molecula m;

-- v_equivalencias: la tabla de traducción nombre antiguo -> nombre nuevo.
-- Es la que usaremos para releer la matriz de miRDeep2 vieja.
DROP VIEW IF EXISTS v_equivalencias;
CREATE VIEW v_equivalencias AS
SELECT
    a.nombre                                   AS nombre_alias,
    f.nombre || ' ' || f.version               AS fuente,
    p.nombre                                   AS nombre_principal,
    a.molecula_id
FROM alias a
JOIN fuente  f ON f.id = a.fuente_id
JOIN alias   p ON p.molecula_id = a.molecula_id AND p.es_principal = 1
WHERE a.es_principal = 0;

-- v_casos: recuento de cada caso de cotejo. Es el "descriptivo" que
-- se publica con cada release.
DROP VIEW IF EXISTS v_casos;
CREATE VIEW v_casos AS
SELECT tipo, COUNT(*) AS n_moleculas
FROM (SELECT DISTINCT tipo, molecula_id FROM marca)
GROUP BY tipo
ORDER BY tipo;
