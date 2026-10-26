# Pipeline de Sima de los Huesos — ADN nuclear

Reproducción del recuento por rama de Meyer et al. 2016 con los datos de ENA PRJEB10597,
y D-statistics sobre los conjuntos de referencia construidos en `../branch_counts/`.

**Preregistro:** `../PREREGISTRO_D_statistics.md`. Hipótesis, umbrales, controles y reglas de
interpretación quedaron fijados antes de tocar estos datos.

## Estado

| Paso | Estado |
|---|---|
| 1. Metadatos de ENA: qué espécimen es cada corrida | **hecho** (`ESPECIMENES.tsv`) |
| 2. Referencia GRCh37 e índice de bwa | `prepare.sh` |
| 3. Descarga de los FASTQ | `prepare.sh` |
| 4. Mapeo con `bwa aln` sin siembra | **escrito y probado** (`map_sima.sh`) |
| 5. Filtro de daño (C→T terminal, enmascarado de 3+3) | pendiente |
| 6. Extracción de alelos en los sitios de referencia | pendiente |
| 7. Recuento por rama y D-statistics | pendiente |

## Preparación

```bash
bash scripts/prepare.sh                 # a ~/sima/sima_data
nohup bash scripts/prepare.sh > results/preparacion.log 2>&1 &
```

Unas 2 h, dominadas por `bwa index` (~90 min a un hilo). Disco ~9 GB. Reanudable: salta
cada paso cuyo resultado ya exista. Requiere `bwa`, `samtools`, `wget` y `curl`.

**No contiene ninguna decisión analítica**, solo descargas e indexado, por lo que puede
dejarse corriendo sin supervisión. Los pasos 4 a 7 son otra cosa y se escriben con sus
pruebas de respuesta conocida, como el resto del proyecto.

## Especímenes (resuelto el 25 oct 2026)

Del filereport de ENA para PRJEB10597. Nueve corridas agrupadas en cinco muestras:

| Espécimen | Muestra | Corridas | Lecturas | Bases útiles | Papel |
|---|---|---|---|---|---|
| Incisivo AT-5482 | SAMEA3518421 | ERR995358, ERR995362 | 727.796 | 2.015.167 | principal |
| Fémur AT-5431 | SAMEA3518420 | ERR995361, ERR995364 | 237.745 | 1.186.442 | principal |
| Fémur XIII | SAMEA3518423 | ERR995357 | 155.538 | 225.329 | confirmación |
| Molar AT-5444 | SAMEA3518424 | ERR995360, ERR995365 | 72.445 | 188.974 | descriptivo |
| Escápula AT-6672 | SAMEA3518422 | ERR995359, ERR995363 | 101.217 | 37.524 | **excluida** |

Bases útiles según Meyer et al. 2016, Tabla 1 (alineadas tras el filtro de daño).

**Las corridas de una misma muestra sí se fusionan**: son librerías del mismo hueso. Las
muestras distintas no se fusionan nunca. Confundir las dos cosas fue uno de los errores de
2026, y en el BAM fusionado de entonces entró la escápula.

La escápula confirma de paso que los metadatos son coherentes con el artículo: 101.217
lecturas y solo 37.524 bases útiles, justo lo que cabe esperar de una muestra contaminada al
98%. El incisivo es el caso opuesto.

El fémur XIII es el del mtDNA, y es el de menos datos nucleares de los tres utilizables: de
ahí que quede como confirmación y con prueba binomial exacta.

## Pieza 1: mapeo

```bash
bash scripts/test_filter_dedup.py          # 20 pruebas de respuesta conocida
bash scripts/map_sima.sh                   # ~15 min
nohup bash scripts/map_sima.sh > results/mapeo.log 2>&1 &
```

Un BAM por espécimen en `~/sima/sima_data/bam/`. Reanudable. Parámetros del apartado 10 ter
del preregistro: `bwa aln -l 16500 -n 0.01 -o 2`, MAPQ ≥ 25, longitud ≥ 35 pb.

`filter_dedup.py` hace el filtrado y la eliminación de duplicados en la propia tubería, sin
BAM intermedios. Elimina duplicados por **(cromosoma, inicio, fin, hebra)**: el fin se calcula
sobre la referencia a partir del CIGAR, no sumando la longitud de la secuencia, porque las
inserciones y deleciones hacen que ambas cosas no coincidan. Agrupar solo por inicio, que es
el error fácil, borraría fragmentos legítimos de distinta longitud que empiezan en el mismo
punto.

La escápula AT-6672 se mapea igualmente, para poder documentar su contaminación, pero queda
fuera del análisis.

## Hallazgo del mapeo (26 oct 2026): los datos de ENA están premapeados

Resultado de la pieza 1 (`results/mapeo.log`):

| Espécimen | Lecturas ENA | Conservadas | % | Duplicadas |
|---|---|---|---|---|
| AT-5482 incisivo | 727.796 | 722.569 | 99,3% | 315 (0,04%) |
| AT-5431 fémur | 237.745 | 233.793 | 98,3% | 219 (0,09%) |
| Fémur XIII | 155.538 | 155.393 | 99,9% | 20 (0,01%) |
| AT-5444 molar | 72.445 | 69.606 | 96,1% | 161 (0,22%) |
| AT-6672 escápula | 101.217 | 98.605 | 97,4% | 65 (0,06%) |

**Entre el 96 y el 99,9% de las lecturas mapean.** En ADN antiguo auténtico lo normal es
entre el 1 y el 10%: la mayoría de lo secuenciado es ADN microbiano del sedimento. Y los
duplicados son del orden del 0,05%, inverosímil en una librería amplificada.

**Conclusión: los FASTQ depositados en ENA no son los datos crudos, sino el subconjunto de
lecturas que ya mapeaban contra hg19, con los duplicados ya eliminados.** Encaja con el
volumen: Meyer secuenció entre 600 y 2.600 millones de lecturas por espécimen, y aquí hay
entre 70.000 y 730.000.

**Limitación que debe constar en cualquier texto que salga de aquí.** El remapeo no es
independiente: estamos realineando lecturas que el pipeline de Meyer ya había aceptado. Se
puede reproducir el recuento de alelos derivados por rama y calcular D-statistics nuevas,
pero **no** se puede verificar el mapeo original ni descartar un sesgo introducido en él. Es
una reproducción del análisis, no del pipeline completo. Presentarla de otro modo sería
inexacto.

Una reproducción completa exigiría los datos crudos, que no están depositados. Esa es otra
petición concreta que se puede hacer a Leipzig.

## Decisiones ya fijadas

- **Referencia `human_g1k_v37`**, cuyos cromosomas se llaman `1`, `2`, `21`, sin prefijo
  `chr`, igual que los VCF de MPI-EVA. Mezclar convenciones arruinaría el cruce de
  posiciones. Meyer 2016 usó hg19, que es el mismo ensamblaje.
- **Cada espécimen por separado.** Fusionar los BAM de las cinco corridas fue uno de los
  errores del trabajo de 2026, sobre todo porque una de ellas es la escápula, contaminada
  al 98%. De ahí que el paso 1 sea identificar qué corrida corresponde a cada espécimen.
- **Análisis principal: AT-5431 y AT-5482, por separado y combinados.** Fémur XIII solo como
  confirmación, con prueba binomial exacta por tener pocos sitios. Escápula excluida.
- **Conjuntos de referencia en modo homocigoto** desde los genomas de alta cobertura;
  muestreo pseudohaploide solo para Sima, donde es obligatorio.
- **Error típico por jackknife**, nunca binomial.

## Sitios esperados

De `../branch_counts/README.md`, apartado 3 bis:

| Espécimen | Rama profunda | Exclusivos D5 | Exclusivos Vi |
|---|---|---|---|
| Fémur XIII | ~126 | ~14 | ~17 |
| AT-5431 | ~664 | ~71 | ~90 |
| AT-5482 | ~1.128 | ~121 | ~153 |
| AT-5431 + AT-5482 | ~1.797 | ~192 | ~244 |
