# Pipeline de Sima de los Huesos — ADN nuclear

## RESULTADO (26 de octubre de 2026)

**D(Denisova 3, Vindija 33.19; Sima, ancestral) = 0,586 ± 0,062 (jackknife), Z = 9,5**
sobre 174 sitios informativos, combinando AT-5431 y AT-5482. Prueba binomial exacta:
p = 2,8 × 10⁻¹⁵. IC 90%: [0,48, 0,69].

Los dos especímenes principales, procesados por separado y sin sitios en común,
coinciden: **0,582** (incisivo AT-5482, n = 91) y **0,590** (fémur AT-5431, n = 83).

**Aplicando la regla preregistrada (D ≥ 0,50): confirma Meyer et al. 2016. Sima de los
Huesos pertenece al linaje neandertal. La hipótesis de un linaje hermano de neandertales y
denisovanos —la que haría defendible una taxonomía propia— queda descartada por la vía
genética.**

La D predicha a partir de los porcentajes publicados por Meyer era 0,68 con IC [0,63, 0,75];
la observada es 0,586 con IC [0,48, 0,69]. Los intervalos se solapan.

### Controles

| Control | Resultado |
|---|---|
| Reproducción de Meyer, rama neandertal (f_N) | 0,46 y 0,33 frente a 0,43 y 0,39 publicados |
| Reproducción de Meyer, rama denisovana (f_D) | 0,08 y 0,11 frente a 0,09 y 0,07 publicados |
| Reproducción de Meyer, rama común (f_ND) | 0,70 y 0,64 frente a 0,87 y 0,68 publicados |
| Sin filtro de daño | D baja a 0,36 y 0,40, **en la dirección predicha** (más contaminación moderna diluye la señal arcaica) |
| Solo criterio de cadena sencilla | D sube a 0,69 y 0,64 |
| Solo criterio de doble cadena | D sube a 0,68 y 0,79 |
| Dispersión entre 20 muestreos pseudohaploides | 0,0000 y 0,0011, frente a un error típico de 0,087 y 0,089: despreciable |

La D es positiva y significativa en **todas** las variantes de sensibilidad, con magnitudes
entre 0,36 y 0,79. El resultado no depende del ajuste elegido.

### Lo que no reproduce bien

El único número que no cuadra dentro del error es **f_ND en AT-5431: 0,70 frente al 0,87
publicado**, a unas tres desviaciones típicas con n = 71. Las otras cinco comparaciones con
Meyer caen dentro del intervalo. Posibles causas: el subconjunto de lecturas depositado en
ENA, o que Meyer definiera la rama común con criterios algo distintos. **Debe constar como
discrepancia no resuelta.**

### Limitaciones

- ENA contiene en torno al 60% de las lecturas que Meyer analizó, y están premapeadas: el
  remapeo no es independiente del suyo (ver más abajo).
- 174 sitios informativos bastan para rechazar H2 con holgura, pero no para sostenerla: la
  prueba de equivalencia con esta n solo habría podido confirmar H2 si la D hubiera salido
  prácticamente nula.
- Mide pertenencia al linaje, no distancia taxonómica.



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
| 4. Mapeo con `bwa aln` sin siembra | **EJECUTADO** |
| 5. Filtro de daño (C→T terminal, enmascarado de 3+3) | **EJECUTADO** |
| 6. Extracción de alelos en los sitios de referencia | **EJECUTADO** |
| 7. Recuento por rama y D-statistics | **EJECUTADO** (`results/RESULTADO.txt`) |

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

## Pieza 2: filtro de daño y enmascarado

```bash
python3 scripts/test_damage_filter.py      # 30 pruebas de respuesta conocida
bash scripts/damage_sima.sh 2>&1 | tee results/dano.log
```

**Corrección al preregistro (26 oct 2026).** Meyer 2016 usó librerías de **cadena sencilla**
(protocolo de Gansauge y Meyer), que es lo que permitió recuperar ADN tan degradado. En esas
librerías la desaminación aparece como C→T en **ambos** extremos de la cadena original, no
como C→T en el 5' y G→A en el 3' como ocurre en las de doble cadena. Al alinear, una lectura
en hebra inversa muestra ese mismo daño como G→A, porque el SAM guarda el complemento
inverso.

Esto explica la frase del artículo «T enmascaradas en las 3 primeras y 3 últimas posiciones»,
que está escrita en la orientación de la cadena original. El preregistro decía «A en las 3
últimas», que es la convención de doble cadena y es incorrecta para estos datos.

La implementación es por tanto **dependiente de la hebra**, leída del FLAG. Se conserva la
convención de doble cadena como variante (`--library ds`) y se ejecuta siempre, de modo que
la desviación queda documentada con su comparación. **Pendiente de verificar en el artículo.**

Otras decisiones de implementación, tomadas aquí y probadas:

- Las posiciones se cuentan sobre la porción **alineada**, descartando el recorte blando: las
  bases recortadas no se alinean y no contribuyen a ninguna llamada de alelo.
- Se enmascara la base **resultante** del daño (T en hebra directa, A en inversa), no
  cualquier base, para no perder información sin motivo.

Salidas en `~/sima/sima_data/bam_dmg/`: el principal, más `.nofilt` y `.ds` de los análisis
de sensibilidad preregistrados.

**Control de validación:** las bases útiles resultantes deben acercarse a las de Meyer 2016
(Tabla 1). Si difieren mucho, el filtro está mal implementado.

## Resultado de la pieza 2 (26 oct 2026, ejecución final con MAPQ 30)

| Espécimen | Lecturas con daño | Tasa | Bases útiles | Meyer | Cociente |
|---|---|---|---|---|---|
| AT-5482 incisivo | 33.409 | 4,6% | 1.290.935 | 2.015.167 | 0,64 |
| AT-5431 fémur | 18.132 | 7,8% | 696.912 | 1.186.442 | 0,59 |
| Fémur XIII | 3.396 | 2,2% | 136.780 | 225.329 | 0,61 |
| Molar AT-5444 | 2.925 | 4,2% | 114.055 | 188.974 | 0,60 |
| Escápula AT-6672 | 652 | 0,7% | 40.799 | 37.524 | 1,09 |

**1.987.847 bases útiles entre los dos especímenes principales**, el 60% de los 3.201.609 de
Meyer.

### Comparación de criterios (ejecución anterior con MAPQ 25, cocientes prácticamente iguales)

Criterio de Meyer (C→T o G→A en la primera o última posición alineada), enmascarado de todas
las timinas en 3+3. Cociente entre nuestras bases útiles y las de Meyer 2016, Tabla 1:

| Criterio | AT-5482 | AT-5431 | Fémur XIII | Molar | Escápula | Media endógenos |
|---|---|---|---|---|---|---|
| **any (principal)** | 0,64 | 0,59 | 0,61 | 0,60 | 1,09 | **0,61 ± 0,02** |
| ss (sensibilidad) | 0,58 | 0,54 | 0,48 | 0,52 | 0,42 | 0,53 ± 0,04 |
| ds (sensibilidad) | 0,36 | 0,34 | 0,34 | 0,35 | 0,74 | 0,35 ± 0,01 |

**Interpretación.** El criterio principal da un cociente de 0,61 con desviación de 0,02 en
los cuatro especímenes endógenos. Esa constancia descarta un error de implementación, que
produciría dispersión entre especímenes, y apunta a un factor común: **ENA contiene en torno
al 60% de las lecturas que Meyer analizó**, coherente con que el depósito esté premapeado.

El artículo respalda este criterio para el análisis nuclear: dice que para el mitocondrial
usaron el filtro dependiente de hebra, pero que aquí usaron el enmascarado de timinas porque
«this filter is less strict than that based on alignment orientation described above and thus
retains more data».

La escápula se desvía (1,09) porque con un 98% de contaminación el criterio permisivo recoge
falsos positivos de daño en ADN moderno. Refuerza su exclusión, ya decidida.

**Consecuencia para el análisis.** Una pérdida uniforme y no selectiva del 40% reduce el
número de sitios, y con él la potencia, pero **no desplaza el valor esperado de la D**, que es
un cociente. Los sitios informativos esperados bajan de los 436 proyectados a unos 260
combinando AT-5431 y AT-5482. Sigue bastando para rechazar D = 0 si la D real ronda 0,68,
pero **la prueba de equivalencia del apartado 9 bis del preregistro queda más débil**: con 260
sitios el margen mínimo alcanzable sube de 0,079 a unos 0,10, de modo que solo se podría
sostener H2 con una D observada prácticamente nula.

**Control de contaminación (criterio principal).** Tasa de daño: escápula 0,7%, fémur XIII
2,2%, molar 4,2%, incisivo 4,6%, AT-5431 7,8%. La escápula tiene once veces menos daño que
AT-5431, como predice su contaminación.

## Reproducibilidad del muestreo pseudohaploide

La semilla está fijada (`--seed 20261025`), se pasa por argumento y se guarda dentro del
`npz`. Una prueba comprueba que dos ejecuciones con la misma semilla coinciden.

**Pero fijar la semilla da reproducibilidad, no estabilidad.** Garantiza que otro obtenga el
mismo número, no que ese número sea robusto. Si con una semilla saliera D = 0,68 y con otra
0,61, el resultado sería reproducible y frágil a la vez, que es peor que irreproducible
porque aparenta solidez.

Por eso el extractor guarda además **20 muestreos independientes** (`--repeticiones`), de una
sola pasada por el BAM. La repetición 0 es la principal; las demás permiten medir cuánto se
mueve la D por el azar del muestreo, que es una fuente de variación **distinta** de la que
mide el jackknife (variación entre regiones del genoma).

El muestreo solo interviene en los sitios cubiertos por más de una lectura. Con cobertura muy
por debajo de 1x se espera que sean pocos, y el extractor informa de cuántos son y de cuántos
cambian de valor entre repeticiones. **Si la dispersión entre muestreos resultara comparable
al error típico del jackknife, el muestreo estaría mandando sobre los datos y habría que
decirlo en cualquier resultado.**

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
