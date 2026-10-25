# Ficheros necesarios para la validación del recuento por rama

Verificado en cdna.eva.mpg.de el 25 de octubre de 2026.

## Resumen de la decisión

**Validamos con el cromosoma 21, no con el genoma completo.** Con unos 100.000 sitios el error
de muestreo es de ±0,2 puntos, quince veces menor que nuestra tolerancia de ±3 puntos
(`README.md`). Bajar los 23 cromosomas multiplicaría por 40 el volumen sin mejorar la decisión.

**No usamos los ficheros de ENA.** Son lecturas crudas; habría que remapear y llamar genotipos.
Los VCF de cdna.eva.mpg.de ya están procesados y, salvo Denisova 17, todos sobre hg19.

**Denisova 17 queda fuera.** Ver el apartado 3.

## 1. Ficheros a descargar (cromosoma 21, hg19)

Base: `https://ftp.eva.mpg.de/neandertal/`
(`https://cdna.eva.mpg.de/neandertal/` es el mismo servidor; el script prueba las dos.)

| Genoma | Ruta relativa a la base | Tamaño |
|---|---|---|
| Denisova 5 (Altai) | `Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz` | 634 M |
| Vindija 33.19 | `Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz` | 622 M |
| Denisova 3 | `Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz` | 596 M |
| Chagyrskaya 8 | `Chagyrskaya/VCF/chr21.noRB.vcf.gz` | 576 M |
| Mbuti (HGDP00982) | `altai/ModernHumans/vcf/SS6004471.hg19_1000g.21.mod.vcf.gz` | — |

De cada uno, **también el `.tbi`** del mismo nombre. Y el `md5sums.txt` de cada directorio,
para comprobar las descargas.

Total aproximado: **2,5 GB**.

### Descarga

```bash
bash scripts/fetch_chr21.sh            # a data/chr21/
bash scripts/fetch_chr21.sh /ruta/que/sea
```

El script comprueba primero que las 10 URL existen, y si alguna falla no descarga nada: así un
cambio en el servidor se detecta antes de bajar gigabytes. Es reanudable (`wget -c`), elige
host automáticamente y renombra los ficheros al espécimen, porque todos los directorios usan
el mismo nombre `chr21_mq25_mapab100.vcf.gz`.

El directorio `Vindija/VCF/` contiene además `LBK/`, `Loschbour/`, `Mez1/`, `Ust_Ishim/` e
`indels/`, que no hacen falta. Ust'-Ishim en concreto **no** se usa: fue el grupo externo
equivocado del preprint retirado.

El identificador de Mbuti está confirmado en `altai/ModernHumans/vcf/README_Bteam`:
SS6004471 = Mbuti = HGDP00982, que es exactamente el individuo que cita Massilani.

### Estado ancestral

Massilani define el ancestral con chimpancé, bonobo, gorila, orangután y macaco rhesus, a
partir de los alineamientos del navegador UCSC. La alternativa práctica es la secuencia
ancestral de Ensembl para GRCh37 (`homo_sapiens_ancestor_GRCh37`, alineamiento EPO de seis
primates), que viene ya resuelta por posición y para el cromosoma 21 ocupa unos pocos MB.

**No son el mismo conjunto de primates.** Esa diferencia es justamente la clase de
discrepancia de implementación que cubre la tolerancia de ±3 puntos. Hay que documentar cuál
se usó. Comprobar el nombre exacto del fichero en el FTP de Ensembl antes de bajarlo.

## 2. Qué validamos ahora que Denisova 17 no está disponible

Las cifras 91/45/20 están confirmadas literalmente en el artículo (apartado *Sexing and
Lineage Assignment*): «The Denisova 17 genome carries 91% of the derived alleles found in both
Neandertals, as well as 45% and 20% of those unique to Denisova 5 and Vindija 33.19».

Pero esas tres cifras exigen el genoma de D17, que no se puede usar (apartado 3). La
validación se traslada a lo que sí es reproducible con los VCF de arriba:

**a) Prueba de topología con Chagyrskaya 8.** Aplicar los tres estadísticos con
X = Chagyrskaya 8, A = Denisova 5, B = Vindija 33.19. La topología establecida (Mafessoni et
al. 2020; Fig. 1 de Massilani) exige **f_B > f_A**: Chagyrskaya está más emparentada con
Vindija que con el neandertal del Altái. Es la misma comprobación que ya pasa en simulación
(`scripts/sim_check.py`), ahora con datos reales.

**b) Reproducir los paneles D, E y F de la Fig. 1 de Massilani** en las barras que no
involucran a D17. Son D-statistics de la forma D(ind1, ind2; ind3, Mbuti) con ind3 = D5,
Chag8 y Vi33.19. El criterio es reproducir el **signo y la significación** (|Z| ≥ 3) de cada
barra, no su valor exacto.

**c) Coherencia del conjunto de referencia.** El número de sitios de cada conjunto (AB, A, B)
debe ser estable al cambiar de cromosoma, en proporción a su longitud.

Si las tres pasan, el pipeline está validado para pasar a Sima.

## 3. Por qué Denisova 17 queda fuera

`cdna.eva.mpg.de/neandertal/Denisova17/` contiene **únicamente** `rawBAM/`, con un fichero:

```
A26088_merged.T2T_CHM13v2.0_evan.bam    166 G
```

Dos problemas, y el segundo es el grave:

1. **No hay VCF.** Habría que llamar genotipos sobre 166 GB de lecturas.
2. **Está mapeado contra T2T-CHM13v2.0**, no contra hg19. Todos los demás genomas, y todo
   Meyer 2016, están en hg19.

Para usarlo habría que remapear las lecturas a hg19 o convertir coordenadas entre ensamblajes.
Cualquiera de las dos cosas **cambia el resultado respecto al de Massilani**, que trabajó con
su propio mapeo. Es decir: el esfuerzo no produciría una comparación limpia con las cifras que
queremos validar. Se estaría comparando contra un blanco que nosotros mismos hemos movido.

**Conclusión: no se remapea D17.** Si en algún momento Leipzig publica el VCF sobre hg19, o si
la conversación con Kelso y Prüfer prospera, la comprobación se hace entonces. Ellos ya tienen
esos datos procesados.

## 4. Accesiones citadas por Massilani (para referencia)

| Genoma | Accesión |
|---|---|
| Denisova 17 | ENA PRJEB108287 (estudio ERP189149, público 2026-03-09) |
| Goyet 1 (GN1) | ENA PRJEB98484 |
| Denisova 25 | Peyrégne et al.; comprobar si ya está publicado |

Ninguno de los tres hace falta para la validación descrita arriba.
