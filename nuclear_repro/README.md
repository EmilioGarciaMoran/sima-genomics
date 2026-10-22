# ADN nuclear reproducible — Sima de los Huesos

El flujo del análisis nuclear se ejecuta de principio a fin sobre **BAM recortados** a las regiones que el análisis consulta: las posiciones del catálogo ± 2 kb y las CDS de los genes de interés. Los recortes ocupan pocos MB y se versionan en `slices/`, así que cualquiera puede reproducir los resultados en minutos, sin AWS y sin descargar genomas completos. Las coordenadas y las cabeceras `@SQ` no cambian: un recorte se analiza exactamente igual que el BAM completo.

Los análisis:

1. **Genotipado directo** de las posiciones funcionales del catálogo, con llamada pseudohaploide.
2. **Verosimilitudes genotípicas** sitio a sitio sobre las CDS fusionadas.
3. **Distancias por pares** entre Sima y otros homínidos, humanos antiguos y actuales, y simios, a partir de un vector de alelos por espécimen generado en la nube (ver abajo).

## Cómo reproducirlo

```bash
cd nuclear_repro
micromamba create -f environment.yml      # python, pysam, numpy, samtools, angsd (opcional)
micromamba activate nucrepro
bash run_all.sh                           # usa slices/ y catalog/; escribe results/
git diff --stat results/                  # si no cambia nada, se ha reproducido
```

Prueba sin datos reales (caso sintético con resultados conocidos; es lo que ejecuta GitHub Actions):

```bash
python3 -m pytest -v tests
```

## Generar los recortes (una sola vez, donde están los BAM completos)

```bash
# 1. Catálogo: rellenar catalog/positions.tsv y generar las CDS
python3 scripts/make_cds_bed.py --gtf Homo_sapiens.GRCh37.87.gtf.gz --ucsc-chroms \
    --genes FOXP2 CNTNAP2 SCN9A SLC24A5 MC1R SLC45A2 TKTL1 NOVA1 ASPM TLR1 OAS1 ADSL BAZ1B ... \
    > catalog/cds_regions.bed
# 2. Recortar BAM y referencia
python3 scripts/00_make_slices.py --ref /ruta/hg19.fa --bam /ruta/Sima_*.hg19.bam
# 3. Ejecutar y versionar
bash run_all.sh
git add slices/ catalog/ results/
```

`slices/MANIFEST.tsv` guarda la ruta y el md5 de cada BAM completo de origen, el número de lecturas conservadas y el md5 de cada recorte.

## Pasos

| Paso | Script | Qué hace | Salidas |
|---|---|---|---|
| 0 | `scripts/00_make_slices.py` | Recorta BAM y referencia a las regiones del catálogo (se ejecuta una vez, fuera del flujo) | `slices/` |
| 1 | `scripts/01_check_inputs.py` | **Aborta** si un recorte no coincide con su md5, si un BAM está en otro ensamblaje (longitud de chr1) o si el REF de alguna posición no coincide con la referencia, que es el síntoma de mezclar hg19 y hg38 | — |
| 2 | `scripts/02_direct_genotype.py` | Por muestra y en conjunto (POOL): recuentos A/C/G/T, llamada pseudohaploide ANC/DER, la misma sin los extremos de lectura y lecturas a ± 2 kb | `results/direct_genotypes.tsv`, `results/direct_summary.tsv` |
| 3 | `scripts/03_genotype_likelihoods.py` | GL de cada sitio CDS cubierto (modelo de GATK, el de ANGSD `-GL 2`), en Python legible y probado | `results/gl_sites.tsv`, `results/gl_summary.tsv` |
| 3b | `scripts/03b_angsd.sh` | Lo mismo con ANGSD y los parámetros exactos del manuscrito (`-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -minMapQ 20 -minQ 20`), también con `-noTrans 1`; se omite si ANGSD no está instalado | `results/angsd/` |
| 5 | `scripts/05_pairwise.py` | Compone la matriz sitios × especímenes desde los vectores de todas las ventanas procesadas y calcula distancias, alelos derivados compartidos y D (ver abajo); se omite si aún no hay ventanas procesadas | `results/pairwise/` |
| 6 | `scripts/06_report.py` | Tablas resumen leídas de los TSV | `results/REPORT.md` |

Los parámetros están en `config.env` (ensamblaje, MAPQ, calidad de base, ventana, recorte de extremos, semilla 12345). Cualquiera se puede sobrescribir con una variable de entorno.

## Distancias por pares: ventanas → nube → vectores → matriz → distancias

En lugar de recorrer genomas completos, el análisis trabaja por **ventanas de coordenadas fijas** (p. ej. 100.000–200.000 de cada cromosoma). Cada ventana es una unidad independiente: se procesa por separado en la nube, sus ficheros no cambian al añadir otras, y en el jackknife es un bloque. Así el análisis escala y se acumula: se empieza con pocas ventanas y se van añadiendo ventanas o especímenes sin rehacer nada.

```text
pairwise/windows.bed  (cloud/make_windows.py, se amplía cuando se quiera)
        │
NUBE (cloud/run_vectors.sh, solo lo que falta)              REPOSITORIO (run_all.sh, siempre)
  por ventana:                                              pairwise/windows/<id>/panel.tsv.gz
    09_site_panel.py --region  → panel de esa ventana   ─►  pairwise/windows/<id>/vectors/*.vec.gz
    10_specimen_vector.py      → un vector por espécimen           │
                                                                    └─► 05_pairwise.py: matriz (ventanas
                                                                        concatenadas) → distancias, D
```

```bash
# 1. Ventanas: una por cromosoma, 100.000-200.000
python3 cloud/make_windows.py --fai hg19.fa.fai --ref hg19.fa --start 100000 --size 100000
# 2. En la instancia (ver "Primera ejecución en la nube")
BUCKET=s3://sima-egarmo-2026 bash cloud/run_vectors.sh DATA_DIR hg19.fa 8
# 3. De vuelta en local: traer, ejecutar y versionar
aws s3 sync s3://sima-egarmo-2026/pairwise/windows/ pairwise/windows/
bash run_all.sh && git add pairwise/ results/pairwise/

# Más adelante: ampliar (se suman a las anteriores) y volver a lanzar 2 y 3
python3 cloud/make_windows.py --fai hg19.fa.fai --ref hg19.fa --start 5000000 --size 100000 --every 10000000
```

### Primera ejecución en la nube (lanzada a mano)

La instancia se arranca a mano. El script se encarga de que, una vez lanzado, el coste esté acotado y no se pierda nada:

1. **Instancia bajo demanda, no spot**, con el rol IAM de S3 (sin claves en la máquina). Como red de seguridad independiente del script, conviene una alerta de AWS Budgets.
2. **Plan, sin calcular nada:** `PLAN=1 bash cloud/run_vectors.sh DATA_DIR hg19.fa` muestra cuántas ventanas y trabajos hay pendientes y qué límites se aplicarán.
3. **Piloto de una ventana:** `MAX_WINDOWS=1 BUCKET=s3://... bash cloud/run_vectors.sh DATA_DIR hg19.fa 8`. Revisa `pairwise/runs/<fecha>/SUMMARY.txt` y los tiempos de `done`. Con eso sabes cuánto tarda cada ventana antes de lanzar el resto.
4. **El resto,** con límites ajustados a lo visto en el piloto: `MAX_HOURS=3 BUCKET=... bash cloud/run_vectors.sh DATA_DIR hg19.fa 8`. Añade `SHUTDOWN=1` si no vas a estar delante cuando termine.

Lo que garantiza el script:

- **Antes de calcular:** comprueba que existen todos los ficheros de entrada y que se puede escribir en S3. Si algo falta, se detiene sin gastar.
- **Una sola ejecución a la vez** sobre el mismo directorio: una segunda se rechaza (cerrojo).
- **Fallos aislados:** cada trabajo tiene tiempo máximo (`JOB_TIMEOUT`, 1 h) y un reintento. Si vuelve a fallar se apunta en `runs/<fecha>/failed`, los demás siguen, y la siguiente ejecución lo reintenta.
- **Hora límite** (`MAX_HOURS`, 4 h por defecto): pasado ese tiempo no se empieza nada nuevo, y lo que queda pendiente se recoge en la siguiente ejecución.
- **Pase lo que pase** (fallos, hora límite, Ctrl-C): al salir escribe el resumen y el `MD5SUMS`, sube a S3 todo lo terminado (con 3 intentos) y, si `SHUTDOWN=1`, apaga la instancia. El código de salida es 1 si hubo fallos.
- **Latido** cada 5 minutos en `s3://.../pairwise/runs/<fecha>/HEARTBEAT.txt`, para ver el progreso desde fuera.

`tests/test_cloud_run.py` prueba cada una de estas garantías con un `aws` simulado: plan sin efectos, piloto, entrada corrupta, hora límite, fichero ausente, doble Ctrl-C y ejecución simultánea.

Como solo se leen las ventanas, `path` en `specimens.tsv` puede ser una URL a un BAM o VCF indexado (por ejemplo, los VCF públicos de Leipzig): htslib descarga solo los tramos necesarios, sin bajar el genoma entero.

- **Sesgo de las ventanas, y cómo diluirlo.** Las posiciones 100.000–200.000 están junto al telómero. Son regiones subteloméricas, con más GC, más duplicaciones segmentarias, peor mapabilidad y tasas de mutación y recombinación atípicas, que pueden inflar las distancias y el error de Sima. En hg19, los cromosomas acrocéntricos (13, 14, 15, 21, 22) son solo N en ese tramo; `--ref` descarta esas ventanas y lo avisa. Al ampliar con ventanas repartidas (`--every`), el sesgo se diluye. `windows_summary.tsv` deja comparar resultados por ventana.
- **Unidad de acumulación.** Una ventana procesada no se rehace aunque se añada un comparador: el panel de cada ventana queda fijo y el nuevo comparador solo aporta su vector. Las escrituras son atómicas, así que un proceso cortado no deja ficheros que parezcan terminados. La matriz usa las ventanas procesadas en el orden de `windows.bed`, y un espécimen sin vector en una ventana cuenta como sin dato en ella.
- **Panel.** Transversiones bialélicas variables en algún comparador de alta cobertura o en el chimpancé. Sima no interviene en la selección, así que sus errores no pueden crear sitios.
- **Vector.** Un alelo por sitio y por espécimen, también en los de alta cobertura (uno de los dos alelos al azar), para que todas las distancias se midan igual. Un alelo que no es REF ni ALT (error, daño, tercer alelo) queda como N. Cada vector registra el md5 del panel con el que se hizo, y la matriz rechaza vectores de otro panel.
- **Las distancias dependen del panel.** Un espécimen que interviene en la selección de sitios sale más lejos de todos. En el caso simulado, el gorila queda más cerca de Sima que el chimpancé porque el panel incluye los sitios donde el chimpancé difiere. Las comparaciones válidas son entre especímenes con el mismo papel en el panel (p. ej. Sima frente a Neandertal y frente a Denisova), no entre filas de simios.
- **Distancia p** con error estándar por jackknife, un bloque por ventana. Es relativa al panel, no una divergencia por pb. En un espécimen de baja cobertura incluye su tasa de error: la fila Sima–chimpancé frente a neandertal–chimpancé mide ese exceso.
- **D(Sima; X, Y)** cuenta los alelos derivados (polarizados con el chimpancé) que Sima comparte con X y no con Y, y al revés. Los errores de Sima caen por igual en los dos lados, así que es el estadístico robusto para decir de quién está más cerca.
- **Prueba.** `tests/test_pairwise.py` comprueba la acumulación: procesar chr1, añadir chr2 y volver a lanzar deja chr1 intacto byte a byte y da resultados idénticos a hacerlo todo de una vez. También simula datos sobre un árbol conocido: Sima a ~1x con un 2% de error y daño, arcaicos en VCF y simios en FASTA. Comprueba que se recupera la topología y que D y los alelos derivados de la Sima de baja cobertura coinciden, dentro del error, con los de su genoma verdadero sin errores.

Panel de especímenes propuesto en `pairwise/specimens.tsv`: Sima (FemurXIII y los 4 BAM juntos), Altai, Vindija 33.19, Chagyrskaya 8, Denisova 3, Ust-Ishim, cinco humanos actuales del SGDP (Mbuti, yoruba, francés, han, papú) y chimpancé (grupo externo), gorila y orangután proyectados a hg19.

## Decisiones que conviene conocer

- **Daño post mortem.** Las posiciones de transición (C/T, G/A) se marcan, y toda llamada se repite descartando las 3 bases de cada extremo de lectura. En las GL se separan las transversiones, que son las menos afectadas por el daño.
- **Pseudohaploide reproducible.** La lectura elegida depende de una semilla fija, la muestra y la posición, no del orden de ejecución.
- **Sitios con base no-REF ≠ variantes.** A ~1x, `gl_summary.tsv` cuenta candidatos. La interpretación depende de cuántos sobreviven a los filtros de transversión y extremos.
- **Cromosoma X.** El modelo es diploide. Para genes ligados a X, como TKTL1, en un individuo masculino la llamada correcta es haploide.

## Estado

- [x] Flujo completo probado con datos sintéticos (`tests/`, 24 pruebas: detección de mezcla hg19/hg38 y de recortes alterados, recuperación de un árbol simulado, vectores deterministas, acumulación de ventanas, rechazo de vectores de otro panel y modos de fallo de la ejecución en la nube)
- [ ] `catalog/positions.tsv` con las 66 posiciones y `catalog/cds_regions.bed` con las 189 regiones CDS
- [ ] Recortes reales en `slices/` y resultados versionados en `results/`
- [ ] Primeras ventanas reales (100.000–200.000 de cada cromosoma) procesadas en la nube y versionadas en `pairwise/windows/`
- [ ] Fijar versiones exactas del entorno (`micromamba env export`) tras la primera ejecución real
