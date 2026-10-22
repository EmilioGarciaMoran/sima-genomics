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
| 5 | `scripts/05_pairwise.py` | Compone la matriz sitios × especímenes desde los vectores y calcula distancias, alelos derivados compartidos y D (ver abajo); se omite si aún no hay vectores | `results/pairwise/` |
| 6 | `scripts/06_report.py` | Tablas resumen leídas de los TSV | `results/REPORT.md` |

Los parámetros están en `config.env` (ensamblaje, MAPQ, calidad de base, ventana, recorte de extremos, semilla 12345). Cualquiera se puede sobrescribir con una variable de entorno.

## Distancias por pares: nube → vectores → matriz → distancias

Recorrer genomas completos no cabe en un recorte, así que este análisis se parte en dos:

```text
NUBE (cloud/run_vectors.sh, una vez)                       REPOSITORIO (run_all.sh, siempre)
BAM / VCF / FASTA completos                                pairwise/panel.tsv.gz
  ├─ 09_site_panel.py  → panel de sitios informativos  ─►  pairwise/vectors/*.vec.gz  ─►  05_pairwise.py
  └─ 10_specimen_vector.py → un vector por espécimen                                         ├─ matriz
                                                                                             └─ distancias, D
```

```bash
# En la instancia (con los datos en DATA_DIR; ver pairwise/specimens.tsv)
BUCKET=s3://sima-egarmo-2026 bash cloud/run_vectors.sh DATA_DIR hg19.fa 8
# De vuelta en local: traer panel y vectores, ejecutar y versionar
aws s3 sync s3://sima-egarmo-2026/pairwise/ pairwise/ --exclude "logs/*"
bash run_all.sh && git add pairwise/ results/pairwise/
```

- **Panel.** Transversiones bialélicas variables en algún comparador de alta cobertura o en el chimpancé. Sima no interviene en la selección, así que sus errores no pueden crear sitios.
- **Vector.** Un alelo por sitio y por espécimen, también en los de alta cobertura (uno de los dos alelos al azar), para que todas las distancias se midan igual. Un alelo que no es REF ni ALT (error, daño, tercer alelo) queda como N. Cada vector registra el md5 del panel con el que se hizo, y la matriz rechaza vectores de otro panel.
- **Las distancias dependen del panel.** Un espécimen que interviene en la selección de sitios sale más lejos de todos. En el caso simulado, el gorila queda más cerca de Sima que el chimpancé porque el panel incluye los sitios donde el chimpancé difiere. Las comparaciones válidas son entre especímenes con el mismo papel en el panel (p. ej. Sima frente a Neandertal y frente a Denisova), no entre filas de simios.
- **Distancia p** con error estándar por jackknife de bloques de 5 Mb. Es relativa al panel, no una divergencia por pb. En un espécimen de baja cobertura incluye su tasa de error: la fila Sima–chimpancé frente a neandertal–chimpancé mide ese exceso.
- **D(Sima; X, Y)** cuenta los alelos derivados (polarizados con el chimpancé) que Sima comparte con X y no con Y, y al revés. Los errores de Sima caen por igual en los dos lados, así que es el estadístico robusto para decir de quién está más cerca.
- **Prueba.** `tests/test_pairwise.py` simula datos sobre un árbol conocido: Sima a ~1x con un 2% de error y daño, arcaicos en VCF y simios en FASTA. Comprueba que se recupera la topología y que D y los alelos derivados de la Sima de baja cobertura coinciden, dentro del error, con los de su genoma verdadero sin errores.

Panel de especímenes propuesto en `pairwise/specimens.tsv`: Sima (FemurXIII y los 4 BAM juntos), Altai, Vindija 33.19, Chagyrskaya 8, Denisova 3, Ust-Ishim, cinco humanos actuales del SGDP (Mbuti, yoruba, francés, han, papú) y chimpancé (grupo externo), gorila y orangután proyectados a hg19.

## Decisiones que conviene conocer

- **Daño post mortem.** Las posiciones de transición (C/T, G/A) se marcan, y toda llamada se repite descartando las 3 bases de cada extremo de lectura. En las GL se separan las transversiones, que son las menos afectadas por el daño.
- **Pseudohaploide reproducible.** La lectura elegida depende de una semilla fija, la muestra y la posición, no del orden de ejecución.
- **Sitios con base no-REF ≠ variantes.** A ~1x, `gl_summary.tsv` cuenta candidatos. La interpretación depende de cuántos sobreviven a los filtros de transversión y extremos.
- **Cromosoma X.** El modelo es diploide. Para genes ligados a X, como TKTL1, en un individuo masculino la llamada correcta es haploide.

## Estado

- [x] Flujo completo probado con datos sintéticos (`tests/`, 14 pruebas: detección de mezcla hg19/hg38 y de recortes alterados, recuperación de un árbol simulado, vectores deterministas y rechazo de vectores de otro panel)
- [ ] `catalog/positions.tsv` con las 66 posiciones y `catalog/cds_regions.bed` con las 189 regiones CDS
- [ ] Recortes reales en `slices/` y resultados versionados en `results/`
- [ ] Panel y vectores reales en `pairwise/` (ejecutar `cloud/run_vectors.sh` en la nube)
- [ ] Fijar versiones exactas del entorno (`micromamba env export`) tras la primera ejecución real
