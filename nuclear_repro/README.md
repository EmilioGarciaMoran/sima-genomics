# ADN nuclear reproducible — Sima de los Huesos

El flujo del análisis nuclear se ejecuta de principio a fin sobre **BAM recortados** a las regiones que el análisis consulta: las posiciones del catálogo ± 2 kb y las CDS de los genes de interés. Los recortes ocupan pocos MB y se versionan en `slices/`, así que cualquiera puede reproducir los resultados en minutos, sin AWS y sin descargar genomas completos. Las coordenadas y las cabeceras `@SQ` no cambian: un recorte se analiza exactamente igual que el BAM completo.

Los dos análisis del manuscrito:

1. **Genotipado directo** de las posiciones funcionales del catálogo, con llamada pseudohaploide.
2. **Verosimilitudes genotípicas** sitio a sitio sobre las CDS fusionadas.

## Cómo reproducirlo

```bash
cd nuclear_repro
micromamba create -f environment.yml      # python, pysam, samtools, angsd (opcional)
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
| 4 | `scripts/04_report.py` | Tablas resumen leídas de los TSV | `results/REPORT.md` |

Los parámetros están en `config.env` (ensamblaje, MAPQ, calidad de base, ventana, recorte de extremos, semilla 12345). Cualquiera se puede sobrescribir con una variable de entorno.

## Decisiones que conviene conocer

- **Daño post mortem.** Las posiciones de transición (C/T, G/A) se marcan, y toda llamada se repite descartando las 3 bases de cada extremo de lectura. En las GL se separan las transversiones, que son las menos afectadas por el daño.
- **Pseudohaploide reproducible.** La lectura elegida depende de una semilla fija, la muestra y la posición, no del orden de ejecución.
- **Sitios con base no-REF ≠ variantes.** A ~1x, `gl_summary.tsv` cuenta candidatos. La interpretación depende de cuántos sobreviven a los filtros de transversión y extremos.
- **Cromosoma X.** El modelo es diploide. Para genes ligados a X, como TKTL1, en un individuo masculino la llamada correcta es haploide.

## Estado

- [x] Flujo completo probado con datos sintéticos (`tests/`, 7 pruebas, incluida la detección de mezcla hg19/hg38 y de recortes alterados)
- [ ] `catalog/positions.tsv` con las 66 posiciones y `catalog/cds_regions.bed` con las 189 regiones CDS
- [ ] Recortes reales en `slices/` y resultados versionados en `results/`
- [ ] Fijar versiones exactas del entorno (`micromamba env export`) tras la primera ejecución real
