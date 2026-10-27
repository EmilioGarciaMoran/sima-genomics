# Nuclear DNA: independent reproduction of Meyer et al. 2016
# ADN nuclear: reproducción independiente de Meyer et al. 2016

## Result / Resultado

**D(Denisova 3, Vindija 33.19; Sima, ancestral) = 0.586 ± 0.062** (weighted block jackknife,
563 blocks of 5 Mb), **Z = 9.5**, over 174 informative sites pooling the two best-preserved
specimens. Exact binomial p = 2.8 × 10⁻¹⁵. 90% CI [0.48, 0.69].

The two specimens, processed separately and sharing no sites, agree: **0.582** (incisor
AT-5482, n = 91) and **0.590** (femur fragment AT-5431, n = 83).

**This confirms Meyer et al. 2016: the Sima de los Huesos hominins fall on the Neanderthal
lineage.** The alternative hypothesis tested — that Sima is a sister lineage to Neanderthals
and Denisovans combined, which is what a separate taxon would require — is rejected. The
decision rule (D ≥ 0.50) was fixed in [`PREREGISTRO_D_statistics.md`](PREREGISTRO_D_statistics.md)
before any Sima read was mapped.

Meyer's per-branch fractions reproduce within error for the Neanderthal branch (0.46 and 0.33
against 0.43 and 0.39 published) and the Denisovan branch (0.08 and 0.11 against 0.09 and
0.07). The common branch reproduces for AT-5482 (0.64 against 0.68) but **not** for AT-5431
(0.70 against 0.87, about three standard deviations); this is recorded as an unresolved
discrepancy.

Full output: [`sima_pipeline/results/RESULTADO.txt`](sima_pipeline/results/RESULTADO.txt).

## Limitations / Limitaciones

1. **The ENA deposit is pre-mapped.** Between 96% and 99.9% of the deposited reads align to
   the human reference, against the 1–10% expected of raw ancient-DNA libraries, and
   duplicates are already removed. PRJEB10597 therefore holds the subset of reads that
   Meyer's own pipeline accepted, roughly 60% of what he analysed. **The remapping here is
   not independent of his:** the allele counts and the new D-statistics are reproducible, the
   original mapping decisions are not verifiable. A full reproduction would need the raw
   data, which is not deposited.
2. **174 informative sites suffice to reject the sister-lineage hypothesis, not to sustain
   it.** With this sample size the pre-registered equivalence test (δ = 0.10) could only have
   supported it had D come out essentially zero. See §9bis of the pre-registration for the
   sample sizes that stricter margins would require: about 8 M usable bases for δ = 0.05,
   against the 3.2 M available.
3. **The statistic measures lineage membership, not taxonomic distance.** It places Sima on
   the Neanderthal branch; how far along that branch is a separate question, and the
   specimen-specific sets (about 200 sites each) have little resolution for it.

## Estructura

| Carpeta | Contenido |
|---|---|
| [`PREREGISTRO_D_statistics.md`](PREREGISTRO_D_statistics.md) | Hipótesis, umbrales, controles obligatorios, reglas de interpretación y parámetros del preprocesado. **Fijado antes de los datos**; las correcciones posteriores llevan fecha y motivo. |
| [`power_sim/`](power_sim/) | Estudio de potencia por simulación coalescente. Qué preguntas tienen respuesta con los datos de Sima y cuáles no. |
| [`branch_counts/`](branch_counts/) | Construcción de los conjuntos de referencia a partir de los genomas de alta cobertura: 1.822.345 sitios informativos sobre los 22 autosomas. Incluye la validación del software. |
| [`sima_pipeline/`](sima_pipeline/) | Mapeo, filtro de daño, extracción de alelos y análisis de los cinco especímenes de Sima. |

## Reproducir

El orden importa: los conjuntos de referencia son la entrada del pipeline de Sima.

```bash
# 1. Validar el software de recuento por rama (no necesita datos de Sima)
cd branch_counts
python3 scripts/test_branch_counts.py
python3 scripts/test_read_vcfs.py
python3 scripts/sim_check.py

# 2. Construir los conjuntos de referencia (~8 h, ~200 GB de tránsito, 21 GB de disco)
bash scripts/run_genome.sh
python3 scripts/check_genome.py results/genome
python3 scripts/merge_genome.py results/genome --out results/genoma_inf.npz

# 3. Preparar y ejecutar el pipeline de Sima (~2 h la preparación, ~30 min el resto)
cd ../sima_pipeline
bash scripts/prepare.sh
bash scripts/map_sima.sh
bash scripts/damage_sima.sh
bash scripts/extract_sima.sh
python3 scripts/analyze_sima.py --sitios ../branch_counts/results/genoma_inf.npz \
        --alelos results/alelos --out results/RESULTADO.txt
```

Cada paso es reanudable y verifica sus entradas. Las pruebas de respuesta conocida
(`scripts/test_*.py`, unas 120 en total) se ejecutan en segundos y no necesitan datos.

## Cómo se trabajó

El preregistro se cerró antes de tocar los datos de Sima y las desviaciones posteriores
quedaron documentadas con su fecha y su motivo, incluidas tres correcciones hechas al leer el
artículo original: la calidad de mapeo (30, no 25), el criterio de daño y el enmascarado.

Cada pieza lleva pruebas de respuesta conocida escritas antes de ejecutarla con datos reales.
Encontraron seis fallos, entre ellos uno silencioso que descartaba el primer registro de cada
VCF y otro en el que la propia prueba enmascaraba un error de rutas del módulo.

Las cifras nucleares del preprint de 2026 (D = 0,7425, Z = 30,6 sobre 218.757.031 sitios)
quedaron retiradas al no encontrarse ninguna salida que las respaldara. El estudio de
potencia muestra además que un Z de ese orden es inalcanzable con los volúmenes de datos de
Sima.
