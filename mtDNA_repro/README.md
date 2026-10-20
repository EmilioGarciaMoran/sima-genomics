# mtDNA reproducible — Sima de los Huesos

Re-ejecución completa y documentada de la filogenia mitocondrial (6 genomas, *Pan troglodytes* como grupo externo), desde las secuencias hasta los tests de topología. Cada paso es un script y cada resultado está versionado en el commit que lo produjo.

## Cómo reproducirlo

```bash
cd mtDNA_repro
micromamba create -f environment.yml      # mafft 7.526, iqtree 3.0.1, raxml-ng 1.2.2
micromamba activate mtphylo
bash run_all.sh genbank 2                 # ~16 min con 2 núcleos
git diff --stat                           # si no cambia nada en results/, se ha reproducido
```

Semilla fija (12345) en todos los análisis estocásticos.

## Pasos

| Paso | Script | Qué hace | Salidas |
|---|---|---|---|
| 1 | `scripts/01_get_sequences.sh` | Obtiene las 6 secuencias de `accessions.tsv` (GenBank, o sin huecos del alineamiento legacy) | `sequences/` |
| 2 | `scripts/02_align.sh` | MAFFT G-INS-i; versión sin D-loop (posiciones 600–16000 de Vindija); comparación con el alineamiento legacy | `alignments/` |
| 3 | `scripts/03_trees.sh` | Árboles ML enraizados en Pan: IQ-TREE (modelo BIC + UFBoot 1000; HKY+G4 + bootstrap estándar 1000) y RAxML-NG (GTR+G + bootstrap 1000). Recuento de agrupaciones en las réplicas | `results/trees/`, `results/CLADE_SUPPORT.tsv`, `results/ML_TREES.tsv` |
| 4 | `scripts/04_topology_tests.sh` | H0/H1/H2 evaluadas juntas: KH, SH, AU (10.000 RELL) | `results/topology/`, `results/TOPOLOGY_TESTS.tsv` |
| 5 | `scripts/05_legacy_check.sh` | Relee los resultados originales de la fase 3a con el mismo recuento | `results/LEGACY_PHASE3A_CHECK.txt` |

Hipótesis: **H0** Sima basal a (Denisova + Neandertales) · **H1** Sima hermano de Denisova · **H2** Sima hermano de los Neandertales.

## Resultados

**Alineamiento.** 17.150 columnas y 292 sitios informativos, igual que el alineamiento legacy de la fase 3a; 17.138 de 17.150 columnas emparejan los mismos residuos. Las 12 diferencias son huecos desplazados una posición en el D-loop. Sin D-loop: 15.409 columnas, 256 sitios informativos.

**Soporte bootstrap** (% de réplicas que contienen cada agrupación; `results/CLADE_SUPPORT.tsv`):

| Alineamiento | Análisis | Sima+Denisova | Sima+Neandertales | Denisova+Neandertales |
|---|---|---|---|---|
| completo | IQ-TREE TPM2u+F+G4, UFBoot | 98.0 | 2.0 | 0.0 |
| completo | IQ-TREE HKY+G4, bootstrap estándar | 92.9 | 7.1 | 0.0 |
| completo | RAxML-NG GTR+G, bootstrap estándar | 97.9 | 2.1 | 0.0 |
| sin D-loop | IQ-TREE HKY+F+G4, UFBoot | 97.4 | 1.2 | 1.4 |
| sin D-loop | IQ-TREE HKY+G4, bootstrap estándar | 98.8 | 0.4 | 0.8 |
| sin D-loop | RAxML-NG GTR+G, bootstrap estándar | 99.6 | 0.3 | 0.1 |

**Tests de topología** (p-AU; `results/TOPOLOGY_TESTS.tsv`):

| Alineamiento | Modelo | H0 Sima basal | H1 Sima+Denisova | H2 Sima+Neandertales |
|---|---|---|---|---|
| completo | HKY+G4 | 0.011 | 0.946 | 0.071 |
| completo | TPM2u+F+G4 | 0.009 | 0.962 | 0.049 |
| sin D-loop | HKY+G4 | 0.038 | 0.963 | 0.038 |
| sin D-loop | HKY+F+G4 | 0.038 | 0.963 | 0.038 |

**Conclusión de los datos:** con enraizamiento explícito en *Pan*, el mtDNA de Sima es hermano del de Denisova en todos los análisis. H0 (Sima basal) queda rechazada por el test AU en los cuatro casos. Esto coincide con Meyer et al. (2014), que usaron un análisis bayesiano (MrBayes) enraizado en chimpancé y bonobo.

## Relación con la fase 3a (marzo 2026)

Los cálculos de la fase 3a eran correctos; el error fue de lectura (`results/LEGACY_PHASE3A_CHECK.txt`):

- IQ-TREE se ejecutó sin `-o Pan_outgro`; el propio `.iqtree` avisa de que el árbol no tiene raíz y Sima solo aparece dibujada en ella.
- El soporte 93 de ese árbol corresponde a la bipartición {Sima, Denisova} | {Neandertales, Pan}, es decir, al clado Sima+Denisova.
- Las 1.000 réplicas de RAxML-NG de la fase 3a apoyan Sima+Denisova en un 98,4% y nunca recuperan Denisova+Neandertales.
- Los "KH tests" de la fase 3a evaluaban cada hipótesis por separado (un árbol por ejecución), por eso nunca dieron p-valores.
- `inputMitocondria/remove_dloop.py` no quitaba el D-loop (eliminaba las primeras 1.122 columnas); marcado como obsoleto.

## Limitaciones

- Solo 6 taxones, sin humanos modernos ni más denisovanos (Meyer 2014 usó ~120 genomas).
- Las secuencias versionadas proceden del modo `recover` porque NCBI no era accesible desde el entorno de construcción. Sus longitudes coinciden con GenBank. Para cerrar la verificación, hay que ejecutar `bash run_all.sh genbank` en local y comprobar que `sequences/SEQUENCES.tsv` no cambia.

## Verificación de reproducibilidad (20 oct 2026)

Ejecución completa de `bash run_all.sh recover 2` en un clon limpio de la rama `a74368f` (15 min 38 s, 2 núcleos), comparada con `git diff` contra los resultados versionados:

- **Idénticos bit a bit:** `sequences/`, `alignments/`, `results/CLADE_SUPPORT.tsv`, `results/LEGACY_PHASE3A_CHECK.txt`.
- **Diferencias numéricas mínimas:** la log-verosimilitud varía en la 4ª–5ª cifra decimal (p. ej. −30700.65589 frente a −30700.65562) y las longitudes de rama en la 8ª. Es redondeo de coma flotante del optimizador multihilo de IQ-TREE.
- **Sin cambios** en ningún soporte bootstrap, p-valor (KH, SH, AU) ni topología.
