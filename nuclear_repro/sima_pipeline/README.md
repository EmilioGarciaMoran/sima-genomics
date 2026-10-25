# Pipeline de Sima de los Huesos — ADN nuclear

Reproducción del recuento por rama de Meyer et al. 2016 con los datos de ENA PRJEB10597,
y D-statistics sobre los conjuntos de referencia construidos en `../branch_counts/`.

**Preregistro:** `../PREREGISTRO_D_statistics.md`. Hipótesis, umbrales, controles y reglas de
interpretación quedaron fijados antes de tocar estos datos.

## Estado

| Paso | Estado |
|---|---|
| 1. Metadatos de ENA: qué espécimen es cada corrida | `prepare.sh` |
| 2. Referencia GRCh37 e índice de bwa | `prepare.sh` |
| 3. Descarga de los FASTQ | `prepare.sh` |
| 4. Mapeo con `bwa aln` sin siembra | pendiente |
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
