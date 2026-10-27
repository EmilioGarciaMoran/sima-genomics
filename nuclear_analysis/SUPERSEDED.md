# Material superado — análisis nuclear de 2026

> **Esta carpeta se conserva solo por trazabilidad. Nada de lo que contiene respalda
> ninguna conclusión vigente.** El análisis nuclear actual está en
> [`../nuclear_repro/`](../nuclear_repro/).

## Qué hay aquí y por qué no vale

| Fichero | Estado |
|---|---|
| `GoodAnalysisPostAWS.md`, `GoodPracticeVGconstruct.md` | Notas operativas del trabajo en AWS de febrero de 2026. Útiles como registro de cómo se trabajó; **sus afirmaciones sobre resultados obtenidos no están respaldadas** (ver abajo). |
| `fastqs/Sima_de_los_Huesos/Sima_FemurXIII.hg19.*` | `flagstat` y `depth` de un mapeo cuyos parámetros no son adecuados para ADN antiguo. |
| `fastqs/descargar.sh` | Script de descarga de aquella época. Sustituido por `../nuclear_repro/sima_pipeline/scripts/prepare.sh`. |

## Por qué quedó retirado (octubre de 2026)

1. **No existe ninguna salida que respalde las cifras publicadas.** Una búsqueda exhaustiva
   en S3 y en las dos instancias EC2, antes de darlas de baja el 24 de octubre de 2026, no
   encontró ningún fichero de ANGSD (`.arg`, `.abbababa2`, `estAvgError`), ni llamadas a esa
   herramienta en el historial de órdenes. Las cifras del preprint —D = 0,7425, Z = 30,6
   sobre 218.757.031 sitios— quedan **retiradas**.

2. **Un Z de ese orden es inalcanzable con los datos de Sima.** El estudio de potencia
   ([`../nuclear_repro/power_sim/`](../nuclear_repro/power_sim/)) sitúa el máximo alcanzable
   en torno a 7,8 con los volúmenes reales. El análisis preregistrado posterior obtuvo
   Z = 9,5 sobre 174 sitios, no sobre 218 millones.

3. **Problemas de método en aquel trabajo**, documentados para que no se repitan:
   - mapeo con `bwa mem` sin parámetros de ADN antiguo y **sin filtro de daño**;
   - las cinco corridas de ENA fusionadas en un solo BAM, mezclando especímenes distintos,
     incluida la escápula AT-6672, contaminada al 98%;
   - llamada diploide sobre datos de cobertura muy inferior a 1x;
   - grupo externo Ust'-Ishim, un humano moderno con ascendencia neandertal, en lugar del
     chimpancé.

Nada de esto invalida el esfuerzo ni las notas operativas, que describen con precisión cómo
se montó la infraestructura. Invalida los resultados.
