# ¿Tienen potencia las D-statistics con el volumen de datos de Sima? — Simulación

**Estado: diseño y criterio fijados ANTES de ejecutar** (este fichero se versiona en un commit anterior al de los resultados).

## Pregunta
Con el volumen real de ADN nuclear de Sima de los Huesos (Meyer et al. 2016, Tabla 1: bases alineadas tras el filtro de daño), ¿detecta una D-statistic con chimpancé como grupo externo la afinidad neandertal de un individuo de hace ~430 ka, sin dar falsos positivos?

No se pregunta si las D-statistics "funcionan" en general (es un método establecido), sino si **este volumen de datos** les da potencia.

## Método
Simulación coalescente de genomas completos (msprime) bajo una demografía conocida. Después, degradación de los genomas hasta el volumen de Sima y cálculo de la D en 100 réplicas por condición. Como la verdad se conoce, se mide con qué frecuencia la D la recupera.

### Demografía (supuestos; valores redondeados de la literatura)
| Parámetro | Valor principal | Fuente orientativa |
|---|---|---|
| Tiempo generacional | 29 años | convención habitual |
| Tasa de mutación | 1,25 × 10⁻⁸ /pb/generación | convención habitual |
| Tasa de recombinación | 1 × 10⁻⁸ /pb/generación | convención habitual |
| Separación humanos modernos–arcaicos | 600 ka | Prüfer et al. 2014; Meyer et al. 2016 |
| **Separación neandertal–denisovano (T_ND)** | **470 ka** (sensibilidad: 440 ka y 550 ka) | Meyer et al. 2016 (anterior a ~430 ka); Prüfer et al. 2017 |
| Ne africanos actuales / ancestro humano-arcaico | 15.000 / 10.000 | — |
| Ne ancestro arcaico / neandertales / denisovanos | 5.000 / 3.000 / 3.000 | heterocigosidad baja de los arcaicos (Prüfer 2014; Massilani 2026) |

Simplificaciones: sin flujo génico entre linajes, poblaciones de tamaño constante entre eventos, genoma autosómico de longitud hg19, modelo de sitios infinitos.

**Grupo externo:** el chimpancé no se simula. Con sitios infinitos, el alelo ancestral de cada mutación se conoce con exactitud, así que se usa directamente como O. Equivale a un chimpancé que polariza sin errores. En los datos reales la polarización con chimpancé tiene una pequeña tasa de error que aquí no se modela. Esta decisión se tomó tras una prueba de tiempo: con la rama del chimpancé, el cromosoma 21 tardaba 227 s; sin ella, 11 s.

### Muestras (diploides)
| Muestra | Población | Edad de muestreo |
|---|---|---|
| Mbuti | africanos actuales | 0 |
| Denisova 3 | denisovanos | 65 ka |
| Vindija 33.19 | neandertales | 50 ka |
| Chagyrskaya 8 | neandertales | 80 ka |
| **Sima** | **linaje neandertal** (después de T_ND) | **430 ka** |

### Estadístico
D(H1 = Denisova 3, H2 = Vindija; X, O = estado ancestral) = (ABBA − BABA) / (ABBA + BABA), con **D > 0 si X comparte más alelos derivados con Vindija (neandertal)**. Todas las muestras son pseudohaploides: un alelo al azar por sitio. Error típico por jackknife de bloques ponderado (Busing et al. 1999), con bloques de 5 Mb. Z = D / EE.

### Volúmenes de datos (bases por individuo X)
| Nivel | Bases | Espécimen de referencia |
|---|---|---|
| femurXIII | 225.329 | Fémur XIII |
| AT5431 | 1.186.442 | fémur AT-5431 |
| AT5482 | 2.015.167 | incisivo AT-5482 |
| full | genoma completo | referencia sin degradar |

Cada base del individuo X se conserva con probabilidad p = bases / longitud del genoma simulado, y se toma un alelo al azar. Las otras muestras (H1, H2 y O) no se degradan, como corresponde a genomas de alta cobertura.

### Condiciones
- **Prueba (Sima):** X = individuo del linaje neandertal muestreado a 430 ka. Verdad: afinidad neandertal (D > 0).
- **Control positivo fácil:** X = Chagyrskaya 8. Verdad: D > 0 y grande.
- **Control negativo:** X = Mbuti. Verdad: D = 0 (sin flujo génico simulado).

## Criterio de decisión (fijado de antemano)
Para el escenario principal (T_ND = 470 ka), en cada volumen de datos:
- **Las D-statistics son apropiadas** si la prueba (Sima) da **Z > 3 con signo correcto en ≥ 90%** de las 100 réplicas **y** el control negativo (Mbuti) da **|Z| < 3 en ≥ 95%**.
- **Potencia insuficiente** si la prueba da Z > 3 en < 90%. En ese caso hay que usar el recuento de alelos derivados por rama (Meyer 2016) o combinar datos.
- Los escenarios de sensibilidad (T_ND = 440 y 550 ka) se informan siempre; si cambian la conclusión, se dice expresamente.

## Nota de transparencia
Antes de este commit se ejecutó una prueba de funcionamiento con un genoma de juguete (22 cromosomas de 2 Mb, escenario T_ND = 470 ka, 5 réplicas) para depurar el código. En ella el Fémur XIII mostraba potencia insuficiente y AT-5482 suficiente. El criterio de decisión de arriba se redactó antes de esa prueba y **no se ha modificado** tras verla.

## Limitaciones conocidas (se informarán con los resultados)
Mide **potencia estadística**: cuántos sitios informativos hay y si bastan. No modela el daño post mortem, la contaminación moderna, el sesgo de mapeo ni la estructura poblacional real de los arcaicos. Un resultado favorable aquí es necesario pero no suficiente; la comprobación con datos reales vendrá después.

## Reproducir
```bash
cd nuclear_repro/power_sim
micromamba create -f environment.yml && micromamba activate powersim
bash run_power.sh 2        # 2 = hilos
```
