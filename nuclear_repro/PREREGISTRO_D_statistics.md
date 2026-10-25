# Preregistro: D-statistics con los datos nucleares de Sima de los Huesos

**Fecha de cierre de este documento:** 25 de octubre de 2026
**Estado:** fijado ANTES de ejecutar cualquier análisis con los datos reales de PRJEB10597.
**Commit:** este fichero se sube antes que cualquier script o resultado del análisis real. Si la
fecha de un resultado precede a la de este commit, ese resultado no cuenta como preregistrado.

---

## 1. Pregunta

¿Pertenece el linaje de Sima de los Huesos al linaje neandertal, o es un linaje hermano de
neandertales y denisovanos en conjunto?

La segunda posibilidad es la que haría defendible una taxonomía propia (*Homo atapuercae*).
La primera es la publicada por Meyer et al. 2016.

## 2. Estadístico

D(Denisova 3, Vindija 33.19; X, chimpancé), con X = cada espécimen de Sima.

- Alelos pseudohaploides: uno al azar por posición y por individuo.
- Error típico por jackknife de bloques ponderado (Busing et al. 1999), bloques de 5 Mb.
- Grupo externo: chimpancé (*Pan troglodytes*). **No** un humano moderno: Ust'-Ishim tiene
  mezcla neandertal y su uso fue uno de los errores del preprint retirado.

## 3. Hipótesis y predicciones

| | Posición de Sima | Predicción |
|---|---|---|
| **H1** (Meyer 2016) | Dentro del linaje neandertal | D > 0, hacia Vindija |
| **H2** (linaje aparte) | Hermano de neandertales + denisovanos | D ≈ 0 |

**H2 es la hipótesis que predice D = 0.** Por tanto una D significativa hacia Vindija no apoya
una taxonomía propia: la descarta. Queda constancia de que el resultado "positivo" del análisis
va en contra de la propuesta que motivó el proyecto.

## 4. Umbral de significación

**|Z| > 3** (p ≈ 0,0027 bilateral). No 1,96. Razones:

1. El error típico del jackknife puede subestimarse si los bloques de 5 Mb no son del todo
   independientes.
2. Se harán varias pruebas (tres especímenes, combinaciones, varios neandertales de referencia).
   Equivale aproximadamente a una corrección de Bonferroni para 10–20 pruebas.
3. Hay sesgos sistemáticos (daño, contaminación moderna, sesgo hacia la referencia) que un Z
   alto no corrige.

**Con n < 50 sitios informativos** la aproximación normal no es fiable: se usará además una
prueba binomial exacta sobre ABBA frente a BABA, y se reportarán las dos.

## 5. Controles obligatorios

Un análisis cuyos controles fallen se descarta íntegro, sin interpretar el resultado de Sima.

| Control | Criterio |
|---|---|
| **Negativo:** X = Mbuti (humano moderno) | |Z| < 3. Si da señal, hay sesgo de mapeo o contaminación. |
| **Positivo:** X = Chagyrskaya 8 submuestreado al volumen de Sima | Z > 3. Si no da señal, no hay potencia y el resultado de Sima no es interpretable. |
| **Transversiones:** repetir solo con transversiones | La D no debe cambiar de signo ni perder más de la mitad de su magnitud. |
| **Por espécimen:** cada espécimen por separado | No fusionar BAM de distintos especímenes. Fue uno de los errores de 2026. |

## 6. Volumen de datos y potencia esperada

Del estudio de potencia (`nuclear_repro/power_sim/`, commits `b60b391`, `369477b`, `569949a`,
`b29c846`). Bases alineadas tras filtro de daño según Meyer 2016, Tabla 1.

| Espécimen | Bases alineadas | Sitios informativos esperados | Potencia (escenario calibrado) |
|---|---|---|---|
| Fémur XIII | 225.329 | ~37 | 94% |
| AT-5431 (fémur) | 1.186.442 | ~194 | 100% |
| AT-5482 (incisivo) | 2.015.167 | ~331 | 100% |
| AT-5431 + AT-5482 | 3.201.609 | ~526 | 100% |
| Molar | 188.974 | ~31 | no se usa como prueba principal |
| Escápula | 37.524 | ~6 | excluida (contaminación >98%) |

**Análisis principal:** AT-5431 y AT-5482, por separado y combinados.
**Confirmación:** Fémur XIII, con prueba binomial exacta.
**No se usa:** escápula; molar solo descriptivo.

## 7. D mínima detectable

D = Z / √(n + Z²), con Z = 3.

| Sitios informativos | D mínima detectable |
|---|---|
| 37 | 0,44 |
| 194 | 0,21 |
| 331 | 0,16 |
| 526 | 0,13 |

Solo con AT-5431 y AT-5482 combinados podríamos afirmar D < 0,13, que es el único resultado
que dejaría H2 en pie.

## 8. Predicción cuantitativa a partir de Meyer 2016

Derivada de los porcentajes publicados (rama neandertal 43% en AT-5431 y 39% en AT-5482; rama
denisovana 9% y 7%) combinados con la cobertura por rama de la simulación.

| Espécimen | ABBA | BABA | n | D predicha | Z predicho |
|---|---|---|---|---|---|
| Fémur XIII | 32 | 6 | 38 | 0,68 | 5,8 |
| AT-5431 | 171 | 35 | 206 | 0,66 | 12,6 |
| AT-5482 | 264 | 46 | 310 | 0,70 | 17,4 |
| AT-5431 + AT-5482 | 435 | 80 | 515 | 0,69 | 21,6 |

**Esto es una predicción, no un resultado.** Su valor está en que fija el blanco antes de
disparar: el análisis real es una comprobación genuina y no un ejercicio de confirmación.

### 8.1 Coherencia interna de los porcentajes de Meyer

Comprobación hecha el 25 oct 2026, antes del análisis real:

| Espécimen | f_N/f_D | IC 95% |
|---|---|---|
| AT-5431 | 4,78 | 3,41 – 6,70 |
| AT-5482 | 5,57 | 4,15 – 7,48 |
| Combinado | **5,21** | 4,17 – 6,51 |

Los dos especímenes no difieren (Z = −0,67; p = 0,50), pese a tener coberturas distintas y
haber sido procesados por separado. Contrastado con los escenarios simulados: 440 ka
incompatible (Z = 13,3), 470 ka incompatible (Z = 9,6), **550 ka compatible (Z = 0,3)**.
Los escenarios de poca potencia quedan descartados por los propios datos publicados.

## 9. Reglas de interpretación

Fijadas aquí para que un resultado intermedio no se interprete a conveniencia.

| D real observada | Lectura | Acción |
|---|---|---|
| **≥ 0,50** | Confirma Meyer 2016: Sima dentro del linaje neandertal | Publicar como reproducción independiente. H2 descartada por la vía genética. |
| **0,13 – 0,50** (significativa) | Significativa pero muy por debajo de lo predicho | **No interpretar biológicamente.** Investigar primero contaminación, sesgo hacia la referencia y la traducción de los porcentajes de Meyer al estadístico con Vindija. |
| **< 0,13** (no significativa) | Incompatible con la predicción | Revisar el razonamiento y los datos de partida antes de reclamar nada. No se presenta como apoyo a H2 sin haber descartado fallos técnicos. |
| Signo invertido (hacia Denisova) | No predicho por ninguna hipótesis | Error técnico casi seguro. Revisar antes de interpretar. |

## 10. Reconfirmación de los porcentajes de Meyer

Los números del apartado 8 se toman del artículo. El paso que los reconfirma es
independiente y previo al análisis nuevo:

1. **Reproducir el recuento por rama de Meyer 2016** con su propio método, y comprobar si
   salen 87/43/9 (AT-5431) y 68/39/7 (AT-5482).
   - Mapeo con `bwa aln` con siembra desactivada, contra hg19.
   - Filtro de daño: solo fragmentos con C→T terminal; T enmascaradas en las 3 primeras y
     3 últimas posiciones.
   - Sitios informativos definidos con los VCF de alta cobertura y los alineamientos lastz de
     UCSC, no comparando BAM entre sí.
   - Cada espécimen por separado.
2. **Criterio de validación del pipeline:** los porcentajes reproducidos deben caer dentro del
   IC 95% de los publicados. Si no, el pipeline está mal y no se ejecuta nada más.
3. **Solo si pasa 2** se calculan las D-statistics nuevas.

**Limitación conocida:** los FASTQ de ENA (PRJEB10597) son pequeños (p. ej. ERR995360:
69.404 lecturas) y casi seguro contienen solo las lecturas que ya mapeaban contra hg19. Los
porcentajes reproducidos saldrán de un subconjunto. Un cociente f_N/f_D es relativamente
robusto a eso, pero los porcentajes absolutos pueden no coincidir. Si ocurre, se reporta la
discrepancia y se usa el cociente, no los absolutos.

**Anclaje independiente del método:** Massilani et al. 2026 aplica el mismo recuento por rama
a Denisova 17 (91% de alelos derivados compartidos por ambos neandertales; 45% de los
exclusivos de D5; 20% de los de Vindija). Reproducir esas cifras valida el método sin tocar
los datos de Sima, y es un control previo más barato.

## 11. Qué NO responde este análisis

La D ubica a Sima en el árbol. **No** mide una distancia que por sí sola justifique una especie
nueva. Para eso haría falta particionar la rama neandertal en profunda y reciente usando los
cinco genomas neandertales de alta cobertura (D5, D17, Chagyrskaya 8, Vindija 33.19, Goyet 1),
y preguntar en qué punto de esa rama se separa Sima. Eso mide divergencia dentro del linaje,
no pertenencia a él, y sigue siendo una pregunta de genética de poblaciones, no de taxonomía.

## 12. Cifras retiradas

D = 0,7425 con Z = 30,6 sobre 218.757.031 sitios (preprint de 2026): sin respaldo, retiradas.
El Z máximo alcanzable con los volúmenes reales de Sima, según simulación, fue 7,8. Un Z de
30 con D = 0,74 requeriría unos 800 sitios informativos, cuatro veces más de lo disponible.
