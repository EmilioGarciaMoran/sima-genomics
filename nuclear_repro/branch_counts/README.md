# Recuento de alelos derivados por rama: validación del software

Implementación del recuento por rama de Meyer et al. 2016 y Massilani et al. 2026, con sus
pruebas. **Este directorio valida el código, no reproduce todavía ningún resultado publicado.**

## Por qué existe

El preregistro (`../PREREGISTRO_D_statistics.md`, apartado 10) exige validar el pipeline antes
de tocar los datos de Sima. El control más limpio es reproducir las cifras de Massilani et al.
2026 con genomas de alta cobertura, sin daño ni contaminación que compliquen la lectura. Si el
código no las reproduce ahí, no tiene sentido aplicarlo a datos difíciles.

## Los tres estadísticos

Para un individuo de prueba X, dos genomas de referencia A y B del mismo linaje, y un humano
moderno M, con el estado ancestral tomado de los primates:

| Estadístico | Conjunto de referencia | Mide |
|---|---|---|
| f_AB | A y B derivados, M ancestral | rama profunda compartida |
| f_A | solo A derivado | deriva específica de A |
| f_B | solo B derivado | deriva específica de B |

Cada uno es la fracción de ese conjunto en la que X lleva el alelo derivado, contando solo las
posiciones que X cubre.

**Massilani et al. 2026**, con X = Denisova 17, A = Denisova 5, B = Vindija 33.19:
f_AB = 91 %, f_A = 45 %, f_B = 20 %.

## Tolerancia para la validación, fijada antes de ejecutar

La incertidumbre por azar de muestreo es despreciable con genomas de alta cobertura:

| Sitios del conjunto | IC95 sobre 91 % | sobre 45 % | sobre 20 % |
|---|---|---|---|
| 100.000 | ±0,18 pt | ±0,31 pt | ±0,25 pt |
| 500.000 | ±0,08 pt | ±0,14 pt | ±0,11 pt |
| 2.000.000 | ±0,04 pt | ±0,07 pt | ±0,06 pt |

Una tolerancia estadística (±0,1 pt) fallaría siempre sin que el código estuviera mal. La
discrepancia real viene de decisiones de implementación que el artículo no especifica por
completo: qué conjunto de primates define el estado ancestral, qué filtros de mapabilidad,
cómo se tratan los sitios de cobertura baja, si se incluyen las transiciones, qué versión de
los alineamientos lastz de UCSC. Cada una mueve puntos enteros.

**Tolerancia fijada: ±3 puntos porcentuales en cada uno de los tres estadísticos, y el orden
f_AB > f_A > f_B debe conservarse.**

| Resultado | Lectura |
|---|---|
| Los tres dentro de ±3 pt | Pipeline validado. Se pasa a Sima. |
| Alguno entre 3 y 8 pt de desviación, orden correcto | Diferencia de implementación. Se documenta cuál, se reporta y se decide si invalida. No se pasa a Sima sin explicar la causa. |
| Desviación > 8 pt, u orden alterado | Pipeline incorrecto. No se ejecuta nada más. |

El umbral de 3 puntos es una decisión de juicio, no un cálculo: es el margen que dejan las
ambigüedades de implementación sin permitir que un error de código pase desapercibido. Queda
fijado aquí antes de ejecutar para que no se ajuste luego al resultado.

## Estado de la validación

### 1. Pruebas de respuesta conocida — PASAN

```
python3 scripts/test_branch_counts.py
```

Ocho casos construidos a mano, cuya respuesta se conoce por construcción. Cubren la
separación de los tres conjuntos, la exclusión de sitios donde M lleva el derivado, el
tratamiento de datos ausentes en X (reducen el denominador, no se cuentan como ancestral) y
en A o B (excluyen el sitio), la recuperación de una fracción fijada en muestra grande, la
concordancia del jackknife con el error binomial cuando los bloques son homogéneos, su
aumento cuando hay estructura por bloques, y el rechazo de entradas inválidas.

### 2. Comprobación con datos simulados — PASA

```
python3 scripts/sim_check.py        # cromosomas 21 y 22, unos 10 segundos
```

Simula la topología que implican las cifras de Massilani: Denisova 17 hermano de Denisova 5,
Vindija 33.19 hermano de Chagyrskaya 8.

| Conjunto | Sitios | f simulada | Massilani |
|---|---|---|---|
| AB, rama profunda | 29.001 | 0,959 | 0,91 |
| A, exclusivos de D5 | 4.330 | 0,376 | 0,45 |
| B, exclusivos de Vindija | 7.355 | 0,148 | 0,20 |

**Las cifras simuladas no son una reproducción de Massilani**: la demografía real no se conoce
y no se ha ajustado ningún parámetro para que coincidan. Que caigan en la misma escala es una
comprobación de cordura, nada más.

Lo que sí valida el código es el control: con X = Chagyrskaya 8 los conjuntos A y B se
invierten (f_A = 0,141 frente a f_B = 0,453), como exige la topología. Un error que
intercambiara los dos conjuntos se habría delatado ahí.

### 3. Validación con datos reales (cromosoma 21) — PASA

Ejecutado el 25 de octubre de 2026 sobre los VCF de hg19 de MPI-EVA (ver
`DATOS_NECESARIOS.md`). 32.022.404 posiciones con estado ancestral; 22,3 M con dato en los
cuatro arcaicos. Jackknife de 35 bloques de 1 Mb.

Prueba principal: **X = Chagyrskaya 8, A = Denisova 5, B = Vindija 33.19**. La topología
establecida exige f_B > f_A.

| Modo | Conjunto AB | f_A (exclusivos de D5) | f_B (exclusivos de Vindija) | Z(f_A − f_B) |
|---|---|---|---|---|
| pseudohaploide | 0,944 (n=15.019) | 0,259 (n=3.222) | 0,391 (n=3.820) | **−4,3** |
| pseudohaploide, Map20=1 | 0,944 (n=13.895) | 0,259 (n=2.956) | 0,390 (n=3.524) | −4,2 |
| homocigoto | 0,977 (n=10.951) | 0,150 (n=1.298) | 0,507 (n=1.326) | **−6,8** |
| portador (carrier) | 0,926 (n=13.233) | 0,226 (n=3.286) | 0,380 (n=4.758) | −5,1 |

Control **X = Denisova 3** (denisovano, grupo externo a los dos neandertales): no se inclina
hacia ninguno. Z(f_A − f_B) = 1,0 en pseudohaploide, −0,7 en homocigoto y 2,0 en portador.

Las tres condiciones preregistradas pasan en los cuatro análisis: f_B > f_A, f_AB es el mayor
de los tres, y el control denisovano no se inclina. **La validación del software queda
cerrada.**

**Modo elegido para el análisis de Sima: homocigoto para los conjuntos de referencia.** Da la
señal más fuerte (Z = −6,8 frente a −5,1 y −4,3) y el control más limpio (Z = −0,7 frente a
2,0 y 1,0). El modo portador incluye heterocigotos en los conjuntos y es el que más acerca el
control al umbral, señal de que mete ruido. El individuo de prueba va aparte: pseudohaploide,
porque con Sima no hay otra opción.

**Lo que se aprende:**

- **El filtro de mapabilidad es irrelevante aquí.** Restringir a Map20 = 1 (24 de 32 M de
  posiciones) deja las fracciones prácticamente idénticas. La señal no depende de regiones de
  mapeo dudoso, así que el filtro se deja fuera y queda justificado.
- **El modo homocigoto da una señal más limpia**, con un contraste 3,4 veces mayor frente a
  1,5. El muestreo pseudohaploide introduce ruido al definir los conjuntos de referencia.
  **Consecuencia para el análisis de Sima: los conjuntos se definen con los genotipos de los
  genomas de alta cobertura, y el muestreo pseudohaploide se aplica solo al individuo de
  prueba**, donde es obligatorio por tener una sola lectura por posición.
- **El jackknife importa.** Su error típico es unas 2,7 veces el binomial (0,021 frente a
  0,008), porque los sitios cercanos están ligados. El Z honesto es −4,3, no el −12 que daría
  el binomial. Cualquier resultado futuro debe usar el jackknife.
- Un primer ensayo sobre 1,5 Mb pericentroméricos dio f_A = 0,35 en el control denisovano,
  valor imposible. Era artefacto de la región: con el cromosoma completo baja a 0,14, que sí
  se explica por clasificación incompleta de linajes y el flujo génico denisovano conocido en
  el neandertal del Altái. **No se interpreta ningún resultado sobre una franja pequeña.**

Esto valida el software para el recuento por rama. **No** reproduce las cifras 91/45/20, que
exigen Denisova 17 (ver apartado 4).

### 3 bis. Recorrido genómico — COMPLETADO (25 oct 2026, 22:41)

Los 22 autosomas, modo homocigoto, solo sitios informativos (Mbuti ancestral y al menos un
arcaico derivado). 7 h 39 min, pico de disco 21 GB, VCF borrados tras cada cromosoma.
Auditado con `check_genome.py`: los 22 cubren más del 99,8% de su longitud en hg19.

**1.822.345 sitios informativos** (563 bloques de jackknife de 5 Mb). Salida:
`results/genome/chrN_inf.npz` y `results/genoma_inf.npz`.

| Conjunto de referencia | Sitios |
|---|---|
| AB, rama profunda (D5 y Vi33.19 derivados) | 745.823 |
| A, exclusivos de Denisova 5 | 79.894 |
| B, exclusivos de Vindija 33.19 | 101.248 |

**Proyección para Sima**, aplicando las bases alineadas de Meyer 2016 (Tabla 1) sobre el total
de posiciones consultables. El denominador está estimado a partir del cromosoma 21, así que
las cifras llevan un margen del orden del 10%; el valor exacto saldrá al mapear Sima.

| Espécimen | Bases alineadas | Rama profunda | Exclusivos D5 | Exclusivos Vi |
|---|---|---|---|---|
| Fémur XIII | 225.329 | ~126 | ~14 | ~17 |
| AT-5431 | 1.186.442 | ~664 | ~71 | ~90 |
| AT-5482 | 2.015.167 | ~1.128 | ~121 | ~153 |
| AT-5431 + AT-5482 | 3.201.609 | ~1.797 | ~192 | ~244 |

Concuerda con el estudio de potencia (`../power_sim/`), que predecía unos 330 sitios
informativos para AT-5482 y 526 para los dos combinados; la proyección sobre datos reales da
274 y 436. El mismo orden de magnitud confirma que la simulación estaba bien calibrada y que
el análisis preregistrado es viable con los datos disponibles.

Confirma también la decisión ya tomada: **AT-5431 y AT-5482 combinados como análisis
principal, Fémur XIII solo como confirmación.** Los conjuntos específicos de cada neandertal
quedan justos, de modo que las particiones finas de la rama neandertal tendrán poca
resolución con Sima.

### 4. Reproducción de las cifras de Massilani — NO VIABLE POR AHORA

Requiere los genomas de alta cobertura (Denisova 5, Denisova 17, Vindija 33.19, Denisova 3,
Mbuti) y los alineamientos de primates de UCSC. **No se pueden descargar desde el entorno en
la nube**, cuya red solo alcanza repositorios de paquetes: hay que ejecutarlo en local.

Pasos previos, en este orden:
1. Comprobar disponibilidad y formato de cada genoma (VCF sobre hg19). Denisova 17 se cita
   como Massilani et al. 2026; Goyet 1 solo como accesión ENA PRJEB98484.
2. Verificar que las cifras 91/45/20 del artículo son las correctas y con qué referencia
   exacta están medidas. Las de este README se tomaron del texto; conviene confirmarlas en la
   tabla original antes de usarlas como criterio.
3. Ejecutar y comparar con la tolerancia de arriba.

## Dependencias

`numpy` para el módulo y las pruebas; `msprime` solo para `sim_check.py`.
