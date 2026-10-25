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

Control **X = Denisova 3** (denisovano, grupo externo a los dos neandertales): no se inclina
hacia ninguno. f_A = 0,140 y f_B = 0,120 (Z = 1,0) en pseudohaploide; 0,104 y 0,126
(Z = −0,7) en homocigoto.

Las tres condiciones preregistradas pasan en los tres modos: f_B > f_A, f_AB es el mayor de
los tres, y el control denisovano no se inclina.

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
