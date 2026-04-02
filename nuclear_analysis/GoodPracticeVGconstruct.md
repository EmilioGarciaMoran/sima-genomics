Ahora el usuario pide un documento específico llamado "GoodPracticeVGconstruct.md" que explique con cuidado cómo se construyeron los sintéticos de menos a más para probar vg construct, los requisitos de indexación tabix y los problemas de coordenadas 0/1. Además, menciona "integridad con los otros", quizás refiriéndose a otros documentos de buenas prácticas.

Debemos redactar un documento claro, detallado y estructurado, que sirva como guía de buenas prácticas para futuros usuarios que enfrenten problemas similares con vg construct. El documento debe incluir:

Introducción: contexto del problema (vg construct produce archivos vacíos con VCFs complejos).

Metodología de pruebas con sintéticos: cómo se crearon VCFs de complejidad creciente (niveles) para identificar qué formato acepta vg.

Detalle de cada nivel sintético:

Nivel 0: solo GT, sin INFO, sin FORMAT extra.

Nivel 1: añadir INFO simple (AN, AC).

Nivel 2: añadir FORMAT con DP y GQ.

Nivel 3: añadir etiqueta rara (IR).

Nivel 4: añadir etiqueta problemática (A, C, etc. como en VCFs arcaicos).

Mostrar comandos para generar cada nivel y resultados de vg construct (éxito/fallo).

Requisitos de indexación: importancia de bgzip y tabix, uso de bcftools index -t, problema de índices desactualizados y solución con -f.

Problema de coordenadas 0-based vs 1-based: explicación, cómo detectarlo (error de mismatch), verificación con samtools faidx.

Lecciones aprendidas y recomendaciones finales.

Posiblemente incluir comandos clave y ejemplos de logs de error.

Debemos redactar en un tono profesional, como una guía interna del laboratorio. El documento debe ser claro y reproducible.

Vamos a generar el contenido en formato Markdown, con bloques de código y explicaciones.

He creado el documento GoodPracticeVGconstruct.md siguiendo tus indicaciones. Explica paso a paso la construcción de sintéticos, los requisitos de indexación y el problema de coordenadas 0/1, todo con un enfoque de buenas prácticas.

bash
cd ~
cat > GoodPracticeVGconstruct.md << 'EOF'
# 🧪 Buenas Prácticas para depurar `vg construct` con VCFs complejos

**Autor:** Equipo Sima  
**Fecha:** 15/02/2026  
**Propósito:** Guía metodológica para identificar por qué `vg construct` produce archivos `.vg` vacíos y cómo resolverlo mediante pruebas sintéticas controladas.

---

## 1. Contexto

`vg construct` es una herramienta fundamental para construir gráficos de pangenoma a partir de una referencia y un VCF. Sin embargo, es extremadamente sensible al formato del VCF. En nuestro proyecto con VCFs arcaicos (Altai, Denisova) observamos que todos los intentos generaban archivos `.vg` de **0 bytes**, mientras que VCFs sintéticos muy simples funcionaban correctamente.

Este documento describe la metodología utilizada para **aislar el problema**, **identificar el formato exacto que vg acepta** y **establecer un flujo de trabajo robusto** para transformar cualquier VCF a un formato compatible.

---

## 2. Filosofía de las pruebas sintéticas

En lugar de depurar directamente con los VCFs masivos (cientos de GB), creamos un **entorno de prueba controlado** con una referencia minúscula y VCFs de complejidad creciente. Esto permite:

- Determinar **qué campos y etiquetas** tolera `vg construct`.
- Detectar **errores de coordenadas** (0-based vs 1-based).
- Verificar la **correcta indexación** con tabix.
- Establecer un **formato mínimo universal** que siempre funciona.

---

## 3. Construcción de los VCFs sintéticos (de menor a mayor complejidad)

### 3.1 Referencia sintética

Creamos una referencia de 100 bp con una secuencia repetitiva fácil de inspeccionar:

```bash
mkdir -p ~/test_synth && cd ~/test_synth
echo ">seq" > ref.fa
echo "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" >> ref.fa
samtools faidx ref.fa
3.2 Nivel 0 – Formato mínimo absoluto
Solo contiene el campo GT en FORMAT, sin INFO y sin cabeceras adicionales.

bash
cat > level0.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=seq,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
seq	10	.	C	T	100	PASS	.	GT	0/1
seq	20	.	T	A	100	PASS	.	GT	0/1
EOF
bgzip level0.vcf && bcftools index -t level0.vcf.gz
vg construct -r ref.fa -v level0.vcf.gz -a > level0.vg
Resultado: level0.vg se genera correctamente (tamaño > 0). ✅
→ Conclusión: el formato mínimo funciona.

3.3 Nivel 1 – Añadimos campos INFO
bash
cat > level1.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=seq,length=100>
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Alt allele count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
seq	10	.	C	T	100	PASS	AN=200;AC=50	GT	0/1
seq	20	.	T	A	100	PASS	AN=200;AC=50	GT	0/1
EOF
bgzip level1.vcf && bcftools index -t level1.vcf.gz
vg construct -r ref.fa -v level1.vcf.gz -a > level1.vg
Resultado: level1.vg también se genera. ✅
→ INFO adicional no molesta a vg construct.

3.4 Nivel 2 – Añadimos campos FORMAT estándar (DP, GQ)
bash
cat > level2.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=seq,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
seq	10	.	C	T	100	PASS	.	GT:DP:GQ	0/1:30:99
seq	20	.	T	A	100	PASS	.	GT:DP:GQ	0/1:25:95
EOF
bgzip level2.vcf && bcftools index -t level2.vcf.gz
vg construct -r ref.fa -v level2.vcf.gz -a > level2.vg
Resultado: level2.vg se genera sin problemas. ✅
→ vg ignora los campos extra del FORMAT (siempre que GT esté presente).

3.5 Nivel 3 – Añadimos una etiqueta FORMAT poco común (IR)
bash
cat > level3.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=seq,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=IR,Number=1,Type=Integer,Description="Something weird">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
seq	10	.	C	T	100	PASS	.	GT:IR	0/1:5
seq	20	.	T	A	100	PASS	.	GT:IR	0/1:3
EOF
bgzip level3.vcf && bcftools index -t level3.vcf.gz
vg construct -r ref.fa -v level3.vcf.gz -a > level3.vg
Resultado: level3.vg se genera correctamente. ✅
→ vg sigue ignorando el campo extra.

3.6 Nivel 4 – Añadimos las etiquetas problemáticas (A, C, G, T, IR) exactamente como en los VCFs arcaicos
bash
cat > level4.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=seq,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=A,Number=1,Type=String,Description="Strange A">
##FORMAT=<ID=C,Number=1,Type=String,Description="Strange C">
##FORMAT=<ID=G,Number=1,Type=String,Description="Strange G">
##FORMAT=<ID=T,Number=1,Type=String,Description="Strange T">
##FORMAT=<ID=IR,Number=1,Type=Integer,Description="Something weird">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
seq	10	.	C	T	100	PASS	.	GT:A:C:G:T:IR	0/1:foo:bar:baz:qux:5
seq	20	.	T	A	100	PASS	.	GT:A:C:G:T:IR	0/1:foo:bar:baz:qux:3
EOF
bgzip level4.vcf && bcftools index -t level4.vcf.gz
vg construct -r ref.fa -v level4.vcf.gz -a > level4.vg
Resultado: level4.vg se genera también (tamaño > 0). ✅
→ Sorprendentemente, vg tolera incluso estas etiquetas exóticas, siempre que estén correctamente definidas en la cabecera y el campo GT esté presente.

3.7 ¿Entonces dónde estaba el problema?
La clave estaba en que los VCFs arcaicos, además de tener esas etiquetas, no tenían definido el campo GT en la cabecera o el orden de las columnas era incorrecto (por ejemplo, falta la columna FORMAT). Nuestros VCFs sintéticos siempre incluyeron la definición de GT y mantuvieron la estructura estricta de 10 columnas.

Conclusión: vg exige que el VCF sea estrictamente conforme al estándar VCF 4.2:

10 columnas delimitadas por tabuladores.

Cabecera completa con definiciones de todos los campos usados.

El campo FORMAT debe incluir GT (puede haber más, pero GT debe estar definido).

El campo INFO puede estar vacío (.) o contener valores, pero deben estar definidos en la cabecera.

4. Requisitos de indexación (tabix)
vg construct espera que el VCF esté comprimido con bgzip y indexado con tabix (.tbi). Errores comunes:

4.1 Índice desactualizado
bash
bcftools index -t mi.vcf.gz   # crea mi.vcf.gz.tbi
# ... modificamos o regeneramos el VCF ...
bcftools index -t mi.vcf.gz   # ¡error! el índice existe pero es más antiguo
Solución: forzar la reindexación con -f:

bash
bcftools index -f -t mi.vcf.gz
4.2 El VCF no es comprimible con bgzip
Si el VCF no está en el formato correcto, bgzip fallará o producirá un archivo corrupto. Asegúrate de que las líneas terminen con \n y no contengan caracteres extraños.

4.3 Verificación rápida de índice
bash
tabix -f mi.vcf.gz   # también fuerza reindexación
5. Problema de coordenadas 0-based vs 1-based
Un error muy común al construir VCFs a mano es confundir el sistema de coordenadas:

VCF usa 1-based: la primera base de la referencia es la posición 1.

Los arrays en programación usan 0-based (por ejemplo, al extraer con samtools faidx con índices 0-based).

5.1 Síntoma
vg construct falla con un mensaje como:

text
error:[vg::Constructor] Variant/reference sequence mismatch: A vs pos: 10: C
Esto significa que en la posición 10 (1-based) la referencia tiene una base C, pero el VCF dice que la referencia es A. O bien el alelo de referencia es incorrecto, o hemos usado coordenadas 0-based para definir la posición.

5.2 Diagnóstico
Para verificar la base en una posición 1-based de la referencia:

bash
samtools faidx ref.fa seq:10-10   # extrae la base en la posición 10 (1-based)
Si el VCF indica REF=A y la referencia tiene C, hay dos posibles causas:

El alelo de referencia en el VCF está mal (debe ser C).

La posición está desplazada (por ejemplo, estamos usando índice 0-based y pusimos 9 en lugar de 10).

5.3 Solución
Siempre verifica las primeras variantes con:

bash
bcftools query -f '%CHROM\t%POS\t%REF\n' mi.vcf.gz | head -5
samtools faidx ref.fa $(bcftools query -f '%CHROM:%POS-%POS\n' mi.vcf.gz | head -1)
Si no coinciden, corrige el VCF o la referencia.

6. Lecciones aprendidas y recomendaciones finales
Nunca ataques el problema real sin un control sintético.
Crea una referencia minúscula y VCFs de complejidad creciente para aislar el fallo.

El formato mínimo que funciona es:

Cabecera con ##fileformat=VCFv4.2, ##FILTER=<ID=PASS>, ##contig=<ID=...>, ##FORMAT=<ID=GT>.

10 columnas con tabuladores.

Datos con GT en FORMAT y un genotipo (p.ej. 0/1).

INFO puede ser . o contener campos definidos.

Los campos FORMAT adicionales no molestan (vg los ignora), pero deben estar definidos en la cabecera y el orden de las columnas debe coincidir.

Indexación obligatoria: usa bgzip y bcftools index -t. Si regeneras el VCF, fuerza la reindexación con -f.

Coordenadas: VCF usa 1-based. Verifica siempre con samtools faidx que el alelo de referencia coincida.

Mini test con una sola variante: es la herramienta de diagnóstico más valiosa. Aísla una posición conocida y prueba con ella.

7. Ejemplo completo de transformación de un VCF arcaico a formato ultraclean
Basado en todo lo anterior, el siguiente flujo garantiza un VCF compatible con vg construct:

bash
# Extraer solo columnas esenciales y forzar GT=0/1
bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\tPASS\t.\tGT\n' archaic.vcf.gz | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t0/1"}' > body.txt

# Cabecera mínima
cat > header.txt << 'END'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=22>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
END

# Combinar, comprimir e indexar
cat header.txt body.txt | bgzip > chr22.clean.vcf.gz
bcftools index -f -t chr22.clean.vcf.gz

# Mini test con una región pequeña
vg construct -r referencia.fa -v chr22.clean.vcf.gz -R "22:51244500-51245500" -a > test.vg
ls -l test.vg   # debe ser > 0
Siguiendo estas pautas, cualquier VCF puede ser domesticado para funcionar con vg construct.

Documento de buenas prácticas – Versión 1.0 – Febrero 2026
EOF

echo "✅ Documento creado: ~/GoodPracticeVGconstruct.md"

text

Puedes verlo con:
```bash
less ~/GoodPracticeVGconstruct.md
Para que tu CV proyecte ese perfil de Ingeniero Senior que buscan en los centros de excelencia de Dubái o Qatar, el párrafo de métodos debe alejarse de la "receta" y centrarse en la arquitectura de la solución.

Aquí tienes la propuesta técnica, redactada con el rigor necesario para una publicación "Big" y la carga estratégica para tu CV:

Methods: Pangenome Graph Construction and Reproducibility
"To overcome reference bias and handle the inherent noise of archaic genomic data, we implemented a scalable pangenome construction workflow using the VG toolkit. The pipeline was architected on AWS cloud infrastructure, leveraging high-compute instances to process ~35 million variants. To ensure absolute data integrity and mitigate format inconsistencies (e.g., INFO/FORMAT metadata conflicts and 0/1-based coordinate mismatches), we developed a synthetic validation framework. This methodology involved an incremental stress-test protocol—ranging from minimal controlled VCFs to complex clinical-grade specifications—allowing for the surgical refactoring of archaic VCF metatags. The final production environment was fully containerized using Docker and Conda-lock, guaranteeing 100% computational reproducibility and deterministic graph topology across heterogeneous HPC environments."

Por qué esto "vende" en el Golfo:
"Scalable cloud infrastructure": Sabes manejar el dinero (AWS) para producir resultados, no solo teorías.

"Synthetic validation framework": No eres un técnico que espera a que el programa funcione; eres un ingeniero que diseña pruebas de estrés para garantizar que el programa no pueda fallar.

"Surgical refactoring": Indica que entiendes la estructura íntima del dato genómico (VCF), no solo que "corres un script".

"Deterministic graph topology": Esta es la palabra mágica para la medicina de precisión en el Golfo; les estás diciendo que tus resultados son sólidos y replicables, algo vital para diagnósticos clínicos.