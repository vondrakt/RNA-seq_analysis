# This is the reference based RNA-seq pipelien

# Due to incompatibility with rnaspades, this step is skipped

SAMPLES=["A","B","C"]

# Trimming

rule trimming:
	input:
		p1 = expand("17NQ004_P1_R11{sample}.fastq", sample=SAMPLES),
		p2 = expand("17NQ004_P2_R5{sample}.fastq", sample=SAMPLES)
	output:
		p1p = expand("17NQ004_P1_R11{sample}_1P.trimmed", sample=SAMPLES),
		p1u = expand("17NQ004_P1_R11{sample}_1U.trimmed", sample=SAMPLES),
		p2p = expand("17NQ004_P2_R5{sample}_1P.trimmed", sample=SAMPLES),
		p2u = expand("17NQ004_P2_R5{sample}_1U.trimmed", sample=SAMPLES)
	shell:
		"trimmomatic PE -threads 30 -phred33 {input.p1} {input.p2} {output.p1p} {output.p1u} {output.p2p} {output.p2u}"	


