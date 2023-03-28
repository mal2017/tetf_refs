
rule bt1_unmasked_index:
    """
    Currently used only for CSEM/MOSAICS pipeline, where we want to align to the unmasked genome.
    """
    input:
        ref =rules.get_flybase_resources.output.GENOME_FASTA
    output:
        directory('results/bt1_unmasked-idx/')
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/bowtie:1.3.1--py310h4070885_4"
    shell:
        """
        mkdir -p {output} &&
        cp {input.ref} {output} &&
        bowtie-build --threads {threads} {input.ref} {output}/idx
        """