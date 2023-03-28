
rule dna_bwa_mem2_index:
    input:
        rules.combined_references.output.genome,
    output:
        expand("results/bwa_mem2-idx/idx.{suf}",suf=["0123","amb","ann","bwt.2bit.64","pac"])
    params:
        pfx = "results/bwa_mem2-idx/idx"
    threads:
        24
    resources:
        time=60,
        mem=48000,
        cpus=24
    singularity:
        "docker://quay.io/biocontainers/bwa-mem2:2.2.1--h9a82719_1"
    shell:
        """
        bwa-mem2 index -p {params.pfx} {input}
        """