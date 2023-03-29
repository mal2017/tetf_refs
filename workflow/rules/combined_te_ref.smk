rule combined_references:
    """
    note the rpm genome isn't gzipped
    """
    input:
        genome = "results/repeatmasker/genome.fasta.masked",
        transcripts = rules.get_flybase_resources.output.FULL_TRANSCRIPT_FASTA,
        tes = rules.bgzip_tes.output.fasta,
    conda:
        "../envs/te_gtf.yaml"
    output:
        genome = "results/combined-anno/genome-plus-tes.fasta.gz",
        transcripts = "results/combined-anno/transcripts-plus-tes.fasta.gz",
    shell:
        """
        # combine the genome and tes
        gunzip -c {input.tes} | cat {input.genome} - | bgzip -c > {output.genome}

        # combine the transcripts and tes
        cat {input.transcripts} {input.tes} | gunzip -c | bgzip -c > {output.transcripts}
        """

rule make_transcripts_and_consensus_tes_gtf:
    """
    A few small cleaning steps to make the host txome and TE gtf files compatible with
    salmon other downstream steps..s
    """
    input:
        host_gtf = rules.get_flybase_resources.output.TRANSCRIPTOME_GTF,
        te_gtf = rules.get_te_gtf.output.gtf,
    output:
        gtf = "results/combined-anno/transcripts-plus-tes.gtf",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.52.0--r41hd029910_0"
    script:
        "../scripts/make_transcripts_and_consensus_tes_gtf.R"




rule make_transcripts_and_consensus_tes_tx2gene:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        gtf = rules.make_transcripts_and_consensus_tes_gtf.output.gtf,
    output:
        tx2id = "results/combined-anno/transcripts-plus-tes.tx2id.tsv",
        tx2symbol = "results/combined-anno/transcripts-plus-tes.tx2symbol.tsv",
        tx2txsymbol = "results/combined-anno/transcripts-plus-tes.tx2txsymbol.tsv",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.52.0--r41hd029910_0"
    script:
        "../scripts/make_transcripts_and_consensus_tes_tx2gene.R"