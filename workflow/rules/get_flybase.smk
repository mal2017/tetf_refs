rule get_flybase_resources:
    params:
        GENOME_FASTA = config.get("GENOME_FASTA"),
        FULL_TRANSCRIPT_FASTA = config.get("FULL_TRANSCRIPT_FASTA"),
        TRANSCRIPTOME_GTF = config.get("TRANSCRIPTOME_GTF"),
        MISCRNA_FASTA = config.get("MISCRNA_FASTA"),
        NCRNA_FASTA = config.get("NCRNA_FASTA"),
        TRNA_FASTA = config.get("TRNA_FASTA"),
    output:
        GENOME_FASTA = "results/flybase-anno/genome.fasta.gz",
        FULL_TRANSCRIPT_FASTA = "results/flybase-anno/transcripts.fasta.gz",
        TRANSCRIPTOME_GTF= "results/flybase-anno/transcriptome.gtf.gz",
        MISCRNA_FASTA = "results/flybase-anno/miscRNA.fasta.gz",
        NCRNA_FASTA = "results/flybase-anno/ncRNA.fasta.gz",
        TRNA_FASTA = "results/flybase-anno/tRNA.fasta.gz",
    shell:
        """
        curl {params.GENOME_FASTA} > {output.GENOME_FASTA}
        curl {params.FULL_TRANSCRIPT_FASTA} > {output.FULL_TRANSCRIPT_FASTA}
        curl {params.TRANSCRIPTOME_GTF} > {output.TRANSCRIPTOME_GTF}
        curl {params.MISCRNA_FASTA} > {output.MISCRNA_FASTA}
        curl {params.NCRNA_FASTA} > {output.NCRNA_FASTA}
        curl {params.TRNA_FASTA} > {output.TRNA_FASTA}
        """
