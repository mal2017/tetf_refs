rule get_tidal:
    output:
        temp("results/tmp_tidal_anno.tar.gz"),
        multiext("results/tidal-anno/","dm6.chr.len","fly_virus_structure_repbase.fa","gem_mappability_dm6_100mer.mappability","refflat_dm6.txt","repmasker_dm6_track.txt","Tidalbase_Dmel_TE_classifications_2015.txt","Tidalbase_transposon_sequence.fa")
    params:
        uri = config.get("TIDAL_ANNOTATION"),
    shell:
        """
        mkdir -p results &&
        curl -LJ {params.uri} > {output[0]} &&
        tar -xzf {output[0]} -C results/ &&

        mv results/annotation/* $(dirname {output[1]}) &&
        rm -r results/annotation/
        """


rule filter_for_dmel_only_tes:
    """
    For Tidalbase TEs, non-*Dmel* seqs can be removed with the `seqkit` command below. 
    This works for Tidalbase because *Dmel* TEs don't have 'Dmel' in the name, but other TEs have the
    4 letter species abbreviation.
    """
    input:
        "results/tidal-anno/Tidalbase_transposon_sequence.fa"
    output:
        "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fa"
    singularity:
        "docker://quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0"
    shell:
        """
        seqkit grep -r -n -v -p '\\\\|D.{{3}}\\\\\\\\' {input} > {output}
        """

rule fix_te_names:
    input:
        fasta = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fa"
    output:
        fasta = temp("results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.fa")
    script:
        "../scripts/fix_te_names.R"

rule bgzip_tes:
    input:
        fasta = rules.fix_te_names.output.fasta
    output:
        fasta = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.fa.gz"
    conda:
        "../envs/te_gtf.yaml"
    shell:
        """
        bgzip -f -c {input.fasta} > {output.fasta}
        """


rule get_te_gtf:
    input:
        fasta = rules.bgzip_tes.output.fasta
    output:
        fai = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.fa.gz.fai",
        bed = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.consensus.bed",
        gp = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.consensus.gp",
        gtf = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.consensus.gtf.gz"
    conda:
        "../envs/te_gtf.yaml"
    shell:
        """
        samtools faidx {input.fasta} &&

        bedtools makewindows -g {output.fai} -n 1 -i src > {output.bed} &&

        bedToGenePred {output.bed} stdout | awk 'BEGIN{{FS=OFS="\t"}} {{$3="+"}} 1' > {output.gp}

        genePredToGtf file {output.gp} stdout -source='TIDAL' | gffread -T | grep 'exon' | tr 'exon' 'mRNA' | gzip > {output.gtf}
        """

