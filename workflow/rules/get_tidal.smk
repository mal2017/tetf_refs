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


# BEGIN ANTISENSE CHECKING RULES
# workflow to doublecheck that the TIDAL/FLYBASE tes have the appropriate sense
# i observed that a few TE consensus seqs are antisense relative to the equivalent DFAM annotation.
# relevance is most clear for autonomous TEs
rule get_repeat_peps:
    params:
        peps = config.get("REPEAT_PEPS")
    output:
        "results/test_repeat_strandedness/repeat_peps.faa"
    shell:
        """
        curl -LJ {params.peps} > {output}
        """

rule align_repeat_peps_to_tidal_repeats:
    input:
        reps = rules.fix_te_names.output.fasta,
        peps = rules.get_repeat_peps.output
    output:
        "results/test_repeat_strandedness/repeat_peps_aligned_to_tidal_repeats.gff"
    singularity:
        "docker://quay.io/biocontainers/miniprot:0.9--h7132678_0"
    shell:
        """
        miniprot --gff {input.reps} {input.peps} | grep -v "##PAF" | grep "mRNA" > {output}
        """

rule get_antisensed_tes:
    input:
        gff = rules.align_repeat_peps_to_tidal_repeats.output,
    output:
        "results/test_repeat_strandedness/antisensed_tes.txt"
    shell:
        """
        cat {input.gff} | awk '$7=="-"{{print $1}}' | sort | uniq > {output}
        """

rule fix_antisensed_tes:
    """
    takes TEs that don't match the patterns for antisensed tes and just puts them into the output.
    them, takes the tes that do match the patterns for antisensed TEs and reverses and complements (-p and -r)
    them and concatenates them to the output.
    checked by aligning copia (negcontrol) from before and after the antisense fixing step and accord2 (poscontrol) from
    before and after the antisense fixing step to themselves. Inspection of the dotplots showed that accord2 flipped and copia didn't.
    """
    input:
        fa = rules.fix_te_names.output.fasta,
        antisensed = rules.get_antisensed_tes.output
    output:
        fasta = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa"
    singularity:
        "docker://quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0"
    shell:
        """
        seqkit grep -v -f {input.antisensed} {input.fa} > {output.fasta} &&

        seqkit grep -f {input.antisensed} {input.fa} | seqkit seq -p -r >> {output.fasta}
        """
# END ANTISENSE CHECKING RULES

rule bgzip_tes:
    input:
        fasta = rules.fix_antisensed_tes.output.fasta
    output:
        fasta = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz"
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
        fai = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.fa.gz.fai",
        bed = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.consensus.bed",
        gp = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.consensus.gp",
        gtf = "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.consensus.gtf.gz"
    conda:
        "../envs/te_gtf.yaml"
    shell:
        """
        samtools faidx {input.fasta} &&

        bedtools makewindows -g {output.fai} -n 1 -i src > {output.bed} &&

        bedToGenePred {output.bed} stdout | awk 'BEGIN{{FS=OFS="\t"}} {{$3="+"}} 1' > {output.gp}

        genePredToGtf file {output.gp} stdout -source='TIDAL' | gffread -T | grep 'exon' | bgzip -c > {output.gtf}
        """

