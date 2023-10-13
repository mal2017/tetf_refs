rule salmon_decoy_mask_exons:
    """
    We generated a partial decoy-aware transcriptome reference for salmon as described
    in the [salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).
    Briefly, we began by extracting exon genomic coordinates and masking them with bedtools
    v2.30.0 maskfasta.
    """
    input:
        gtf = rules.get_flybase_resources.output.TRANSCRIPTOME_GTF,
        genome = rules.get_flybase_resources.output.GENOME_FASTA,
    output:
        exons = "results/salmon-idx/decoy/exons.bed",
        genome = "results/salmon-idx/decoy/genome.fasta",
        masked = "results/salmon-idx/decoy/txome.masked.genome.fasta"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    priority: 50
    shell:
        """
        zcat {input.gtf} | awk -v OFS='\\t' '{{if ($3=="exon") {{print $1,$4,$5}}}}' > {output.exons} &&
        gunzip -c {input.genome} > {output.genome} &&
        bedtools maskfasta -fi {output.genome} -bed {output.exons} -fo {output.masked}
        """

rule salmon_decoy_mashmap_align_txome:
    """
    Next, we applied mashmap v2.0 with options "--pi 80 -s 500" to identify non-genic regions with similarity
    to the transcriptome.
    """
    input:
        masked = rules.salmon_decoy_mask_exons.output.masked,
        txpfile = rules.get_flybase_resources.output.FULL_TRANSCRIPT_FASTA,
    output:
        mm = "results/salmon-idx/decoy/mashmap.out",
        txlike = "results/salmon-idx/decoy/genome_found.sorted.bed",
    threads:
        8
    singularity:
        "docker://quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    log:
        "results/logs/mashmap/transcripts_and_consensus_tes.txt"
    resources:
        time=20,
        mem=15000,
        cpus=8
    priority: 50
    shell:
        """
        mashmap -r {input.masked} -q {input.txpfile} -t {threads} --pi 80 -s 500 -o {output.mm} 2> {log} &&
        awk -v OFS='\\t' '{{print $6,$8,$9}}' {output.mm} | sort -k1,1 -k2,2n - > {output.txlike}
        """

rule salmon_decoy_finalize:
    """
    We extracted the sequence of the transcriptome-similar regions with bedtools getfasta.
    Finally, we concatenated the custom TE and transcriptome reference sequences with the decoy
    sequences to generate the decoy-aware transcriptome reference which we provided to salmon index.
    """
    input:
        txlike = rules.salmon_decoy_mashmap_align_txome.output.txlike,
        masked = rules.salmon_decoy_mask_exons.output.masked,
        txpfile = rules.combined_references.output.transcripts,
        miscrna = rules.get_flybase_resources.output.MISCRNA_FASTA,
        trna = rules.get_flybase_resources.output.TRNA_FASTA,
    output:
        txpfile = "results/salmon-idx/decoy/txpfile.fasta",
        txlike_merged = "results/salmon-idx/decoy/genome_found_merged.bed",
        txlike_fa = "results/salmon-idx/decoy/genome_found.fasta",
        decoy = "results/salmon-idx/decoy/decoy.fasta",
        gentrome = "results/salmon-idx/decoy/gentrome.fasta",
        masked_fai = "results/salmon-idx/decoy/txome.masked.genome.fasta.fai"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    priority: 50
    shell:
        """
        gunzip -c {input.txpfile} > {output.txpfile} &&
        bedtools merge -i {input.txlike} > {output.txlike_merged} &&
        bedtools getfasta -fi {input.masked} -bed {output.txlike_merged} -fo {output.txlike_fa} &&
        awk '{{a=$0; getline;split(a, b, ":");  r[b[1]] = r[b[1]]""$0}} END {{ for (k in r) {{ print k"\\n"r[k] }} }}' {output.txlike_fa} > {output.decoy} &&
        zcat {input.miscrna} {input.trna} | cat - {output.txpfile} {output.decoy} > {output.gentrome}
        """

rule salmon_decoy_get_ids:
    input:
        decoy = rules.salmon_decoy_finalize.output.decoy
    output:
        ids = "results/salmon-idx/decoy/decoy.txt"
    priority: 50
    shell:
        """
        grep ">" {input.decoy} | awk '{{print substr($1,2); }}' > {output.ids}
        """

rule salmon_index_transcripts_and_consensus_tes:
    """
    We indexed the custom transcriptome with salmon v1.5.2 index using "-k <insert k here>."
    """
    input:
        fa = rules.salmon_decoy_finalize.output.gentrome,
        decoyids = rules.salmon_decoy_get_ids.output.ids,
    output:
        directory("results/salmon-idx/index/")
    params:
        k = config.get("SALMON_K_PARAM"),
        dec = "--decoys results/salmon-idx/decoy/decoy.txt" 
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    threads:
        8
    log:
        "results/logs/salmon_index/transcripts_and_consensus_tes.txt"
    resources:
        time=60,
        mem=20000,
        cpus=8
    priority: 50
    shell:
        """
        salmon index -t {input.fa} \
            --index {output} \
            -k {params.k} \
            -p {threads} \
            {params.dec} 2> {log}
        """

rule make_salmon_te_aux_target_file:
    """
    TEs, miscRNAs, and tRNAs were provided as auxiliary targets to salmon to avoid applying bias correction
    to TE expression estimates.
    """
    input:
        tes = rules.bgzip_tes.output.fasta,
        miscrna = rules.get_flybase_resources.output.MISCRNA_FASTA,
        trna = rules.get_flybase_resources.output.TRNA_FASTA,
    output:
        "results/salmon-idx/transcripts_and_consensus_tes.aux.txt"
    priority: 50
    shell:
        """
        zcat {input.tes} | grep ">" | tr -d ">" > {output} &&
        zcat {input.miscrna} {input.trna} | grep ">" | cut -f 1 -d " " | tr -d ">" >> {output}
        """