#!/usr/bin/env python

from math import log10

configfile: "config.yaml"

localrules: make_motif_database,

COMPARISONS = config["comparisons"]

#get all motif names from motif databases, cleaning nasty characters in some motif names
MOTIFS = set(subprocess.run(args="meme2meme " + " ".join(config["motif_databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split())

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        "motifs/allmotifs.bed",
        expand("comparisons/{comparison}/{comparison}_{group}_unmergedFIMOresults.tsv.gz", comparison=COMPARISONS, group=["condition", "control"]) if COMPARISONS else [],
        expand("comparisons/{comparison}/{comparison}_motif-enrichment.svg", comparison=COMPARISONS) if COMPARISONS else []

rule make_motif_database:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = config["motif_databases"],
    output:
        "motifs/allmotifs.meme"
    log:
        "logs/make_motif_database.log"
    shell: """
        (meme2meme -bg <(fasta-get-markov {input.fasta}) {input.motif_db} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' > {output}) &> {log}
        """

#run fimo in parallel for each motif
rule fimo:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = "motifs/allmotifs.meme"
    output:
        bed = temp("motifs/.{motif}.bed") # a BED6+2 format
    params:
        alpha = config["fimo-pval"],
        find_rc = [] if config["find-revcomp"] else "--norc"
    log:
        "logs/fimo/fimo_{motif}.log"
    shell: """
        (fimo --motif {wildcards.motif} --bgfile <(fasta-get-markov {input.fasta}) {params.find_rc} --thresh {params.alpha} --text {input.motif_db} {input.fasta} | awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $3, $4-1, $5, $1, -log($8)/log(10), $6, $2, $10}}' > {output.bed}) &> {log}
        """

rule cat_fimo_motifs:
    input:
        bed = expand("motifs/.{motif}.bed", motif=MOTIFS)
    output:
        bed = "motifs/allmotifs.bed",
    threads: config["threads"]
    shell: """
        cat {input.bed} | sort -k1,1 -k2,2n --parallel={threads} > {output.bed}
        """

rule get_motifs_unmerged:
    input:
        annotation = lambda wc: COMPARISONS[wc.comparison][wc.group]["path"],
        motifs = "motifs/allmotifs.bed",
        fasta = config["genome"]["fasta"]
    output:
        "comparisons/{comparison}/{comparison}_{group}_unmergedFIMOresults.tsv.gz"
    params:
        upstr = lambda wc: COMPARISONS[wc.comparison][wc.group]["upstream"],
        dnstr = lambda wc: COMPARISONS[wc.comparison][wc.group]["dnstream"],
    log:
        "logs/get_overlapping_motifs/get_overlapping_motifs-{comparison}-{group}.log"
    shell: """
        (cut -f1-6 {input.annotation} | \
         bedtools slop -l {params.upstr} -r {params.dnstr} -s -i stdin -g <(faidx {input.fasta} -i chromsizes) | \
         sort -k1,1 -k2,2n | \
         bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
         cut -f15 --complement | \
         cat <(echo -e "chrom\tregion_start\tregion_end\tregion_id\tregion_score\tregion_strand\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
         pigz -f > {output}) &> {log}
        """

#bedtools intersect regions with fimo motifs
#0. with the region as reference, extend annotation to upstream and 'downstream' distances
#1. merge overlapping (but not book-ended) features
#2. intersect with motif file
rule get_motifs_merged:
    input:
        annotation = lambda wc: COMPARISONS[wc.comparison][wc.group]["path"],
        motifs = "motifs/allmotifs.bed",
        fasta = config["genome"]["fasta"]
    output:
        "comparisons/{comparison}/{comparison}_{group}_mergedFIMOresults.tsv.gz"
    params:
        upstr = lambda wc: COMPARISONS[wc.comparison][wc.group]["upstream"],
        dnstr = lambda wc: COMPARISONS[wc.comparison][wc.group]["dnstream"],
    log:
        "logs/get_overlapping_motifs/get_overlapping_motifs-{comparison}-{group}.log"
    shell: """
        (cut -f1-6 {input.annotation} | \
         bedtools slop -l {params.upstr} -r {params.dnstr} -s -i stdin -g <(faidx {input.fasta} -i chromsizes) | \
         sort -k1,1 -k2,2n | \
         bedtools merge -s -d -1 -c 4,5,6 -o collapse,max,first -i stdin | \
         sort -k1,1 -k2,2n | \
         bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
         cut -f15 --complement | \
         cat <(echo -e "chrom\tregion_start\tregion_end\tregion_id\tregion_score\tregion_strand\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
         pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "comparisons/{comparison}/{comparison}_condition_mergedFIMOresults.tsv.gz",
        fimo_neg = "comparisons/{comparison}/{comparison}_control_mergedFIMOresults.tsv.gz",
    output:
        tsv = "comparisons/{comparison}/{comparison}_motif-enrichment.tsv",
        plot = "comparisons/{comparison}/{comparison}_motif-enrichment.svg",
        # heatmap = "results/{condition}-v-{control}/{condition}-v-{control}_motif-heatmaps.svg",
        # meta = "results/{condition}-v-{control}/{condition}-v-{control}_motif-metagenes.svg",
    params:
        fdr_cutoff = -log10(config["enrichment-fdr"]),
        cond_label = lambda wc: COMPARISONS[wc.comparison]["condition"]["label"],
        ctrl_label = lambda wc: COMPARISONS[wc.comparison]["control"]["label"],
        # upstream = lambda wc: config["comparisons"][wc.condition]["upstream"]
    conda:
        "envs/tidyverse.yaml"
    script:
        "scripts/motif_enrichment.R"

