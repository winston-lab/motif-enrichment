#!/usr/bin/env python

configfile: "config.yaml"

localrules: make_motif_database,
    get_peak_sequences,
    fimo,
    cat_fimo_results,
    intersect_fimo_results,
    test_motif_enrichment

ANNOTATIONS = config["annotations"]
CONDITIONS = [k for k,v in config["comparisons"].items()]
CONTROLS = [v["control"] for k,v in config["comparisons"].items()]

rule all:
    input:
        expand("data/{annotation}_allFIMOresults.tsv", annotation = ANNOTATIONS),
        expand("results/{condition}-v-{control}/{condition}-v-{control}_motif-enrichment.tsv", zip, condition=CONDITIONS, control=CONTROLS)

rule make_motif_database:
    input:
        motif_db = config["motif_databases"],
        fasta = config["genome"]["fasta"],
    output:
        "allmotifs.meme"
    log: "logs/make_motif_database.log"
    shell: """
        (meme2meme -bg <(fasta-get-markov {input.fasta}) {input.motif_db} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' > {output}) &> {log}
        """

# to avoid 'double-counting' regions that may result e.g. from
# a poorly called peak where one initiation event is called as
# two peaks, overlapping (but not book-ended) regions are merged
# and the qvalue of the most significant peak is kept
rule get_peak_sequences:
    input:
        annotation = lambda wc: ANNOTATIONS[wc.annotation]["path"],
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"],
    output:
        "data/{annotation}.fa"
    params:
        upstr = lambda wc: ANNOTATIONS[wc.annotation]["upstream"],
        dnstr = lambda wc: ANNOTATIONS[wc.annotation]["dnstream"],
    log: "logs/get_peak_sequences/get_peak_sequences_{annotation}.log"
    shell: """
        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.annotation} -g {input.chrsizes} | sort -k1,1 -k2,2n | bedtools merge -s -d -1 -c 4,5,6 -o collapse,max,first -i stdin | cut -f1-6 | bedtools getfasta -name+ -s -fi {input.fasta} -bed stdin > {output}) &> {log}
        """

rule fimo:
    input:
        fasta = "data/{annotation}.fa",
        motif_db = "allmotifs.meme"
    output:
        "data/{annotation}_allmotifs.bed"
    params:
        alpha = config["fimo-pval"],
        find_rc = [] if config["find-revcomp"] else "--norc"
    log: "logs/fimo/fimo_{annotation}.log"
    shell: """
        (fimo --bgfile <(fasta-get-markov {input.fasta}) --parse-genomic-coord {params.find_rc} --thresh {params.alpha} --text {input.motif_db} {input.fasta} | awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $3, $4, $5+1, $1, -log($8)/log(10), $6, $2, $10}} ' > {output}) &> {log}
        """

# NOTE: for a motif to match a region, I require that the entire motif lies within the region.
# NOTE: the results report every overlap if there are multiple matches in a region.
#       In the enrichment testing step, region matches are made binary
# TODO: the making of the regions should really be made into a script to ensure that the
#       regions in 'get_peak_sequences' are really the same as the regions being intersected
rule intersect_fimo_results:
    input:
        annotation = lambda wc: ANNOTATIONS[wc.annotation]["path"],
        results = "data/{annotation}_allmotifs.bed",
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"],
    output:
        "data/{annotation}_allFIMOresults.tsv",
    params:
        upstr = lambda wc: ANNOTATIONS[wc.annotation]["upstream"],
        dnstr = lambda wc: ANNOTATIONS[wc.annotation]["dnstream"],
        stranded = [] if config["find-revcomp"] else "-s"
    log: "logs/intersect_fimo_results/intersect_fimo_results-{annotation}.log"
    shell: """
        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.annotation} -g {input.chrsizes} |  sort -k1,1 -k2,2n | bedtools merge -s -d -1 -c 4,5,6 -o collapse,max,first -i stdin | cut -f1-6 | bedtools intersect -a stdin -b {input.results} {params.stranded} -F 1 -wao | cut -f15 --complement | cat <(echo -e "chrom\tregion_start\tregion_end\tregion_id\tregion_score\tregion_strand\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "data/{condition}_allFIMOresults.tsv",
        fimo_neg = "data/{control}_allFIMOresults.tsv",
    output:
        tsv = "results/{condition}-v-{control}/{condition}-v-{control}_motif-enrichment.tsv",
        heatmap = "results/{condition}-v-{control}/{condition}-v-{control}_motif-heatmaps.svg",
        meta = "results/{condition}-v-{control}/{condition}-v-{control}_motif-metagenes.svg",
        plot = "results/{condition}-v-{control}/{condition}-v-{control}_motif-enrichment.svg"
    params:
        fimo_pval = config["fimo-pval"],
        fdr_cutoff = config["enrichment-fdr"],
        upstream = lambda wc: config["comparisons"][wc.condition]["upstream"]
    script: "scripts/motif_enrichment.R"

