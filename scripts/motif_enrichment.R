library(tidyverse)
library(forcats)
library(magrittr)
library(broom)
library(ggrepel)
library(ggpmisc)

import = function(path, annotation_id){
    read_tsv(path, col_types = 'ciicdcciicdccc') %>%
        mutate(start = case_when(region_strand=="+" & motif_chrom != "." ~ region_start-motif_start,
                                 region_strand=="-" & motif_chrom != "." ~ motif_start-region_end,
                                 motif_chrom=="." ~ as.integer(NaN)),
               end = case_when(region_strand=="+" & motif_chrom != "." ~ region_start-motif_end,
                               region_strand=="-" & motif_chrom != "." ~ motif_end-region_end,
                               motif_chrom=="." ~ as.integer(NaN)),
               annotation = annotation_id) %>%
        arrange(region_score) %>%
        return()
}

main = function(fdr_cutoff, condition, control,
                in_fimo_pos, in_fimo_neg, upstream,
                out_path, out_plot
                # out_meta, out_heatmap
                ){

    df = in_fimo_pos %>%
        import(annotation_id = condition) %>%
        bind_rows(in_fimo_neg %>% import(annotation_id = control)) %>%
        # mutate(index = paste(chrom, region_start, region_end, region_strand) %>%
        #            fct_inorder(ordered=TRUE) %>%
        #            as.integer())
        mutate(annotation = fct_inorder(annotation, ordered=TRUE))

    ##TODO: facet_grid paginate for multiple motifs
    #heatmap = ggplot(data = df %>% filter(! is.na(start)),
    #       aes(x=start, xend=end, y=index, yend=index, color=motif_logpval)) +
    #    geom_vline(xintercept = 0) +
    #    geom_segment(na.rm=TRUE, size=0.3, alpha=0.8) +
    #    scale_color_distiller(palette="Blues", direction=1,
    #                          name = bquote(bold(-log[10] ( "motif p-value"))),
    #                          guide=guide_colorbar(barwidth=8, barheight=0.3,
    #                                               title.position = "top")) +
    #    scale_x_continuous(limits = c(-upstream, NA),
    #                       labels = function(x) if_else(x==0, "region end", as.character(x)),
    #                       name = NULL) +
    #    facet_grid(annotation~motif_id, scales="free_y", space="free_y", switch="y") +
    #    theme_light() +
    #    theme(text = element_text(size=12, color="black", face="bold"),
    #          axis.text = element_text(size=10, color="black"),
    #          axis.text.y = element_blank(),
    #          axis.title = element_blank(),
    #          strip.placement = "outside",
    #          strip.background = element_blank(),
    #          strip.text = element_text(color="black"),
    #          strip.text.y = element_text(angle=-180, hjust=1, vjust=0.5),
    #          legend.position = "top",
    #          legend.text = element_text(size=10, face="plain"))
    #ggsave(out_heatmap, plot=heatmap, width=14, height=20, units="cm")

    #meta_df = df %>% filter(! is.na(start)) %>%
    #    transmute(annotation = annotation,
    #              motif_id = motif_id,
    #              position = (end+start)/2)

    #meta = ggplot(data = meta_df,
    #       aes(x=position, fill=annotation, color=annotation)) +
    #    geom_density(aes(y=..density..),
    #                 bw=4, na.rm=TRUE, alpha=0.6) +
    #    scale_x_continuous(limits = c(-upstream, NA),
    #                       expand = c(0,0)) +
    #    facet_wrap(~motif_id) +
    #    theme_light() +
    #    theme(text = element_text(size=12, color="black", face="bold"),
    #          axis.text = element_text(size=10, color="black"),
    #          axis.title.x = element_blank(),
    #          strip.placement = "outside",
    #          strip.background = element_blank(),
    #          strip.text = element_text(color="black"),
    #          legend.position = "bottom",
    #          legend.text = element_text(size=10, face="plain"))
    #ggsave(out_meta, plot=meta, width=12, height=12, units="cm")

    #TODO: possible to add in logistic regression based on the region scores using
    # this df before counting
    df %<>% select(region_start, region_end, region_id, region_strand, region_score,
                  motif_start, motif_id, motif_alt_id, annotation) %>%
        complete(nesting(annotation, region_start, region_end, region_id, region_strand, region_score),
                 nesting(motif_id, motif_alt_id)) %>%
        filter(motif_id != ".") %>%
        group_by(annotation, region_start, region_end, region_id, region_strand, region_score,
                 motif_id, motif_alt_id) %>%
        summarise(match = if_else(all(! is.na(motif_start)), TRUE, FALSE)) %>%
        ungroup()

    fisher_df = df %>%
        filter(annotation==condition) %>%
        count(motif_id, motif_alt_id, match) %>%
        spread(match, n) %>%
        magrittr::set_colnames(c("motif_id", "motif_alt_id", "condition_nomotif", "condition_withmotif")) %>%
        left_join(df %>% filter(annotation==control) %>%
                      count(motif_id, motif_alt_id, match) %>%
                      spread(match, n) %>%
                      magrittr::set_colnames(c("motif_id", "motif_alt_id", "control_nomotif", "control_withmotif")),
                  by=c("motif_id", "motif_alt_id")) %>%
        group_by(motif_id, motif_alt_id,
                 condition_withmotif, condition_nomotif, control_withmotif, control_nomotif) %>%
        do(fisher.test(matrix(c(.$condition_withmotif, .$condition_nomotif,
                                .$control_withmotif, .$control_nomotif), nrow=2, ncol=2, byrow=TRUE),
                       alternative="two.sided") %>% tidy()) %>%
        ungroup() %>%
        mutate(fdr = p.adjust(p.value, method="BH")) %>%
        mutate_at(vars(estimate, conf.low, conf.high), funs(log2(.))) %>%
        arrange(fdr, p.value, desc(estimate), desc(condition_withmotif)) %>%
        select(motif_id, motif_alt_id, fdr, log2_odds_ratio=estimate,
               conf_low=conf.low, conf_high=conf.high,
               condition_withmotif, condition_nomotif, control_withmotif, control_nomotif) %>%
        mutate_if(is_double, funs(signif(., digits=3))) %>%
        write_tsv(out_path) %>%
        mutate(label=if_else(is.na(motif_alt_id), motif_id, motif_alt_id))

    plot = ggplot() +
        geom_vline(xintercept=0, color="grey65") +
        geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed") +
        geom_point(data=fisher_df, aes(x=log2_odds_ratio, y=-log10(fdr)),
                   shape=16, size=1, alpha=0.8, stroke=0) +
        xlab(expression(bold(paste(log[2], " odds-ratio")))) +
        ylab(expression(bold(paste(-log[10], " FDR")))) +
        ggtitle(paste0("motif enrichment upstream of\n", condition, " vs. ", control),
                subtitle="Fisher's exact test (two-tailed)") +
        theme_light() +
        theme(text = element_text(size=12, face="bold", color="black"),
              axis.text = element_text(size=10, color="black"),
              plot.subtitle = element_text(face="plain"))
    if (nrow(fisher_df %>% filter(fdr<fdr_cutoff)) > 10){
        plot = plot +
            stat_dens2d_labels(data = fisher_df %>% filter(fdr<fdr_cutoff),
                               aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                               geom="text_repel", keep.number=25, size=10/72*25.4)
    } else {
        plot = plot +
            geom_text_repel(data = fisher_df %>% filter(fdr<fdr_cutoff),
                            aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                            size=10/72*25.4)
    }
    ggsave(out_plot, plot=plot, width=14, height=12, units="cm")
}

main(fdr_cutoff = snakemake@params[["fdr_cutoff"]],
     condition= snakemake@params[["cond_label"]],
     control= snakemake@params[["ctrl_label"]],
     in_fimo_pos = snakemake@input[["fimo_pos"]],
     in_fimo_neg = snakemake@input[["fimo_neg"]],
     # upstream = snakemake@params[["upstream"]],
     out_path = snakemake@output[["tsv"]],
     # out_heatmap = snakemake@output[["heatmap"]],
     # out_meta = snakemake@output[["meta"]],
     out_plot = snakemake@output[["plot"]])
