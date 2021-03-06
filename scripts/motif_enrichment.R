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
    # 0. get total number of unique regions per annotation based on unique coordinates
    # 1. filter out regions matching no motifs
    # 2. group by motif id
    # 3. get number of unique regions matched per motif based on unique coordinates
    df %<>%
        select(chrom, region_start, region_end, region_id, region_strand, region_score,
               motif_start, motif_id, motif_alt_id, annotation) %>%
        group_by(annotation) %>%
        mutate(n_regions = n_distinct(chrom, region_start, region_end, region_strand)) %>%
        filter(motif_start >= 0) %>%
        group_by(annotation, motif_id, motif_alt_id, n_regions) %>%
        summarise(n_matches = n_distinct(chrom, region_start, region_end, region_strand)) %>%
        ungroup()

    fisher_df = df %>% filter(annotation==condition) %>%
        transmute(motif_id, motif_alt_id,
                  condition_withmotif = n_matches,
                  condition_nomotif = n_regions-n_matches) %>%
        left_join(df %>% filter(annotation==control) %>%
                      transmute(motif_id, motif_alt_id,
                                control_withmotif = n_matches,
                                control_nomotif = n_regions-n_matches),
                  by=c("motif_id", "motif_alt_id")) %>%
        filter(! is.na(condition_nomotif) & ! is.na(condition_withmotif) &
                   ! is.na(control_nomotif) & ! is.na(control_withmotif))

    if (nrow(fisher_df) > 0){
        fisher_df %<>%
            group_by(motif_id, motif_alt_id,
                     condition_withmotif, condition_nomotif, control_withmotif, control_nomotif) %>%
            do(fisher.test(matrix(c(.$condition_withmotif, .$condition_nomotif,
                                    .$control_withmotif, .$control_nomotif), nrow=2, ncol=2, byrow=TRUE),
                           alternative="two.sided") %>% tidy()) %>%
            ungroup() %>%
            mutate(log10_fdr = -log10(p.adjust(p.value, method="BH")),
                   condition_fraction = condition_withmotif/(condition_withmotif+condition_nomotif),
                   control_fraction = control_withmotif/(control_withmotif+control_nomotif)) %>%
            mutate_at(vars(estimate, conf.low, conf.high), funs(log2(.))) %>%
            arrange(desc(log10_fdr), p.value, desc(estimate), desc(condition_withmotif)) %>%
            select(motif_id, motif_alt_id, log10_fdr, log2_odds_ratio=estimate,
                   conf_low=conf.low, conf_high=conf.high,
                   condition_withmotif, condition_nomotif, condition_fraction,
                   control_withmotif, control_nomotif, control_fraction) %>%
            mutate_if(is_double, funs(signif(., digits=3))) %>%
            write_tsv(out_path) %>%
            mutate(label=if_else(is.na(motif_alt_id), motif_id, motif_alt_id))
    } else {
        fisher_df = tibble(motif_id = character(),
                           motif_alt_id = character(),
                           log10_fdr = double(),
                           log2_odds_ratio = double(),
                           conf_log = double(),
                           conf_high = double(),
                           condition_withmotif = integer(),
                           condition_nomotif = integer(),
                           condition_fraction = double(),
                           control_withmotif = integer(),
                           control_nomotif = integer(),
                           control_fraction = double(),
                           label = character()) %>%
        write_tsv(out_path)
    }

    plot = ggplot() +
        geom_vline(xintercept=0, color="grey65") +
        geom_hline(yintercept=fdr_cutoff, linetype="dashed") +
        geom_point(data=fisher_df %>%
                        filter(log10_fdr <= fdr_cutoff),
                   aes(x=log2_odds_ratio, y=log10_fdr),
                   shape=16, size=1, alpha=0.4, stroke=0) +
        geom_label_repel(data = fisher_df %>%
                            filter(log10_fdr > fdr_cutoff),
                         aes(x=log2_odds_ratio, y=log10_fdr, label=label,
                             fill=if_else(log2_odds_ratio<0, "depleted", "enriched")),
                             size=8/72*25.4,
                             box.padding=unit(0, "pt"),
                             label.r=unit(0.5, "pt"),
                             label.size=NA,
                             label.padding = unit(0.8, "pt"),
                             ylim = c(fdr_cutoff, NA),
                             force = 0.5,
                             segment.size=0.2,
                             segment.alpha=0.5) +
        xlab(expression(bold(log[2] ~ "odds-ratio"))) +
        ylab(expression(bold(-log[10] ~ "FDR"))) +
        scale_fill_manual(values = c("#CDE7FD", "#FFCEC8")) +
        ggtitle(paste0("motif enrichment upstream of\n", condition, " vs. ", control),
                subtitle="Fisher's exact test (two-tailed)") +
        theme_light() +
        theme(text = element_text(size=12, face="bold", color="black"),
              axis.text = element_text(size=10, color="black"),
              axis.title = element_text(face="bold"),
              plot.subtitle = element_text(face="plain"),
              legend.position="none")
    ggsave(out_plot, plot=plot, width=30, height=17, units="cm")
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

