library(tidyverse)
library(RColorBrewer)
library(phyloseq)

get_cumsums <- function(ps_obj){
  css <- cumsum(sort(taxa_sums(ps_obj)/sum(otu_table(ps_obj))*100,
                     decreasing = T))
  return(css)
}

get_cumabd <- function(ps_obj){
  cabd <- sort(taxa_sums(ps_obj)/nsamples(ps_obj),
               decreasing = T)
  return(cabd)
}


plot_bar_top <- function(ps_obj.comp, taxa = 'Phylum',
                         order_list = NA,
                         retall = FALSE,
                         x_label = NA,
                         blank_axis = F,
                         cutoff = NA,
                         facet = NA){
  if (is.na(x_label)) {
    x_label = 'Sample'
  }
  cutoffs <- c(Phylum=100,
               Class=100,
               Order=100,
               Family=99,
               Genus=90,
               ASV=100)
  if (is.na(cutoff)){
    cutoff <- cutoffs[taxa]
  }
  plot_name <- paste0('barplot_', cutoff, '_', taxa,
                      '_', x_label, '.png')
  csv_name <- paste0('barplot_', cutoff, '_', taxa,
                      '_', x_label, '_data.csv')
  ps_obj.comp.glom <- tax_glom(ps_obj.comp, taxa, NArm = F)
  cumsums <- get_cumsums(ps_obj.comp.glom)
  if (cutoff < 100) {
    taxa_rename <- names(which(cumsums > cutoff)[-1])
    # hline <- cumsums[which(cumsums > 95)[1]]
    renamed_taxa_value <- paste0('Other < ', 100 - cutoff,'%')
    # cumabds <- get_cumabd(ps_obj.comp.glom)
    # taxa_rename <- names(which(cumabds < 1.00))
    renamed_taxa <- tax_table(ps_obj.comp.glom)[taxa_rename, taxa]
    tax_table(ps_obj.comp.glom)[taxa_rename, taxa] <- renamed_taxa_value
    ps_obj.comp.glom <- tax_glom(ps_obj.comp.glom, taxa, NArm = F)
  }
  taxa_count <- n_distinct(tax_table(ps_obj.comp.glom)[,taxa])
  # col_pal <- pomological_palette
  col_pal <- cb_pal.12
  top_abd <- get_cumabd(ps_obj.comp.glom)
  taxa_order <- as.character(tax_table(ps_obj.comp.glom)[names(rev(top_abd)), taxa])
  taxa_order <- unique(taxa_order)
  if (cutoff < 100) {
    taxa_order <- taxa_order[taxa_order != renamed_taxa_value]
    taxa_order <- c(renamed_taxa_value, taxa_order)
  }
  if (any(is.na(taxa_order))) {
    taxa_order <- na.omit(taxa_order)
    taxa_order <- c(NA, taxa_order)
    taxa_count <- taxa_count - 1
  }
  if (taxa_count > length(col_pal)) {
    getPalette <- colorRampPalette(col_pal)
    bar_colors <- getPalette(taxa_count)
  } else {
    bar_colors <- col_pal[1:taxa_count]
  }
  if (length(order_list) > 1) {
    sample_order <- order_list
  } else {
    sample_order <- sample_names(ps_obj.comp.glom)
  }
  top_bar <- plot_bar(ps_obj.comp.glom, fill = taxa) +
    geom_bar(aes(color=top_bar$data[,taxa]),
             stat='identity',
             position='stack') +
    scale_color_manual(name = taxa,
                       values = rev(bar_colors),
                       aesthetics = c('color', 'fill'),
                       na.value = "#000000",
                       breaks = taxa_order) +
    ylim(c(0,101)) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  top_bar <- top_bar + xlab(x_label)
  if (blank_axis) {
    top_bar <- top_bar + theme(axis.text.x = element_blank())
  }
  # if (taxa == 'Genus') {
  #   top_bar <- top_bar + theme(legend.position = 'bottom') +
  #     guides(fill = guide_legend(ncol = 7))
  # }
  aspect_ratio <- 1.5
  height <- 6
  if (ntaxa(ps_obj.comp.glom) > 45) {
    aspect_ratio <- 2
    if (ntaxa(ps_obj.comp.glom) > 100) {
      aspect_ratio <- 2
      height <- 8
      # if (ntaxa(ps_obj.comp.glom) > 300) {
      #   aspect_ratio <- 2.5
      #   height <- 8
      # }
    }
  }
    #geom_hline(yintercept=hline, linetype='dashed')
  top_bar$data[,taxa] <- factor(top_bar$data[,taxa], levels = taxa_order,
                                exclude = NULL)
  top_bar$data[,'Sample'] <- factor(top_bar$data[,'Sample'], 
         levels = sample_order)
  if (!is.na(facet)) {
    top_bar <- top_bar +
      facet_wrap(formula(paste0('~', facet)), scales = 'free')
  }
  save_plot(top_bar, filename = plot_name, 
            base_aspect_ratio = aspect_ratio, base_height = height)
  write.csv(top_bar$data, csv_name)
  if (retall) {
    ret_obj <- list()
    ret_obj$ps <- ps_obj.comp.glom
    ret_obj$plot <- top_bar
    if (cutoff < 100) {
      ret_obj$renamed <- renamed_taxa
    }
    return(ret_obj)
  } else {
    return(top_bar)
  }
}

plot_heatmap_top <- function(ps_obj, taxa = 'Phylum') {
  # hm_colors <- colorRampPalette(brewer.pal(n = 9, name ="Blues"))(100)
  n <- 50
  bc_df <- data.frame(t(otu_table(ps_obj)))
  legend_names <- scales::log_breaks(n = 7, base = 4)(range(bc_df) + 1)
  legend_names <- c(0, legend_names)
  legend_breaks <- log(legend_names+1, base = 4)
  # bc_df <- log(bc_df + 1, base = 2)
  bc_df <- log(bc_df + 1, base = 4)
  # breaks <- axisTicks(range(bc_df, na.rm = T), log = T, n = n)
  # hm_colors <- colorRampPalette(c('#000033', '#66CCFF'))(length(breaks)-1)
  hm_colors <- colorRampPalette(c('#000033', '#66CCFF'))(n)
  # legend_names_log2 <- c(0, 1, 2, 4, 8, 16, 32, 64)
  # legend_breaks_log2 <- log(legend_names+1, base = 2)
  # legend_breaks <- c(log2(1), log2(2),log2(3),log2(5),
  #                    log2(9),log2(17),log2(33),
  #                    log2(65))
  
  tax <- as.data.frame(as(tax_table(ps_obj), 'matrix')[,-1])
  tax <- Filter(function(x)!all(is.na(x)), tax)
  # Set defaults
  ratio = 4.0
  height = 8
  rowfont = 15
  colfont = 9
  
  # Modify for big # taxa
  mod_ratio <- floor(ntaxa(ps_obj)/35)
  height = height + mod_ratio*1.25
  rowfont = rowfont - mod_ratio
  colfont = colfont - mod_ratio
  ratio = ratio - mod_ratio*.1
  # if (mod_ratio) {
    # if ( taxa == 'Genus' ) {
    #   ratio = 0.9
    # }
    # if (mod_ratio == 1){
    #   ratio = 1.6
    # }
  # } 
  # else {
  #   if (ncol(tax) != 1) {
  #     ratio = 1.6 + ncol(tax)*.1
  #   }
  # }
  rownames(bc_df) <- do.call(paste, c(tax, sep=';'))
  otu_dists <- vegdist(otu_table(ps_obj), method = 'bray')
  sample_dists <- vegdist(t(otu_table(ps_obj)), method = 'bray')
  top_heat <- pheatmap(bc_df, clustering_distance_rows = sample_dists,
           clustering_distance_cols = otu_dists, 
           clustering_method = 'ward.D2',
           fontsize_row = rowfont,
           fontsize_col = colfont,
           color = hm_colors,
           legend_breaks = legend_breaks,
           legend_labels = legend_names)
  plot_name <- paste0('heatmap_top_', taxa, '.png')
  save_plot(top_heat, filename = plot_name, 
            base_aspect_ratio = ratio, base_height = height)
  return(top_heat)
}

bar_plot_facet <- function(ps, grouped, primary, secondary, cutoff=0.25,
                           plot_cutoff=5, rank='phylum', base_height=9,
                           base_asp=0.9, col_pal=pomological_palette,
                           reverse_y=FALSE) {
  plot_name <- paste0('facet_plot_', rank, '_', grouped, '.png')
  metadata <- ps %>%
    phyloseq::sample_data() %>%
    data.frame()
  ps.merge <- ps %>% 
    merge_samples(group = grouped) %>%
    suppressWarnings()
  ps.merge.glom <- ps.merge %>%
    tax_glom(taxrank = str_to_title(rank), NArm = FALSE)
  ps.merge.glom.comp <- ps.merge.glom %>% 
    transform_sample_counts(function(x) x/sum(x) * 100)
  plot_data <- ps.merge.glom.comp %>% psmelt()
  nsamples <- plot_data[,grouped] %>% n_distinct()
  taxon_level <- str_to_title(rank)
  primary_levels <- levels(metadata[,primary])
  secondary_levels <- levels(metadata[,secondary])
  frm <- paste('Abundance', '~', taxon_level)
  taxa_order <- plot_data %>%
    aggregate(formula(frm), ., FUN = sum) %>%
    dplyr::arrange(Abundance)
  taxa_rename <- taxa_order %>%
    mutate(Abundance = Abundance/nsamples) %>%
    filter(Abundance < cutoff) %>%
    pull(all_of(taxon_level))
  taxa_factor <- taxa_order %>%
    mutate(Abundance = Abundance/nsamples) %>%
    filter(Abundance >= cutoff) %>%
    pull(all_of(taxon_level))
  if (cutoff > 0) {
    new_label <- paste0('Other < ', cutoff, '%')
    taxa_factor <- c(new_label, taxa_factor)
    plot_data <- 
      plot_data %>%
      mutate("{taxon_level}" := ifelse(!!as.name(taxon_level) %in% taxa_rename,
                                       new_label, 
                                       !!as.name(taxon_level)))
  }
  if (any(is.na(plot_data[,taxon_level]))) {
    taxa_factor <- c(NA, taxa_factor)
  }
  if (reverse_y) {
    taxa_factor <- rev(taxa_factor)
  }
  plot_data <- plot_data %>%
    mutate("{taxon_level}" := factor(!!as.name(taxon_level),
                                     taxa_factor,
                                     exclude = NULL)) %>%
    mutate("{primary}" := factor(!!as.name(primary), labels = primary_levels)) %>%
    mutate("{secondary}" := factor(!!as.name(secondary), labels = secondary_levels))
  plot_data <- plot_data %>%
    group_by(!!as.name(primary), !!as.name(secondary), !!as.name(taxon_level)) %>%
    summarize(Abundance = sum(Abundance), .groups='drop')
  taxa_levels <- plot_data %>% pull(all_of(taxon_level)) %>% nlevels()
  frm.facet <- paste(taxon_level, '~', primary)
  plot_keeps <- plot_data %>% group_by(!!as.name(taxon_level)) %>%
    summarize(mean = mean(Abundance), .groups='drop') %>%
    filter(mean > plot_cutoff) %>% 
    pull(!!as.name(taxon_level)) %>%
    as.character()
  nkeeps <- n_distinct(plot_keeps)
  message(paste('Plotting', nkeeps, 'taxa.'))
  nkeeps <- ifelse(any(is.na(plot_keeps)),
                   nkeeps - 1,
                   nkeeps)
  if (nkeeps > length(col_pal)) {
    col_pal <- colorRampPalette(col_pal)(nkeeps)
  } else {
    col_pal <- col_pal[1:nkeeps]
  }
  breaks = rev(plot_keeps)
  if (reverse_y) {
    values = rev(col_pal)
  } else {
    values = col_pal
  }
  faceted_plot <- 
    plot_data %>% 
    filter(!!as.name(taxon_level) %in% plot_keeps) %>%
    ggplot(aes(x = !!as.name(secondary), y = Abundance,
               fill = !!as.name(taxon_level), color = !!as.name(taxon_level),
               label = round(Abundance, digits = 2))) +
    # geom_bar(stat = 'identity',
    #          position = 'stack') +
    geom_col() +
    facet_grid(formula(frm.facet),
               scales = 'free',
               space = 'free_x',
               drop = F,
               as.table = F) +
    cowplot::theme_cowplot() +
    scale_fill_manual(
      values = values,
      aesthetics = c('color', 'fill'),
      breaks = breaks,
      na.value = '#000000',
      ) +
    ylab('Relative Abundance (%)')+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95),
          strip.text.y = element_blank(),
          strip.background = element_rect(fill = 'white', colour="black",
                                          size = 1),
          strip.placement = 'inside',
          axis.title.x = element_blank(),
          ) +
    ggfittext::geom_bar_text(color = 'black', grow = F)
  cowplot::save_plot(plot_name, faceted_plot, base_height = base_height, 
            base_asp = base_asp)
  return(faceted_plot)
}

box_plot_facet <- function(ps, primary, cutoff=0.25,
                           plot_cutoff=5, rank='phylum', base_height=9,
                           col_pal=pomological_palette) {
  plot_name <- paste0('facet_boxplot_', rank, '_', primary, '.png')
  metadata <- ps %>%
    phyloseq::sample_data() %>%
    data.frame()
  
  boxplot_data <- ps %>%
    tax_glom(taxrank = str_to_title(rank), NArm = FALSE) %>%
    transform_sample_counts(function(x) x/sum(x) * 100) %>% 
    psmelt() %>%
    group_by(!!is.name(primary), !!is.name(rank)) %>%
    mutate(outlier = Abundance > median(Abundance) + 
                 IQR(Abundance)*1.5 | Abundance < median(Abundance) -
                 IQR(Abundance)*1.5
           )
  
  nsamples <- boxplot_data[,primary] %>% n_distinct()
  taxon_level <- str_to_title(rank)
  primary_levels <- levels(metadata[,primary])
  frm <- paste('Abundance', '~', taxon_level)
  taxa_order <- boxplot_data %>%
    aggregate(formula(frm), ., FUN = sum) %>%
    dplyr::arrange(Abundance)
  taxa_rename <- taxa_order %>%
    mutate(Abundance = Abundance/nsamples) %>%
    filter(Abundance < cutoff) %>%
    pull(all_of(taxon_level))
  taxa_factor <- taxa_order %>%
    mutate(Abundance = Abundance/nsamples) %>%
    filter(Abundance >= cutoff) %>%
    pull(all_of(taxon_level))
  if (cutoff > 0) {
    new_label <- paste0('Other < ', cutoff, '%')
    taxa_factor <- c(new_label, taxa_factor)
    boxplot_data <- 
      boxplot_data %>%
      mutate("{taxon_level}" := ifelse(!!as.name(taxon_level) %in% taxa_rename,
                                       new_label, 
                                       !!as.name(taxon_level)))
  }
  if (any(is.na(boxplot_data[,taxon_level]))) {
    taxa_factor <- c(NA, taxa_factor)
  }
  boxplot_data <- boxplot_data %>%
    mutate("{taxon_level}" := factor(!!as.name(taxon_level),
                                     taxa_factor,
                                     exclude = NULL)) %>%
    mutate("{primary}" := factor(!!as.name(primary), labels = primary_levels))
  # boxplot_data <- boxplot_data %>%
  #   group_by(!!as.name(primary), !!as.name(taxon_level)) %>%
  #   summarize(Abundance = sum(Abundance), .groups='keep')
  taxa_levels <- boxplot_data %>% pull(all_of(taxon_level)) %>% nlevels()
  frm.facet <- paste('~', primary)
  plot_keeps <- boxplot_data %>% group_by(!!as.name(taxon_level)) %>%
    summarize(mean = mean(Abundance), .groups='drop') %>%
    filter(mean > plot_cutoff) %>% 
    pull(!!as.name(taxon_level)) %>%
    as.character()
  nkeeps <- n_distinct(plot_keeps)
  message(paste('Plotting', nkeeps, 'taxa.'))
  nkeeps <- ifelse(any(is.na(plot_keeps)),
                   nkeeps - 1,
                   nkeeps)
  facet_boxplot <- boxplot_data %>% 
    filter(!!as.name(taxon_level) %in% plot_keeps) %>%
    ggplot(aes(x = Abundance, y = !!as.name(primary),
               fill = !!as.name(primary), color = !!as.name(primary),
               group = !!as.name(primary),
               label = round(Abundance, digits = 2))) +
    geom_boxplot()  +
    facet_wrap(formula(paste0('~', str_to_title(rank))),
               scales = 'free',
               as.table = TRUE) +
    cowplot::theme_cowplot() +
    scale_x_continuous(limits = c(0, NA),
                       # breaks = scales::pretty_breaks(n = 3),
                       ) +
    scale_y_discrete(limits = rev(primary_order_list)) +
    scale_fill_manual(values = xgfs.6,
                      aesthetics = c('color', 'fill'),
                      guide = guide_legend(reverse = FALSE),
                      breaks = primary_order_list) +
    xlab('Relative Abundance (%)') +
    stat_summary(geom = "crossbar",
                 show.legend = F,
                 width=0.65,
                 fatten=0,
                 color="white",
                 fun.data = function(x) {
                   return(c(y=median(x),
                            ymin=median(x),
                            ymax=median(x)))
                   }) +
    theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = 'white', colour="black",
                                          size = 1),
          strip.placement = 'inside',
          )
  cowplot::save_plot(plot_name, facet_boxplot, base_height = base_height)
  return(facet_boxplot)
}
