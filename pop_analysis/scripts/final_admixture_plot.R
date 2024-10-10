devtools::install_github("an-bui/calecopal")

suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(maps)
  library(scatterpie)
  library(ggrepel)
  library(calecopal)
})


#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

create_plots_tess <- function(coords_file, tess_file, fam, set_k, output_path) {

  coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))
  
  # Load the Tess data
  tess.in <- read_csv(tess_file)
  names(tess.in) <- c("Sample", "pop_assignment", "kval", "qvalue", "total_K", "best_k", "min_k")
  
  tess <- tess.in %>% filter(total_K == set_k)
  df_tess <- left_join(tess, coords, by = "Sample")
  
  wide_data <- df_tess %>%
    pivot_wider(names_from = kval, values_from = qvalue)
  wide_data <- wide_data %>% rowwise() %>%
    mutate(pop_name = names(select(., starts_with("K")))[which.max(c_across(starts_with("K")))])
  
  #get samp num from fam file
  df_fam <- read_delim(fam, col_names = FALSE)
  df_fam$samp_num <- row.names(df_fam)
  df_fam2 <- df_fam %>% select(X2, samp_num)
  wide_data <- left_join(wide_data, df_fam2, by = c("Sample" = "X2")) 
  df_tess <- left_join(df_tess, df_fam2, by = c("Sample" = "X2"))  

  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")
  
  # tess_wide <- tess.in %>% 
  #   pivot_wider(id_cols = c(Sample, pop_assignment, total_K), names_from = kval, values_from = qvalue) %>% 
  #   dplyr::rename(K_value = total_K)
  
  # tess_assign <- df_tess %>% 
  #   select(Sample, samp_num, kval, qvalue) %>% 
  #   group_by(Sample) %>% 
  #   #pivot_wider(id_cols = c("Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%    
  #   mutate(likely_assignment = kval[which.max(qvalue)],
  #         assingment_prob = max(qvalue)) %>% 
  #   arrange(likely_assignment, desc(assingment_prob)) %>%
  #   select(Sample, likely_assignment)%>%
  #   distinct(.keep_all = T) %>%
  #   ungroup() %>%
  #   mutate(sample_order = 1:n())

  tess_assign <- df_tess %>% 
    select(Sample, samp_num, kval, qvalue) %>% 
    group_by(Sample) %>% 
    #pivot_wider(id_cols = c("Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%    
    mutate(likely_assignment = kval[which.max(qvalue)],
          assingment_prob = max(qvalue)) %>% 
    pivot_wider(id_cols = c("likely_assignment", "Sample", "samp_num"), names_from = "kval", values_from = "qvalue") %>%
    arrange(likely_assignment, across(all_of(unique(df_tess$kval)), desc)) %>%
    ungroup() %>%
    mutate(sample_order = row_number())

  sample.order <- tess_assign %>% select(Sample, sample_order) %>% distinct(.keep_all = T)
  df_tess <- left_join(df_tess, sample.order, by = c("Sample"))
  wide_data <- left_join(wide_data, sample.order, by = c("Sample"))

  # plots -------------------------------------------------------------------

  # Plot 3: TESS Admixture Proportions
  kvals <- unique(df_tess$kval)

  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = wide_data, aes(x = Longitude, y = Latitude), cols = kvals, color = NA) +
    coord_fixed() +
    labs(title = "TESS Admixture Proportions", x = NULL, y = NULL, fill = "pop") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8)) +
    scale_fill_manual(values = selected_colors)

  p1_labs <- p1 +  geom_text_repel(data = wide_data, aes(x = Longitude, y = Latitude, label = sample_order), 
                    size = 2, color = "black", 
                    max.overlaps = 10,
                    force = 2,       
                    nudge_y = 0.3,  
                    segment.length = 0.8,  
                    segment.color = "black") +
                    labs(title = "Sample Labels")

   p2 <- df_tess %>% 
      ggplot(aes(x = factor(sample_order), y = qvalue, fill = factor(kval), text = sample_order)) +
      geom_col(aes(fill = factor(kval)), size = 0.1) +
      theme_minimal() +
      labs(x = NULL, y = "Ancestry") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = expansion(add = 1)) +
      theme(
        panel.spacing.x = unit(0.0, "lines"),
        axis.text.x = element_text(angle = 90, size = 4),
        panel.grid = element_blank(),
        strip.text.x = element_text(angle = 0, size = 6),
        legend.position = "none"
      ) + 
      #coord_flip() +
      xlab(NULL)  + theme(legend.position = "none") +
      scale_fill_manual(values = selected_colors)
      
    layout <- "
            AABB
            AABB
            AABB
            CCCC
            "

  combined_plot <- p1 + p1_labs + p2 + plot_annotation(tag_levels = 'A') +
    plot_layout(design = layout)
  
  ggplot2::ggsave(output_path, plot = combined_plot, width = 8.5, height = 6, units = "in", device = cairo_pdf)

  return(tess_assign)

}

create_plots_admixture <- function(coords_file, admixture_path,fam, set_k, output_path) {
  # Load the admixture data for the best K value
  coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))

  num_columns <- set_k
  admix_k <- set_k
  admix_best_k <- set_k 
  
  column_names <- paste0("K", 1:num_columns)  # Generate column names from K1 to K{num_columns}
  
  #for simplicity, just load admixture_data as the relavent Q. Easier to plot
  admixture_data <- read_delim(paste0(admixture_path, ".", admix_best_k, ".42.Q"), 
                               col_names = column_names, 
                               delim = " ")
  df_fam <- read_delim(fam, col_names = FALSE)
  df_fam$samp_num <- row.names(df_fam)

  admixture_data$Sample <- df_fam$X2
  admixture_data$samp_num <- df_fam$samp_num
  admixture_data <- inner_join(coords, admixture_data, by = "Sample")
  
  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")
  
  k_columns <- names(admixture_data)[startsWith(names(admixture_data), "K")]

  edit_admix <- admixture_data
  names(edit_admix) <- gsub("K", "pop", names(edit_admix)) #rename ancestral pops
  numpops <- grep("pop", names(edit_admix), value = T)
  
  cat_admx <- edit_admix %>% 
      pivot_longer(names_to = "popGroup", values_to = "prob", 
        cols = -c(Sample, samp_num, Longitude, Latitude)) %>%
      mutate(popGroup = gsub("pop", "K", popGroup))


  # simple ordering, not so good
  # cat_admx.wide <- cat_admx %>% 
  #     select(Sample, samp_num, popGroup, prob) %>% 
  #     pivot_wider(id_cols = c("Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%
  #     arrange_(.dots = c(numpops)) %>% print(n=100)
  #     mutate(sample_order = 1:n()) %>% 
  #     ungroup()

  #includes pop assignment if you want to use
  cat_admx.wide <- cat_admx %>% 
    select(Sample, samp_num, popGroup, prob) %>% 
    group_by(Sample) %>% 
    #pivot_wider(id_cols = c("Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%    
    mutate(likely_assignment = popGroup[which.max(prob)],
    assingment_prob = max(prob)) %>% 
    pivot_wider(id_cols = c("likely_assignment", "Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%
    arrange(likely_assignment, across(all_of(unique(cat_admx$popGroup)), desc)) %>%
    ungroup() %>%
    mutate(sample_order = row_number())

    # arrange(likely_assignment, desc(assingment_prob)) %>%
    # select(Sample, likely_assignment)%>%
    # distinct(.keep_all = T) %>%
    # ungroup() %>%
    # mutate(sample_order = 1:n())
    

  sample.order <- cat_admx.wide %>% select(Sample, sample_order) %>% distinct(.keep_all = T)
  cat_admx <- left_join(cat_admx, sample.order, by = c("Sample"))
  admixture_data <- left_join(admixture_data, sample.order, by = c("Sample"))

  # cat_admx.wide <- cat_admx %>%
  #   select(Sample, samp_num, popGroup, prob) %>%
  #   pivot_wider(id_cols = c("Sample", "samp_num"), names_from = "popGroup", values_from = "prob") %>%
  #   # Arrange by ancestry proportions: first by pop1, then pop2, and so on
  #   arrange(desc(pop1), desc(pop2), desc(pop3)) %>%
  #   mutate(sample_order = row_number()) # Assign sample order

  # # Merge back the sample order into the original tibble
  # cat_admx <- left_join(cat_admx, cat_admx.wide %>% select(Sample, sample_order), by = "Sample")
  # admixture_data <- left_join(admixture_data, sample.order, by = c("Sample"))


  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = admixture_data, aes(x = Longitude, y = Latitude), cols = k_columns, color = NA) +
    coord_fixed() +
    labs(title = "Admixture Proportions", x = NULL, y = NULL, fill = "pop") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8)) +
    scale_fill_manual(values = selected_colors)
  
  p1_labs <- p1 +  geom_text_repel(data = admixture_data, aes(x = Longitude, y = Latitude, label = sample_order), 
                    size = 2, color = "black", 
                    max.overlaps = 10,
                    force = 2,       
                    nudge_y = 0.3,  
                    segment.length = 0.8,  
                    segment.color = "black") +
                    labs(title = "Sample Labels")

  p2 <- cat_admx %>% 
      ggplot(aes(x = factor(sample_order), y = prob, fill = factor(popGroup), text = sample_order)) +
      geom_col(aes(fill = factor(popGroup)), size = 0.1) +
      theme_minimal() +
      labs(x = NULL, y = "Ancestry") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = expansion(add = 1)) +
      theme(
        panel.spacing.x = unit(0.0, "lines"),
        axis.text.x = element_text(angle = 90, size = 4),
        panel.grid = element_blank(),
        strip.text.x = element_text(angle = 0, size = 6),
        legend.position = "none"
      ) + 
      #coord_flip() +
      xlab(NULL)  + theme(legend.position = "none") +
      scale_fill_manual(values = selected_colors)

      
    layout <- "
            AABB
            AABB
            AABB
            CCCC
            "

  combined_plot <- p1 + p1_labs + p2 + plot_annotation(tag_levels = 'A') +
    plot_layout(design = layout)
  
  ggplot2::ggsave(output_path, plot = combined_plot, width = 8.5, height = 6, units = "in", device = cairo_pdf)
  
  return(cat_admx.wide)
}

set_k <- snakemake@params[["manual_k_assignment"]]
k_path <- snakemake@output[["outdir"]]
taxa <- snakemake@params[["project_id"]]

#select color palette
colors <- cal_palette("chaparral1")
selected_colors <- colors[c(1, 4, 5, 2, 3, 6)]

admx_assign <- create_plots_admixture(
  coords_file = snakemake@input[["coords"]],
  admixture_path = snakemake@params[["admixture_path"]],
  fam = snakemake@input[["fam"]],
  set_k = set_k,
  output_path = snakemake@output[["out1"]]
)

tess_assign <- create_plots_tess(
  coords_file = snakemake@input[["coords"]],
  tess_file = snakemake@input[["tess"]],
  fam = snakemake@input[["fam"]],
  set_k = set_k,
  output_path = snakemake@output[["out2"]]
)


dir.create(k_path)

if(snakemake@params[["method"]] == "admixture") {
  df_out <- admx_assign
} else {
  df_out <- tess_assign
}

df_out <- df_out %>%
  rowwise() %>%
  mutate(likely_assignment = if_else(get(likely_assignment) < snakemake@params[["threshold_qval"]], NA_character_, likely_assignment)) %>%
  ungroup() 

write_tsv(df_out, snakemake@output[["all_qvalues"]])

df_out %>% 
  mutate(likely_assignment = gsub("K","", likely_assignment)) %>%
  split(.$likely_assignment) %>% 
  walk2(., names(.), ~ write_tsv(.x %>% select(Sample), paste0(k_path, "/population_", .y, ".txt"), col_names = FALSE))

#testing
# combined_plot_admx <- create_plots_admixture(
#   coords_file = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/41-Cyanocitta.coords.txt",
#   admixture_path = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/admixture/Q_files/41-Cyanocitta",
#   fam = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/admixture/41-Cyanocitta.fam",
#   set_k = set_k,
#   output_path = "projects/41-Cyanocitta/results/GCA_026167965.1/pop_analysis/41-Cyanocitta_admixture.pdf"
# )

# combined_plot_tess <- create_plots_tess(
#   coords_file = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/41-Cyanocitta.coords.txt",
#   tess_file = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/41-Cyanocitta_TESS_qmatrix.csv",
#   fam = "projects/41-Cyanocitta/results/GCA_026167965.1/algatr/admixture/41-Cyanocitta.fam",
#   set_k = set_k,
#   output_path = "projects/41-Cyanocitta/results/GCA_026167965.1/pop_analysis/41-Cyanocitta_tess.pdf"
# )

#next is to add this to the checkpoint