################################################################################
# load file for common project stuff
################################################################################
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), "/_globalStuff.R"))

qLoad('data_tables/meta_ALL.csv')
meta=meta_ALL
################################################################################
# Description
################################################################################
script_title <- "rarefied_tables"

################################################################################
# Read in raw tables from QIIME2 output
# - create Matrix.names.csv for other scripts
# - ready tables for use in other scripts
# - drop low sample counts, rarefy, filter to minimum counts or relative abundance
################################################################################

################################################################################
# Switches
################################################################################
use_minimum_count_group<-making_filtered_matrix<-making_matrix_names<- "no"
making_raw_matrix<-making_rarefied_matrix<-look_at_counts<-testing<-"no"
testing    

# NOTE: set to 'yes' to run that section
# making_matrix_names     <- "yes"
# making_raw_matrix       <- "yes"
# look_at_counts          <- "yes"

# rarefying
making_rarefied_matrix  <- "yes"
minimum_counts           <- 1000
use_auto_minimum_count='yes'; keep_number=1 # will have a minimum of this many if it started with

# filtering
making_filtered_matrix  <- "yes"
filter_by_counts         <- 0
min_relative_abundance   <- 1 / 100000

# use_minimum_count_group <- "yes"
# testing                 <- "yes"

################################################################################
# Data set selection
################################################################################
# NOTE: these are expected to exist in _globalStuff.R
Set1 <- data_set_order <- data_sets <- "Set1"

################################################################################
# Common filters / settings
################################################################################


data_set_filter_group    <- "sample_type"  

# selected groups
selected_groups <- sample_type <- sample_type_order <- c("Skin", "Feces", "Soil")

# filter groups
filter_group1       <- "sample_type"
filter_group1_order <- sample_type_order

filter_group2       <- "primer_set"
filter_group2_order <- c("Standard")

filter_group3       <- filter_group1
filter_group3_order <- filter_group1_order

# rarefy settings
rarefy_by_group <- "sample_type"


################################################################################
# Paths
################################################################################

original_path            <- "data_tables/original/"
raw_matrix_path          <- "data_tables/raw_matrix"
rarefied_matrix_path     <- "data_tables/rarefied_matrix"
filtered_matrix_path     <- "data_tables/filtered_matrix"

################################################################################
# Sanity prints (optional)
################################################################################
print(data_set_order)

qPrint(data_sets)
for (data_set in data_sets) {
  print(get(sym(data_set)))
}

################################################################################
# make matrix names
# - a file to help work on different data sets
################################################################################
if (making_matrix_names == "yes") { qPrint(making_matrix_names)
  
  create_directory(raw_matrix_path)
  create_directory(rarefied_matrix_path)
  create_directory(filtered_matrix_path)
  
  names(taxa_plural) <- taxa_levels
  
  matrix_names <- expand.grid(taxa_levs = taxa_levels, data_sets = data_sets) %>%
    mutate(taxa_plural = case_when(
      taxa_levs %in% names(taxa_plural) ~ taxa_plural[taxa_levs],
      TRUE ~ taxa_levs
    )) %>%
    mutate(raw_path      = paste0(raw_matrix_path, "/raw_matrix_", taxa_levs, ".csv")) %>%
    mutate(rarefied_path = paste0(rarefied_matrix_path, "/rarefied_matrix_", taxa_levs, "_", data_sets, ".csv")) %>%
    mutate(filtered_path = paste0(filtered_matrix_path, "/filtered_matrix_", taxa_levs, "_", data_sets, ".csv")) %>%
    mutate(file_path = filtered_path) %>%
    arrange(taxa_levs)
  
  write.csv(matrix_names, "data_tables/matrix_names.csv", row.names = FALSE)
  print("created matrix_names.csv")
  message()
  stop("\tdeselect making matrix names to continue")
}

################################################################################
# make raw matrix
# - create a matrix for each data set in a standard form
################################################################################
if (making_raw_matrix == "yes") { qPrint(making_raw_matrix)
  
  sample_reads_df2 <- NULL
  
  ################################################################################
  # Singletons
  ################################################################################
  for (taxa_levs in taxa_levels) { qPrint(taxa_levs)
    
    if (taxa_levs == "ASV") {
      file_name <- "Table_Features"
    } else if (taxa_levs == "Contaminated") {
      file_name <- "collapsed_table"
    } else {
      file_name <- paste0("Table_", taxa_levs)
    }
    
    raw_matrix_combined <- data.frame()
    
    for (type in sample_type_order) {
      
      raw_matrix <- read.delim(
        paste0(original_path, type, "/", file_name, ".tsv"),
        check.names = FALSE,
        row.names = 1,
        skip = 1
      ) %>%
        t() %>%
        as.data.frame()
      
      feature_names <- names(raw_matrix)
      
      # fix sample_name
      raw_matrix2 <- raw_matrix %>%
        rownames_to_column(var = "sample_ID") %>%
        left_join(meta %>% filter(data_class == "singleton")) %>%
        select(any_of(meta_names), everything()) %>%
        filter(primer_set == "Standard") %>%
        filter(sample_type %in% selected_groups) %>%
        filter(!is.na(sample_name)) %>%
        imPlode_sample_name() %>%
        select(which(colSums(.) > 0))
      
      feature_names <- names(raw_matrix)
      
      raw_matrix_long <- raw_matrix2 %>%
        xPlode_sample_name() %>%
        pivot_longer(cols = any_of(feature_names), names_to = "features", values_to = "counts")
      
      raw_matrix_combined <- bind_rows(raw_matrix_combined, raw_matrix_long) %>%
        mutate(across(everything(), ~ replace_na(., 0)))
    }
    
    ################################################################################
    # Doubletons
    ################################################################################
    if (taxa_levs == "ASV") {
      file_name <- "Table_Features_PostRemoval"
    } else if (taxa_levs == "Contaminated") {
      file_name <- "collapsed_table"
    } else {
      file_name <- paste0("Table_", taxa_levs)
    }
    
    doubleton <- read.delim(
      paste0(original_path, "Doubleton/", file_name, ".tsv"),
      check.names = FALSE,
      row.names = 1,
      skip = 1
    ) %>%
      t() %>%
      as.data.frame()
    
    feature_names <- names(doubleton)
    
    # fix sample_name
    doubleton2 <- doubleton %>%
      rownames_to_column(var = "sample_ID") %>%
      left_join(meta %>% filter(data_class == "doubleton")) %>%
      select(any_of(meta_names), everything()) %>%
      filter(primer_set == "Standard") %>%
      filter(sample_type %in% selected_groups) %>%
      filter(experiment %in% c("G4", "G6", "G7")) %>%
      filter(!is.na(sample_name)) %>%
      imPlode_sample_name() %>%
      select(which(colSums(.) > 0))
    
    feature_names <- names(doubleton)
    
    doubleton_long <- doubleton2 %>%
      xPlode_sample_name() %>%
      pivot_longer(cols = any_of(feature_names), names_to = "features", values_to = "counts")
    
    ################################################################################
    # Combine singleton + doubleton, pivot wide, and write
    ################################################################################
    raw_matrix1 <- bind_rows(raw_matrix_combined, doubleton_long) %>%
      mutate(across(everything(), ~ replace_na(., 0))) %>%
      pivot_wider(names_from = "features", values_from = "counts", values_fill = list(counts = 0)) %>%
      imPlode_sample_name()
    
    sample_reads <- rowSums(raw_matrix1) %>%
      as.data.frame() %>%
      setNames(paste("Reads")) %>%
      mutate(taxa_level = paste0(taxa_levs, "_Reads")) %>%
      xPlode_sample_name()
    
    sample_reads_df2 <- rbind(sample_reads_df2, sample_reads)
    
    if (any(is.na(raw_matrix1))) stop("unfiltered_df contains NA values")
    
    # transpose for quicker reading
    raw_matrix1 %>%
      t() %>%
      as.data.frame() %>%
      write.csv(paste0(raw_matrix_path, "/raw_matrix_", taxa_levs, ".csv"))
  }
  
  print("created raw matrixes")
  
  sample_reads_df <- sample_reads_df2 %>%
    pivot_wider(names_from = taxa_level, values_from = Reads) %>%
    imPlode_sample_name()
  
  write.csv(sample_reads_df, paste0("data_tables/sample_reads.csv"))
  print("created raw matrices.csv")
  stop("\tdeselect making raw matrix to continue")
} # end making_raw_matrix

################################################################################
# look at counts
################################################################################
# visualizing cut off effects
################################################################################
if (look_at_counts == "yes") { qPrint(look_at_counts)
  
  taxa_levs <- "Phylum"
  
  matrix_df <- read.csv(
    paste0(raw_matrix_path, "/raw_matrix_", taxa_levs, ".csv"),
    row.names = 1,
    check.names = FALSE
  ) %>%
    t() %>%
    as.data.frame()
  
  feature_names <- names(matrix_df)
  matrix_df$sample_counts <- rowSums(matrix_df)
  
  df <- matrix_df %>%
    select(sample_counts) %>%
    xPlode_sample_name() %>%
    ungroup() %>%
    mutate(short_name = factor(short_name, levels = unique(short_name)))
  
  ggplot(df, aes(x = short_name, y = sample_counts, color = temperature, shape = normalization)) +
    geom_point(aes(color = normalization), size = 5, alpha = 0.2) +
    geom_point() +
    scale_y_log10() +
    geom_hline(yintercept = 1000, linetype = "dashed", color = "red") +
    facet_grid(sample_type ~ data_class + temperature)
  
  ggsave(paste(output_plot, "looking at raw counts", project_folder, ".png"), width = 12, height = 12)
  
  stop()
  
  ggplot(df, aes(x = short_name, y = sample_counts)) +
    geom_point(aes(color = Temperature)) +
    geom_hline(yintercept = 10000, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90))
  
  p <- ggplot(df, aes(x = short_name, y = sample_counts)) +
    geom_boxplot(aes(color = temperature), outlier.size = 1) +
    geom_hline(yintercept = 10000, linetype = "dashed", color = "red") +
    facet_grid(~ plate, scales = "free_x") +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") +
    ggtitle("plates 1,3,5 (24 cycles); plates 2,4,6 (targeted fluorescence)\nplate 7 (QB only: fold change; QB and AMPED: slope and fixed cycles)")
  
  gPlot(p)
  ggsave("counts by plate and temperature.png", width = 10, height = 8)
  
  # losing at most one sample from each primer
  minimum_counts_df <- df %>%
    group_by(Exp, Primer) %>%
    arrange(sample_counts) %>%
    slice(2) %>%
    group_by(Exp) %>%
    arrange(sample_counts) %>%
    slice(1) %>%
    ungroup()
  
  minimum_counts_df
  
  Exp_df <- df %>%
    filter(Exp == "E2")
  
  print(Exp_df)
  
  ggplot(df, aes(x = Exp, y = sample_counts)) +
    geom_violin() +
    geom_point() +
    scale_y_log10() +
    facet_wrap(~ Type, scales = "free", nrow = 1)
  
  minimum_counts_df
  
  stop("choose minimum counts and deselect looking to continue")
} # end look_at_counts

################################################################################
# subsample matrices (rarefy)
################################################################################
if (making_rarefied_matrix == "yes") { qPrint(making_rarefied_matrix)
  
  ################################################################################
  # Functions
  ################################################################################
  
  # auto_min_counts: max of slice1
  auto_min_counts <- function(sample_count_dfx, keep_number) {
    min_count <- sample_count_dfx %>%
      group_by(short_name) %>%
      arrange(desc(sample_count)) %>%
      slice(1:keep_number) %>%
      ungroup() %>%
      pull(sample_count) %>%
      min()
    return(min_count)
  }
  
  # rarefy function
  rarefy_function <- function(dfx, min_sample_count) {
    df <- dfx %>%
      filter(sample_count >= min_sample_count) %>%
      select(-sample_count) %>%
      imPlode_sample_name() %>%
      as.matrix()
    
    set.seed(20250314)
    
    rarefied_df <- rrarefy(df, sample = min_sample_count) %>%
      as.data.frame()
    
    return(rarefied_df)
  }
  
  ################################################################################
  # prepping
  ################################################################################
  
  matrix_names <- read.csv("data_tables/matrix_names.csv", check.names = FALSE)
  
  r <- 1
  n_rows <- 1
  if (testing == "no") n_rows <- nrow(matrix_names)
  
  # load matrix from names csv
  for (r in 1:n_rows) {
    
    taxa_levs <- matrix_names[r, "taxa_levs"]
    data_set  <- matrix_names[r, "data_sets"]
    p_title   <- paste0(script_title, "_", taxa_levs, "_", data_set)
    
    data_set_group <- get(sym(data_set))
    qPrint(data_set)
    qPrint(data_set_group)
    
    matrix_df <- read.csv(matrix_names[r, "raw_path"], check.names = FALSE, row.names = 1) %>%
      t() %>%
      as.data.frame()
    
    feature_names <- names(matrix_df)
    
    ##############################################################################
    message("388 # Filtering for data sets and selected group")
    ##############################################################################
    data_set_df <- matrix_df %>%
      mutate(sample_count = rowSums(.)) %>%
      xPlode_sample_name() %>%
      filter(.data[[filter_group1]] %in% filter_group1_order) %>%
      filter(.data[[filter_group2]] %in% filter_group2_order) %>%
      filter(.data[[filter_group3]] %in% filter_group3_order) %>%
      # filter(.data[[data_set_filter_group]] %in% data_set_group) %>%
      filter(.data[[rarefy_by_group]] %in% selected_groups) %>%
      ungroup()
    
    ##############################################################################
    # looping by rarefy_by_group while rarefying
    ##############################################################################
    rarefied_df <- NULL
    
    if (exists("rarefy_by_group")) {
      
      # loop over unique values in the rarefy_by_group column
      for (t in unique(data_set_df[[rarefy_by_group]])) {
        qPrint(t)
        
        subset_df <- data_set_df %>%
          filter(.data[[rarefy_by_group]] == t)
        
        # Determine minimum sample count
        min_sample_count <- minimum_counts
        if (use_auto_minimum_count == "yes") {
          min_sample_count <- auto_min_counts(subset_df, keep_number)
        }
        
        # NOTE: per-group overrides (keep as-is)
        if (t == "Skin")  min_sample_count <- 35000
        if (t == "Feces") min_sample_count <- 25000
        if (t == "Soil")  min_sample_count <- 10000
        
        rarefied_subset <- rarefy_function(subset_df, min_sample_count)
        
        rarefied_df <- bind_rows(rarefied_df, rarefied_subset) %>%
          mutate(across(everything(), ~ replace_na(., 0)))
      }
      
      print(
        rarefied_df %>%
          mutate(sample_count = rowSums(.)) %>%
          xPlode_sample_name() %>%
          select(any_of(rarefy_by_group), sample_count) %>%
          distinct()
      )
      
    } else {
      
      message("\t No subsetting")
      
      min_sample_count <- minimum_counts
      if (use_auto_minimum_count == "yes") {
        min_sample_count <- auto_min_counts(data_set_df, keep_number)
      }
      
      rarefied_df <- rarefy_function(data_set_df, min_sample_count)
      
      print(table(rowSums(rarefied_df)))
    }
    
    print(paste("writing", Sys.time(), "writing")); cat("\n")
    
    rarefied_df %>%
      t() %>%
      as.data.frame() %>%
      write.csv(paste0(rarefied_matrix_path, "/rarefied_matrix_", taxa_levs, "_", data_set, ".csv"))
  }
} # end making_rarefied_matrix

################################################################################
# filter rarefied matrices
################################################################################
if (making_filtered_matrix == "yes") { qPrint(making_filtered_matrix)
  
  matrix_names <- read.csv("data_tables/matrix_names.csv", check.names = FALSE)
  
  r <- 1
  n_rows <- 1
  if (testing == "no") n_rows <- nrow(matrix_names)
  
  # load matrix from names csv
  for (r in 1:n_rows) {
    
    taxa_levs <- matrix_names[r, "taxa_levs"]
    data_set  <- matrix_names[r, "data_sets"]
    p_title   <- paste0(script_title, "_", taxa_levs, "_", data_set)
    
    qPrint(p_title)
    
    matrix_df <- read.csv(matrix_names[r, "rarefied_path"], check.names = FALSE, row.names = 1) %>%
      t() %>%
      as.data.frame()
    
    feature_names <- names(matrix_df)
    
    df <- matrix_df %>%
      select(which(colSums(.) > filter_by_counts)) %>%
      t() %>%
      as.data.frame()
    
    write.csv(df, paste0(filtered_matrix_path, "/filtered_matrix_", taxa_levs, "_", data_set, ".csv"))
  }
} # end making_filtered_matrix
