setwd("C:/Users/USER/Desktop/YOGHURT/Updated-V-July-2024/stamp/Previous")
#https://github.com/cafferychen777/ggpicrust2?tab=readme-ov-file#installation

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(phyloseq)
library("openxlsx")
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

#abundance_file <- "C:/Users/USER/Desktop/YOGHURT/Updated-V-July-2024/stamp/Previous/ko_metagenome-Y.tsv"

#metadata <- read_delim("Metadata.txt", delim = "/t", escape_double = FALSE, trim_ws = TRUE)

ko_abundance2<- read.xlsx("kegg-pic.xlsx", sheet = "kegg-y", colNames =TRUE, rowNames = TRUE)
metadata2<- read.xlsx("kegg-pic.xlsx", sheet = "m", colNames =TRUE, rowNames = TRUE)

#Example data
#data(ko_abundance)
#data(metadata)

results_file_input <- ggpicrust2(data = ko_abundance2,
                                 metadata = metadata2,
                                 group = "pH",
                                 pathway = "KO",
                                 daa_method = "ALDEx2",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")



results_file_input[[1]]$plot

results_file_input[[1]]$results

kegg_abundance <- ko2kegg_abundance(data = ko_abundance2)

daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata2, group = "Group", daa_method = "LinDA")

pathway_annotation

#daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
daa_annotated_results_df2 <- PA(pathway = "EC", daa_results_df = daa_results_df, ko_to_kegg = TRUE)

p <- pathway_errorbar(abundance = kegg_abundance,
                      daa_results_df = daa_annotated_results_df2,
                      Group = metadata2$Group,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.6,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")

##########export##################################
p1 = data.frame(kegg_abundance)
write.xlsx(x = p1, file = "S1.xlsx", col.names = TRUE, row.names = TRUE)  

p2 = data.frame(daa_annotated_results_df2)
write.xlsx(x = p2, file = "S2.xlsx", col.names = TRUE, row.names = TRUE) 

library("writexl")
install.packages("writexl")

p1 = data.frame(kegg_abundance)
write.xlsx(x = p1, file = "Sample_all_V2.xlsx", col.names = TRUE, row.names = TRUE)  

p2 = data.frame(daa_annotated_results_df2)
write.xlsx(x = df2, file = "SV4.xlsx")

write_xlsx(daa_annotated_results_df2, path = "SV5.xlsx")


df = data.frame(p2$pathway_name)
df2 = t(df)
write_xlsx(df2, path = "SV5.xlsx")
write.xlsx(x = df, file = "df.xlsx", col.names = TRUE, row.names = TRUE)  











############################################################################
PA <- function (file = NULL, pathway = NULL, daa_results_df = NULL, 
                ko_to_kegg = FALSE) 
{
  message("Starting pathway annotation...")
  if (is.null(file) && is.null(daa_results_df)) {
    stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
  }
  if (!is.null(file)) {
    message("Reading the input file...")
    file_format <- substr(file, nchar(file) - 3, nchar(file))
    switch(file_format, .txt = {
      message("Loading .txt file...")
      abundance <- readr::read_delim(file, delim = "\t", 
                                     escape_double = FALSE, trim_ws = TRUE)
      message(".txt file successfully loaded.")
    }, .tsv = {
      message("Loading .tsv file...")
      abundance <- readr::read_delim(file, delim = "\t", 
                                     escape_double = FALSE, trim_ws = TRUE)
      message(".tsv file successfully loaded.")
    }, .csv = {
      message("Loading .csv file...")
      abundance <- readr::read_delim(file, delim = "\t", 
                                     escape_double = FALSE, trim_ws = TRUE)
      message(".csv file successfully loaded.")
    }, stop("Invalid file format. Please input file in .tsv, .txt or .csv format. The best input file format is the output file from PICRUSt2, specifically 'pred_metagenome_unstrat.tsv'."))
    abundance <- abundance %>% tibble::add_column(description = rep(NA, 
                                                                    length = nrow(abundance)), .after = 1)
    switch(pathway, KO = {
      message("Loading KO reference data...")
      load(system.file("extdata", "KO_reference.RData", 
                       package = "ggpicrust2"))
      message("Annotating abundance data with KO reference...")
      for (i in seq_len(nrow(abundance))) {
        abundance[i, 2] <- KO_reference[KO_reference[, 
                                                     1] %in% abundance[i, 1], 5][1]
      }
      message("Abundance data annotation with KO reference completed.")
    }, EC = {
      message("Loading EC reference data...")
      load(system.file("extdata", "EC_reference.RData", 
                       package = "ggpicrust2"))
      message("Annotating abundance data with EC reference...")
      for (i in seq_len(nrow(abundance))) {
        abundance[i, 2] <- EC_reference[EC_reference[, 
                                                     1] %in% abundance[i, 1], 2]
      }
      message("Abundance data annotation with EC reference completed.")
      message("Note: EC description may appear to be duplicated due to shared EC numbers across different reactions.")
    }, MetaCyc = {
      message("Loading MetaCyc reference data...")
      load(system.file("extdata", "MetaCyc_reference.RData", 
                       package = "ggpicrust2"))
      message("Annotating abundance data with MetaCyc reference...")
      for (i in seq_len(nrow(abundance))) {
        abundance[i, 2] <- MetaCyc_reference[MetaCyc_reference[, 
                                                               1] %in% abundance[i, 1], 2]
      }
      message("Abundance data annotation with MetaCyc reference completed.")
    }, stop("Invalid pathway option. Please provide one of the following options: 'KO', 'EC', 'MetaCyc'."))
    return(abundance)
  }
  if (!is.null(daa_results_df)) {
    message("DAA results data frame is not null. Proceeding...")
    if (ko_to_kegg == FALSE) {
      message("KO to KEGG is set to FALSE. Proceeding with standard workflow...")
      daa_results_df$description <- NA
      switch(pathway, KO = {
        message("Loading KO reference data...")
        load(system.file("extdata", "KO_reference.RData", 
                         package = "ggpicrust2"))
        for (i in seq_len(nrow(daa_results_df))) {
          daa_results_df[i, ]$description <- KO_reference[KO_reference[, 
                                                                       1] %in% daa_results_df[i, ]$feature, 5][1]
        }
      }, EC = {
        message("Loading EC reference data...")
        load(system.file("extdata", "EC_reference.RData", 
                         package = "ggpicrust2"))
        for (i in seq_len(nrow(daa_results_df))) {
          daa_results_df[i, ]$description <- EC_reference[EC_reference[, 
                                                                       1] %in% daa_results_df[i, ]$feature, 2]
        }
        message("EC description may appear to be duplicated")
      }, MetaCyc = {
        message("Loading MetaCyc reference data...")
        load(system.file("extdata", "MetaCyc_reference.RData", 
                         package = "ggpicrust2"))
        for (i in seq_len(nrow(daa_results_df))) {
          daa_results_df[i, ]$description <- MetaCyc_reference[MetaCyc_reference[, 
                                                                                 1] %in% daa_results_df[i, ]$feature, 2]
        }
      }, stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway"))
      message("Returning DAA results data frame...")
      return(daa_results_df)
    }
    else {
      message("KO to KEGG is set to TRUE. Proceeding with KEGG pathway annotations...")
      daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < 
                                                  0.999999, ]
      if (nrow(daa_results_filtered_df) == 0) {
        stop("No statistically significant biomarkers found. 'Statistically significant biomarkers' refer to those biomarkers that demonstrate a significant difference in expression between different groups, as determined by a statistical test (p_adjust < 0.05 in this case).\n", 
             "You might consider re-evaluating your experiment design or trying alternative statistical analysis methods. Consult with a biostatistician or a data scientist if you are unsure about the next steps.")
      }
      daa_results_filtered_df$pathway_name <- NA
      daa_results_filtered_df$pathway_description <- NA
      daa_results_filtered_df$pathway_class <- NA
      daa_results_filtered_df$pathway_map <- NA
      keggGet_results <- list()
      message("We are connecting to the KEGG database to get the latest results, please wait patiently.")
      if (nrow(daa_results_filtered_df) > 100) {
        cat("\n")
        message("The number of statistically significant pathways exceeds the database's query limit. Please consider breaking down the analysis into smaller queries or selecting a subset of pathways for further investigation.")
        cat("\n")
      }
      if (nrow(daa_results_filtered_df) <= 10) {
        cat("\n")
        message("Processing pathways individually...")
        cat("\n")
        pb <- txtProgressBar(min = 0, max = nrow(daa_results_filtered_df), 
                             style = 3)
        for (i in seq_len(nrow(daa_results_filtered_df))) {
          cat("\n")
          message("Beginning annotation for pathway ", 
                  i, " of ", nrow(daa_results_filtered_df), 
                  "...")
          cat("\n")
          a <- 0
          start_time <- Sys.time()
          repeat {
            tryCatch({
              keggGet_results[[i]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[i])
              a <- 1
            }, error = function(e) {
              cat("\n")
              message("An error occurred. Retrying...")
              cat("\n")
            })
            if (a == 1) {
              break
            }
          }
          end_time <- Sys.time()
          time_taken <- end_time - start_time
          cat("\n")
          message("Annotated pathway ", i, " of ", nrow(daa_results_filtered_df), 
                  ". Time taken: ", round(time_taken, 2), " seconds.")
          cat("\n")
          daa_results_filtered_df[i, ]$pathway_name <- keggGet_results[[i]][[1]]$NAME
          daa_results_filtered_df[i, ]$pathway_description <- keggGet_results[[i]][[1]]$DESCRIPTION[1]
          daa_results_filtered_df[i, ]$pathway_class <- keggGet_results[[i]][[1]]$CLASS
          daa_results_filtered_df[i, ]$pathway_map <- keggGet_results[[i]][[1]]$PATHWAY_MAP
          setTxtProgressBar(pb, i)
        }
        close(pb)
        cat("\n")
        message("Pathway annotation completed.")
        cat("\n")
      }
      if (nrow(daa_results_filtered_df) > 10 && nrow(daa_results_filtered_df) < 
          300) {
        cat("\n")
        message("Processing pathways in chunks...")
        cat("\n")
        pb <- txtProgressBar(min = 0, max = nrow(daa_results_filtered_df), 
                             style = 3)
        start_time <- Sys.time()
        n <- length(c(seq(10, nrow(daa_results_filtered_df), 
                          10), nrow(daa_results_filtered_df)))
        j <- 1
        seq <- c(seq(10, nrow(daa_results_filtered_df), 
                     10), nrow(daa_results_filtered_df))
        for (i in seq) {
          if (i%%10 == 0) {
            keggGet_results[[j]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[seq(i - 
                                                                                            9, i, 1)])
          }
          else {
            keggGet_results[[j]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[seq(nrow(daa_results_filtered_df)%/%10 * 
                                                                                            10 + 1, i, 1)])
          }
          j <- j + 1
          setTxtProgressBar(pb, i)
        }
        end_time <- Sys.time()
        time_taken <- end_time - start_time
        cat("\n")
        message("Finished processing chunks. Time taken: ", 
                round(time_taken, 2), " seconds.")
        cat("\n")
        message("Finalizing pathway annotations...")
        cat("\n")
        start_time <- Sys.time()
        for (k in 1:n) {
          w <- length(keggGet_results[[k]])
          for (j in 1:w) {
            daa_results_filtered_df[daa_results_filtered_df$feature == 
                                      keggGet_results[[k]][[j]]$ENTRY, ]$pathway_name <- keggGet_results[[k]][[j]]$NAME[1]
            daa_results_filtered_df[daa_results_filtered_df$feature == 
                                      keggGet_results[[k]][[j]]$ENTRY, ]$pathway_description <- keggGet_results[[k]][[j]]$DESCRIPTION[1]
            daa_results_filtered_df[daa_results_filtered_df$feature == 
                                      keggGet_results[[k]][[j]]$ENTRY, ]$pathway_class <- keggGet_results[[k]][[j]]$CLASS[1]
            daa_results_filtered_df[daa_results_filtered_df$feature == 
                                      keggGet_results[[k]][[j]]$ENTRY, ]$pathway_map <- keggGet_results[[k]][[j]]$PATHWAY_MAP[1]
          }
          setTxtProgressBar(pb, k)
        }
        end_time <- Sys.time()
        time_taken <- end_time - start_time
        cat("\n")
        message("Finished finalizing pathway annotations. Time taken: ", 
                round(time_taken, 2), " seconds.")
        cat("\n")
        close(pb)
      }
      daa_results_filtered_annotation_df <- daa_results_filtered_df
      message("Returning DAA results filtered annotation data frame...")
      return(daa_results_filtered_annotation_df)
    }
  }
}
###################################################################
