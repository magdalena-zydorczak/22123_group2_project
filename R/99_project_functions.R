quantile_normalization <- function(raw_counts){
  raw_counts_ranks <- raw_counts |>
    mutate(across(!gene,
                  ~rank(., ties.method = "min")))
  
  raw_counts_sorted_cols <- raw_counts |>
    mutate(across(!gene,
                  sort))
  
  raw_counts_means <- raw_counts_sorted_cols |>
    rowwise() |>
    summarise(m = mean(c_across(!gene))) |>
    ungroup() |>
    rowid_to_column("ID")
  
  lookup_vec <- setNames(raw_counts_means$m,
                         raw_counts_means$ID)
  
  quantile_normalized_counts <- raw_counts_ranks |>
    mutate(across(!gene,
                  ~ lookup_vec[.] |>
                    as.numeric()))
  
  return(quantile_normalized_counts)
}

#Method from: https://pubmed.ncbi.nlm.nih.gov/28650338/ (including quantile normalization and log10 transformation)
#Returns sample metadata with augmented GEP scores
GEP_score_calculation <- function(quantile_normalized_counts, metadata, GEP_genes){
  signature_genes <- tibble(
    GEP_gene = GEP_genes
  )
  
  augmented_metadata <- quantile_normalized_counts |>
    mutate(across(!gene,
                  ~log10(.+1))) |>
    right_join(signature_genes,
               join_by(gene == GEP_gene)) |>
    summarise(across(!gene, mean)) |>
    pivot_longer(everything(),
                 names_to = "pat_id",
                 values_to = "GEP_score") |>
    right_join(metadata,
               join_by(pat_id))
  
  return(augmented_metadata)
}