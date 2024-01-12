tbi_pval <- qs::qread('./tables/tbi_xy/tbiDF_2091_allgcms.qs')
tbi_pval<- tbi_pval %>% filter(gcm == 'CanESM2' & pValue2 == 'significant')
qs::qsave(tbi_pval, './tables/sigpValTBICanESM2_2091.qs')

lu_table <- qs::qread('./tables/sigpValTBICanESM2_2091.qs')
out_dir <- './outputs/Aggregations'
data <- list.files(out_dir, pattern =  'aggs_', full.names = TRUE)
data <- qs::qread('./outputs/Aggregations/aggs_CanESM2_2091.qs')


#this function filters the list to select those lists with more of 5 aggs
# Filter the list based on the maximum value of 'aggs'
filtered_list <- lapply(data, function(df) {
  max_aggs <- max(df$aggs, na.rm = TRUE)
  if (max_aggs > 5) {
    df
  } else {
    NULL
  }
})

# Remove NULL elements from the list
filtered_list <- filtered_list[sapply(filtered_list, function(x) !is.null(x))]

## Create a function to filter data frames based on coordinates 
filter_df <- function(data, lu_table) {
  n <- length(names(data))
  
  result <- lapply(seq_len(n), function(i) {
    # Display progress message
    message(paste("Processing:", i, "of", n))
    
    # Extract numeric values from the df name
    name_parts <- str_extract_all(names(data)[i], "-?\\d+\\.\\d+")[[1]]
    y_value <- name_parts[1]
    x_value <- name_parts[2]
    
    # Convert x and y columns in lu_table to character for comparison
    lu_table_chr <- lu_table %>%
      mutate(x = as.character(Longitude), y = as.character(Latitude))
    
    # Check if the coordinates match the conditions
    match_rows <- filter(lu_table_chr, x == x_value, y == y_value)
    if (nrow(match_rows) > 0) {
      # If there is a match, return the df
      return(data = data[[i]])
    } else {
      # If not match, return NULL
      return(NULL)
    }
  })
  #set names for the resulting list using original names 
  result <- setNames(result, names(data))
  
  return(result)
}


# Call the function
filtered_list <- filter_df(data, lu_table)


# Remove NULL elements from the list
filtered_list <- filtered_list[!sapply(filtered_list, is.null)]
qs::qsave(filtered_list, glue('./tables/sigAggs/aggsSigpValTBI_{gcm}_{yr}.qs'))



