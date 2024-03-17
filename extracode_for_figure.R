Chk <- classes %>% select(Metabolites) %>% left_join(sigfile_to_export, by = c("Metabolites"))
df <- Chk

# Find Metabolites with more than one associated Phthalate
duplicates <- df %>%
  group_by(Metabolites) %>%
  filter(n() > 1)

# Assign a suffix to the Phthalates for these Metabolites
duplicates <- duplicates %>%
  group_by(Metabolites) %>%
  mutate(suffix = row_number()) %>%
  ungroup()

# Spread the values for the identified Metabolites
duplicates_wide <- duplicates %>%
  pivot_wider(
    id_cols = Metabolites,
    names_from = suffix,
    values_from = c(Phthalates, Class_new2, beta_ci, adj_p),
    names_sep = "_"
  )

# Find Metabolites with only one associated Phthalate and assign "_1" suffix
singletons <- df %>%
  group_by(Metabolites) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(suffix = 1) %>%
  pivot_wider(
    id_cols = Metabolites,
    names_from = suffix,
    values_from = c(Phthalates, Class_new2, beta_ci, adj_p),
    names_sep = "_"
  )

# Merge the results
result <- bind_rows(duplicates_wide, singletons)

# Print the resulting data frame
result

# Create a new column based on the conditions for adj_p1 through adj_p3
result <- result %>%
  rowwise() %>%
  mutate(Significance = case_when(
    all(is.na(c(adj_p_1, adj_p_2, adj_p_3))) ~ NA_character_,
    any(c(adj_p_1, adj_p_2, adj_p_3) < 0.05, na.rm = TRUE) ~ "Sig",
    TRUE ~ "Marginal sig"
  )) %>%
  ungroup()

# Print the resulting data frame
result


#select just metabolites and significance column

result <- result %>% select(Metabolites, Significance) 
write.csv(result, "prenatal_metfile_sig_asso.csv")



