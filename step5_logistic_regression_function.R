logistic_regression_analysis <- function(data, phenotype_col, snp_list, covariates, output_file) {
  # Ensure PHENOTYPE column is a factor
  data[[phenotype_col]] <- as.factor(data[[phenotype_col]])

  # Initialize an empty data frame to store results
  results <- data.frame(Variant = character(), Beta = numeric(), SE = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

  # Loop through each SNP and fit logistic regression
  for (snp in snp_list) {
    formula <- as.formula(paste(phenotype_col, "~", snp, "+", paste(covariates, collapse = " + ")))
    model <- glm(formula, data = data, family = binomial)
    summary_model <- summary(model)

    # Extract coefficients
    beta <- summary_model$coefficients[2, "Estimate"]
    se <- summary_model$coefficients[2, "Std. Error"]
    p_value <- summary_model$coefficients[2, "Pr(>|z|)"]

    # Store results
    results <- rbind(results, data.frame(Variant = snp, Beta = beta, SE = se, P_value = p_value))
  }

  # Save results to file
  fwrite(results, output_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

  # Return the results
  return(results)
}
