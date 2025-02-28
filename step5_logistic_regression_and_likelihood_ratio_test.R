# Function to perform logistic regression for SNPs with user-defined covariates
run_logistic_regression <- function(data, phenotype_col, snp_list, covariates, output_file) {
  # Ensure phenotype is a factor
  data[[phenotype_col]] <- as.factor(data[[phenotype_col]])

  # Initialize results dataframe
  results <- data.frame(
    Variant = character(), 
    Beta = numeric(), 
    SE = numeric(), 
    P_value = numeric(), 
    stringsAsFactors = FALSE
  )

  # Loop through each SNP and fit logistic regression
  for (snp in snp_list) {
    full_formula <- as.formula(paste(phenotype_col, "~", snp, "+", paste(covariates, collapse = " + ")))
    model <- glm(full_formula, data = data, family = binomial)
    summary_model <- summary(model)

    # Extract coefficients and p-values
    beta <- summary_model$coefficients[2, "Estimate"]
    se <- summary_model$coefficients[2, "Std. Error"]
    p_value <- summary_model$coefficients[2, "Pr(>|z|)"]

    # Store results
    results <- rbind(results, data.frame(Variant = snp, Beta = beta, SE = se, P_value = p_value))
  }

  # Save results
  fwrite(results, output_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

  return(results)
}

# Function to perform Likelihood Ratio Test (LRT) for joint effects of SNPs
run_likelihood_ratio_test <- function(data, phenotype_col, snp_list, covariates) {
  # Ensure phenotype is a factor
  data[[phenotype_col]] <- as.factor(data[[phenotype_col]])

  # Define full model (with SNPs)
  full_formula <- as.formula(paste(phenotype_col, "~", paste(c(snp_list, covariates), collapse = " + ")))
  full_model <- glm(full_formula, data = data, family = binomial)

  # Define reduced model (without SNPs)
  reduced_formula <- as.formula(paste(phenotype_col, "~", paste(covariates, collapse = " + ")))
  reduced_model <- glm(reduced_formula, data = data, family = binomial)

  # Perform Likelihood Ratio Test
  lr_test <- anova(reduced_model, full_model, test = "LRT")
  joint_p_value <- lr_test[2, "Pr(>Chi)"]

  return(data.frame(Variant = paste(snp_list, collapse = " + "), Beta = NA, SE = NA, P_value = joint_p_value))
}
