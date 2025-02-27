library(data.table)
diseases = c('Jo1','JDM','DM','PM')
for (disease in diseases){
  t2 = fread('c4_isotype_t_dosages_info.csv')
  t1 = fread(paste0('/mount/ictr1/Users/catheriz/TOPMed_liftover_myositis_dmpmall/TOPMED_QC_38_REMOVE_OUTLIER/update/',disease,'_DMPM_Pheno_Plink_PCA_rm_dis_rm_first_degree_rm_dup.sample'))
  t1 = t1[-1,]
  colnames(t2)[1]="ID_1"
  merged_data <- merge(t1, t2, by = "ID_1", all.x = TRUE)
  merged_data$PHENOTYPE = as.factor(merged_data$PHENOTYPE)
  for (i in 6:15){
    merged_data[[i]]=as.numeric(merged_data[[i]])
  }

  # Initialize an empty data frame to store the results
  results <- data.frame(Variant = character(), Beta = numeric(), SE = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

  # Loop through each SNP and fit the model
  for (snp in colnames(merged_data)[16:19]) {
    formula <- as.formula(paste("PHENOTYPE ~", snp, "+ PC1 + PC2 + PC3 + PC4 + PC5"))
    model <- glm(formula, data = merged_data, family = binomial)
    summary_model <- summary(model)

    # Extract coefficients and p-values
    beta <- summary_model$coefficients[2, "Estimate"]
    se <- summary_model$coefficients[2, "Std. Error"]
    p_value <- summary_model$coefficients[2, "Pr(>|z|)"]

    # Append the results
    results <- rbind(results, data.frame(Variant = snp, Beta = beta, SE = se, P_value = p_value))
  }

  full_formula_joint <- as.formula("PHENOTYPE ~ C4A + C4B + PC1 + PC2 + PC3 + PC4 + PC5")

  # Fit the full model including C4A and C4B
  full_model_joint <- glm(full_formula_joint, data = merged_data, family = binomial)
  summary_full_model_joint <- summary(full_model_joint)

  # Extract coefficients and standard errors for C4A and C4B from the full model
  beta_C4A <- summary_full_model_joint$coefficients["C4A", "Estimate"]
  se_C4A <- summary_full_model_joint$coefficients["C4A", "Std. Error"]
  beta_C4B <- summary_full_model_joint$coefficients["C4B", "Estimate"]
  se_C4B <- summary_full_model_joint$coefficients["C4B", "Std. Error"]

  # Define the reduced model formula excluding C4A and C4B
  reduced_formula_joint <- as.formula("PHENOTYPE ~ PC1 + PC2 + PC3 + PC4 + PC5")

  # Fit the reduced model
  reduced_model_joint <- glm(reduced_formula_joint, data = merged_data, family = binomial)

  # Perform a likelihood ratio test for the joint effect of C4A and C4B
  lr_test_joint <- anova(reduced_model_joint, full_model_joint, test = "LRT")

  # Extract the p-value from the likelihood ratio test
  joint_p_value <- lr_test_joint[2, "Pr(>Chi)"]

  # Append the individual effect results for C4A and C4B
  results <- rbind(results, data.frame(Variant = "C4A", Beta = beta_C4A, SE = se_C4A, P_value = summary_full_model_joint$coefficients["C4A", "Pr(>|z|)"]))
  results <- rbind(results, data.frame(Variant = "C4B", Beta = beta_C4B, SE = se_C4B, P_value = summary_full_model_joint$coefficients["C4B", "Pr(>|z|)"]))

  # Append the joint effect result to the results data frame
  results <- rbind(results, data.frame(Variant = "C4A + C4B", Beta = NA, SE = NA, P_value = joint_p_value))

  # Test the interaction effect of C4A and C4B

  # Define the full model formula including C4A, C4B, and their interaction
  full_formula_interaction <- as.formula("PHENOTYPE ~ C4A * C4B + PC1 + PC2 + PC3 + PC4 + PC5")

  # Fit the full model including the interaction term
  full_model_interaction <- glm(full_formula_interaction, data = merged_data, family = binomial)
  summary_full_model_interaction <- summary(full_model_interaction)

  # Extract coefficients and standard errors for the interaction term
  beta_interaction <- summary_full_model_interaction$coefficients["C4A:C4B", "Estimate"]
  se_interaction <- summary_full_model_interaction$coefficients["C4A:C4B", "Std. Error"]
  p_value_interaction <- summary_full_model_interaction$coefficients["C4A:C4B", "Pr(>|z|)"]

  # Append the interaction effect result to the results data frame
  results <- rbind(results, data.frame(Variant = "C4A:C4B Interaction", Beta = beta_interaction, SE = se_interaction, P_value = p_value_interaction))

  # Print the results
  print(results)

  write.table(results,paste0(disease,'_C4_logistic_regression_results.txt'),quote=F,col.names=T,row.names=F,sep='\t')
  
}
