# Pipeline for C4 Imputation, Dosage Calculation, and Association Testing
This pipeline performs C4 imputation, dosage extraction, isotype calculation, and incorporates an additional R function to run logistic regression and likelihood ratio tests.
## Step 1: C4 Imputation
This step is adopted from [imputec4](https://github.com/freeseek/imputec4).
Run the Beagle imputation pipeline for C4:
```
bash step1_run_imputec4.sh -i input.vcf.gz -o output_prefix
```
### Parameters
 - ```-i```: Input VCF file (required)
 - ```-o```: Output prefix for the imputed files (required)
 - ```-b```: Genome build version (default: 38, options: 37 or 38)
 - ```-j```: Path to Beagle JAR file (default: beagle.25Nov19.28d.jar)
## Step 2: Extract imptued C4 Ddata from Imputed VCF
Run the shell script to extract relevant dosage fields:
```
bash step2_extract_C4.sh -i imputed.vcf.gz
```
### Parameters
 - ```-i```: Input imputed VCF file (required)
 - ```-o```: Output filtered VCF file (default: <input>_DS.vcf.gz)
 - ```-d```: Output dosage file (default: dosages.txt)
## Step 3: Process Dosage Data to Calculate Isotype
Run the R script to process dosage data and calculate isotypes:
```
Rscript step3_process_data_to_calculate_isotype.R <input_file> <output_prefix>
```
### Parameters
- ```<input_file>```: Input dosage file
- ```<output_prefix>```: Prefix for output files
## Step 4: C4 Haplotype Dosage Processor
Run the Python script to decode C4 isotypes:
```
python step4_decode_C4_isotype.py -i <input_file> -o <output_prefix>
```
### Parameters
- ```-i```: Input file name generated after step 3 (required)
- ```-o```: Output prefix for the output file (required)
## Step 5: Logistic Regression and Likelihood Ratio tests
R functions to perform statistical association testing using logistic regression and likelihood ratio tests.

## Citation
```
Zhu, Catherine et al. “Genetic Architecture of Idiopathic Inflammatory Myopathies From Meta-Analyses.”
Arthritis & Rheumatology (Hoboken, N.J.), 10.1002/art.43088. 16 Dec. 2024, doi:10.1002/art.43088.
```
