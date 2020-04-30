Data can be downloaded from AMP-AD links as reference in the  Key Resource Table

The proteomic intensities and RNA seq counts were log2 transformed. Protein LFQ (proteomics) and normalized FPKM values (RNA-seq) were assessed for effects from biological covariates (diagnosis, age, gender) and technical variables (batch, brain bank, etc). We used a linear regression model accounting for biological and technical covariates depending upon the cohort to be analyzed. The final model used was implemented in R version 3.6.1 (R Core Team, 2019) as follows:
lm(expression ~ diagnosis + age + gender + batch + brain.bank.batch)

The protein/RNA data and metadata was save in an Rdata file which contained
-datExpr.*cohort* - expression data - sample (row) x expression (columns)
-targets.*cohort* - metadata - sample (row) x metadata (columns)
