This is a pipeline for genotype matrix preprocessing and gene association testing of GBR sampleset with control cohort of Europeans from the DNASCoRe (https://dnascore.net/) platform.
# Instruction of using the script:
1) Follow the pipeline (the scheme is listed below)
   ![image](https://github.com/user-attachments/assets/01a9e02b-74ed-4a5e-a23e-7162a09a3266)
You should also make a binary file, you can do it with the SVDFunctions library (https://github.com/alexloboda/SVDFunctions).
And you should make a /yml file of your samples of interest, you can just follow the tutorial on the plantofrm (ttps://dnascore.net/).
2) In the result you will have: table of all variants you have in your data with annotation, and a binary file of your needed subpopulation.
3) From this step you should use the script (file "Data_processing_GBR.Rmd"). The parameteres described there can be changed according to your wish, depending on what allele frequency you are looking for.
