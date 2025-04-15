This is a pipeline for genotype matrix preprocessing and gene association testing of GBR sampleset with control cohort of Europeans from the DNASCoRe (https://dnascore.net/) platform.
# Instruction of using the script:
1) Follow the pipeline (the scheme is listed below)
   ![image](https://github.com/user-attachments/assets/aa11fb56-e0e5-44d7-8b14-ec523fbf892a)

You should also make a binary file, you can do it with the SVDFunctions library (https://github.com/alexloboda/SVDFunctions).
And you should make a /yml file of your samples of interest, you can just follow the tutorial on the plantofrm (ttps://dnascore.net/).
3) In the result you will have: table of all variants you have in your data with annotation, and a binary file of your needed subpopulation.
4) From this step you should use the script (file "Data_processing_GBR.Rmd"). The parameteres described there can be changed according to your wish, depending on what allele frequency you are looking for.
In the result you should get a qqplot that looks like this, if everything went well:
![image](https://github.com/user-attachments/assets/4669ce7e-073a-4187-b677-89fa165b48ee)

