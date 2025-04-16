This is a pipeline for genotype matrix preprocessing and gene association testing of GBR sampleset with control cohort of Europeans from the DNASCoRe (https://dnascore.net/) platform.
# Instruction of using the script:
Follow the pipeline (the scheme is listed below)
## Prepare the data
   ![image](https://github.com/user-attachments/assets/aa11fb56-e0e5-44d7-8b14-ec523fbf892a)
You should also make a binary file, you can do it with the SVDFunctions library (https://github.com/alexloboda/SVDFunctions) (an example of metadata for a binary file for GBR population can be found at this repo, named "new_bin_meta").
## .YAML file for the platform
And you should make a .yml file of your samples of interest, you can just follow the tutorial on the plantofrm (ttps://dnascore.net/).
# Files to work with using this code
In the result you will have: table of all variants you have in your data with annotation, and a binary file of your needed subpopulation .
   ![image](https://github.com/user-attachments/assets/8bfb77ce-a549-45fd-b3ad-36e7e0cfc5d0)
## Using the filtration script 
5) From this step you should use the script (file "Data_processing_GBR.Rmd"). The parameteres described there can be changed according to your wish, depending on what allele frequency you are looking for.
When you have the genotype matrix for your case group and you want to get control samples from the platfrom, you should take the same setting as you had while filtering your case data. For example, if you took synonymous variants with gnomAD threshold of 1% and internal frequency threshold of 1%, you should do the same on the platfrom. Alpha and std cannot be regulated on the platform, so you do not tune it.
In the result you should get a qqplot that looks like this, if everything went well:
![image](https://github.com/user-attachments/assets/4669ce7e-073a-4187-b677-89fa165b48ee)

