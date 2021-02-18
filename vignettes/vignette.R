## ---- fig.show='hold', message=FALSE------------------------------------------

### Removing any previous versions of the package
#First, ensure there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Downloading and installing Platypus

# Download most recent version from master branch at https://github.com/alexyermanos/Platypus We can install the package using the following command. 
# WARNING: This needs to be replaced with your own directory where the downloaded package is found

# For MacOS users it may look like this:
install.packages("~/Downloads/Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# For windows it will likely look something like this: 
# WARNING: You will need to replace 'YourPCName' with your user name for the windows account in the directory. 
# install.packages("C:\Users\YourPCName\Downloads\Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# Now we can load the installed package into the R environment. In case of problems with installing other R packages that are used in Platypus, please see the README file at the https://github.com/alexyermanos/Platypus, where we outline how to install the other R packages for both Windows and MacOS.
library(Platypus)

# Individual R functions can additionally be found on the github Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment in case not all functions are desired. These functions are actively updated and may include more features than the in original tar.gz file. 


### Downloading the test data
# The COVID-19 data (~136 MB size of the zip file) can be found at the following link https://polybox.ethz.ch/index.php/s/fxQJ3NrRSwiPSSo This dataset contains VDJ (separate libraries for B and T cells) and GEX libraries from two convalescent COVID-19 patients.

# After downloading the zip file named "PlatypusTestData.zip", please unzip the file and find the path to the newly formed folder. Typically this will be in the Downloads folder, so the code below should work on MacOS. For Windows please uncomment the code and change the user name to match your PC.

directory_to_covid_patients_gex <- list()
directory_to_covid_patients_gex[[1]] <- c("~/Downloads/PlatypusTestData/Patient1_GEX/")
directory_to_covid_patients_gex[[2]] <- c("~/Downloads/PlatypusTestData/Patient2_GEX/")

# For Windows: 
#directory_to_covid_patients_gex[[1]] <- c("C:\Users\YourPCName\Downloads\PlatypusTestData\Patient1_GEX")
#directory_to_covid_patients_gex[[2]] <- c("C:\Users\YourPCName\Downloads\PlatypusTestData\Patient2_GEX")



