# WassersteinGoF

## Usage 

This package helps to perform goodness-of-fit tests relying on the Wasserstein (p=1,2) distance. 
The functions have been used in https://arxiv.org/pdf/2003.06684.pdf. 

Have a look at the Temp/Examples.R for a first example of certain functions to test the null hypothesis that the data is Gaussian.


## Installation on Unix/ macOS

To install the package from R you must use the two following lines.
```
library(devtools)

install_github("gmordant/WassersteinGoF", ref = "main")
```
Somtimes, you need to update some packages in the process. Accept it and restart the R session before trying to use 'WassersteinGoF'.
 

Note: This package was succesfully tested on both MacOS and Ubuntu (version 18 and 20).

## Installation on Windows

Download the zip containing the package and unzip it. Open the project by clicking on 'WassersteinGoF.Rproj' file. 
Go to R/Functions and remove the two last functions. The latter call a parallelisation package unsupported by Windows. 
You can now click on 'Install and Restart' from the Build Pane. The first two examples will still work.

## Further

Do not hesitate to reach out for support or comments.
