#!/bin/bash

############################################
#
# usage:
#
# ./knit.sh course1 random_variables
#
############################################

printf "\n  *** knit *** \n\n"

cd ../labs/$1
/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/bin/Rscript --no-init-file -e "library(knitr); knit('$2.Rmd')"

printf "\n  *** done! *** \n\n"
