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
Rscript --no-init-file -e "library(knitr); knit('$2.Rmd')"

printf "\n  *** done! *** \n\n"
