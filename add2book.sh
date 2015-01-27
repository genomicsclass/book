#!/bin/bash

############################################
#
# usage:
#
# ./add2book.sh course1 random_variables
#
############################################

printf "\n  *** make sure BOOK is up-to-date *** \n\n"

git checkout gh-pages
git pull
cd ..

printf "\n  *** copy from LABS to BOOK *** \n\n"

cp labs/$1/$2.md book/pages

imgcount=`ls -1 labs/$1/figure/$2* 2> /dev/null | wc -l`
if [ $imgcount -gt 0 ]
then
cp labs/$1/figure/$2* book/pages/figure
fi

printf "  *** add new files to BOOK, commit and push *** \n\n"

cd book
git add pages/$2.md

if [ $imgcount -gt 0 ]
then
git add pages/figure/$2*
fi

git commit -am "adding $2 to book"
git push origin gh-pages
cd ..

printf "\n  *** done! *** \n\n"
