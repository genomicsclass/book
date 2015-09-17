---
layout: page
title: Finding more help for data analysis
---

Here are some resources for finding more help when performing data analysis for genomics.

### SEQanswers and Biostars

These are both very useful for getting feedback on genomics experiments and getting data analysis tips.

- [SEQanswers](http://seqanswers.com/) is "an information resource and user-driven community focused on all aspects of next-generation genomics".  They say "The site will always attempt to cater to everyone, regardless of scientific background or knowledge."  SEQanswers also has a wiki listing many various [software packages for genomics data analysis](http://seqanswers.com/wiki/Software).

- [Biostars](https://www.biostars.org/) is a website focused on "bioinformatics, computational genomics and biological data analysis" They say: "No question is too trivial or too 'newbie'."

### Bioconductor support site

If your question is about Bioconductor software, you can post a question to the 
[Bioconductor support site](https://support.bioconductor.org).  This
is answered by 100s of statisticians, bioinformaticians and developers
around the world. See the
[posting guide](http://www.bioconductor.org/help/support/posting-guide/)
for instructions on how to optimally write a post for getting the best feedback.
The support site is searchable and extends back to the
[first post ever](https://support.bioconductor.org/p/2/) on 26
November, 2001.  If you directly email a Bioconductor developer with a
support question, they will most likely send an auto-response for you
to post to the support site.

In addition to the Posting Guide, there are four **must-do's** for submitting an email to the Bioconductor mailing list:

1. **Read the help pages and vignette**, especially for the function
   you are using: `?functionName`.  While not all questions are
   answered in the vignette, *many* are. Package maintainers can get
   grumpy (they are human) if you ask a question which is clearly
   answered in the vignette or the help page. Remember, typing
   `browseVignettes("pkgName")` will bring you to the software vignettes.
   More on [how to find help within R.](installing_Bioconductor_finding_help.html)
2. **Explain what you are trying to do**. If the question is not
   simply reporting a bug, i.e., you are not sure exactly what kind of
   analysis to run, then you should explain what biological question
   you are trying to answer with your experiment/data, and what the
   experimental design is.  You don't have to give away secrets, e.g.,
   "we are interested in looking at differences in protein binding for
   2 treatments and control, with 3 biological replicates in each
   condition." 
3. **Always paste the output of sessionInfo()** at the bottom of your
   email. This gives the package maintainers lots of useful
   information, like what version of R you are running, your platform,
   and the versions of all the packages. 
4. **Paste in all of your code**, and any errors or warnings if they
   occurred. 

If you have a general R question, there is the
[R-help mailing list](https://stat.ethz.ch/mailman/listinfo/r-help). 

