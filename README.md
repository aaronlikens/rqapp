# rqapp

This is an R package that performs Recurrence Quantification Analyis on univariate and multivariate time series data. We have included the capability of performing cross, joint, and multidimensional variants of recurrence quantification analysis. Also included are functions for performing windowed versions of these tools. The goal of the project was to create fast, memory efficient implementation of RQA in R that mimics functionality of the excellent 'crqa' package by Coco and Dale (2014). The package is fully functional but still a work in progress. We will continue to clean up the code, optimize some of the slower functions, add some new ones, and develop routines suitable for cluster computing. 

Reference:
Coco, M. I., & Dale, R. (2014). Cross-recurrence quantification analysis of categorical and continuous time series: an R package. Frontiers in psychology, 5, 510.

# installation
rqapp can now be installed from github:

library(devtools)

install_github('aaronlikens/rqapp')

NOTE: This is a source package and requires compilation of C++ code. 
Windows users can install RTools software package to get necessary components:
https://cran.r-project.org/bin/windows/Rtools/

Mac users should consider installing the tools found at:
https://cran.r-project.org/bin/macosx/tools/
