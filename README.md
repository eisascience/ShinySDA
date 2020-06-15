# ShinySDA (SiSDA)

A shiny app to load, browse, and process SDA objects. 

The goal of this app is to make a reproducible processing framework of SDA results including QC and component selection.

In this version we use Seurat data structure to store the results and customize their functions for various needs herein.

## Usage

This app expects the following hierarchy in the dir path it is looking: 

MyExperiment > sda_results > MySDAOuts

The MySdaOuts is the path that has the it#### folders. Also tSNE and other shinySDA analysis are saved here.

sda_results is a folder upstream which stores or or several replicates of the same SDA run with a common Meta data. 

MyExperiment is a folder upstream of sda_results. It should have a dataframe (with the right rownames) saverd as _MetaDF.rds in this folder. Also the SDA input files (the matrix file and _dimnames.rds) need to be here. 




## Refs:

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2019). shiny: Web Application
  Framework for R. R package version 1.4.0. https://CRAN.R-project.org/package=shiny
  
Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R
  package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard
  
Daniel Wells and Victoria Hore (2017). SDAtools: SDAtools: A toolkit for SDA. R package version 1.0.

Stuart and Butler et al. Comprehensive integration of single cell data. bioRxiv (2018).


## Contact: 

Bimber Lab:

Produced by: @eisamahyari
Converted to an R packaged by Kosiso Onwuzu.


## Install : 


    devtools::install_github(repo = 'bimberlabinternal/ShinySDA', dependencies = T, upgrade = 'always')

## Launch : 

    ShinySDA::launchShinySDA()
  
 



