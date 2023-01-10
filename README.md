

# ShinySDA (SiSDA)

A sister-companion shiny app for SDA to load, browse, and process SDA objects. 

The goal of this app is to make a reproducible processing framework of SDA results including QC and component selection.

In this version we use Seurat data structure to store the results and customize their functions for various needs herein.

## Usage

This app expects the following hierarchy in the dir path it is looking: 

MyExperiment > sda_results > MySDAOuts

The MySdaOuts is the path that has the it#### folders. Also tSNE and other shinySDA analysis are saved here.

sda_results is a folder upstream which stores or or several replicates of the same SDA run with a common Meta data. 

MyExperiment is a folder upstream of sda_results. It should have a dataframe (with the right rownames) saverd as _MetaDF.rds in this folder. Also the SDA input files (the matrix file and _dimnames.rds) need to be here. 

## Version Update

### 1.1.0

Several new tabs, figures, and features have been implemented. Briefly, 
A) Load SDA files either directly from Prime-seq using a numerical ID or traditionally, directly from a path.
B) Meta data for prime-seq data automatically is downloaded and loaded in or traditionally, save a data frame as an _MetaDF.rds in the path
C) New tabs and figures.

### 1.0.0

Working stable release that works for loading SDA runs from a path

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


## Install : 


    devtools::install_github(repo = 'bimberlabinternal/ShinySDA', dependencies = T, upgrade = 'always')

## Launch : 

    ShinySDA::launchShinySDA()
  
 



