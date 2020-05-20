FROM bimberlab/oosap

RUN Rscript -e "devtools::install_github(repo = 'bimberlabinternal/ShinySDA', dependencies = T, upgrade = 'always')" \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds


                
# select port
EXPOSE 3838
