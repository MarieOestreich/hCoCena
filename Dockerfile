FROM bioconductor/bioconductor_docker:RELEASE_3_18

COPY ./STAR_protocol/ /home/rstudio/STAR_protocol/
COPY ./reference_files/ /home/rstudio/reference_files/

ADD install_hcocena.R /tmp/
RUN R -f /tmp/install_hcocena.R
