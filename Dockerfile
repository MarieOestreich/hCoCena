FROM bioconductor/bioconductor_docker:3.18-R-4.3.2

COPY ./STAR_protocol/ /home/rstudio/STAR_protocol/
COPY ./reference_files/ /home/rstudio/reference_files/

ADD install_hcocena.R /tmp/
RUN R -f /tmp/install_hcocena.R