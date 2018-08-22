FROM rocker/tidyverse

COPY ./traj-converters /home/traj-converters

RUN R -e 'install.packages(c("lle", "SLICER", "optparse","gam"))'
