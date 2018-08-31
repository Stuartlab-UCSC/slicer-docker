FROM rocker/tidyverse

COPY ./traj-converters /home/traj-converters

COPY ./src /home/src

ENV PATH="/home/src:${PATH}"

RUN R -e 'install.packages(c("lle", "SLICER", "optparse","gam"))'
