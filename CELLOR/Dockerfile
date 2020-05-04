FROM rocker/tidyverse:3.6.1
RUN apt-get update && apt-get install -y --no-install-recommends libcurl4-openssl-dev libproj-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev
RUN R -e "library(BiocManager); BiocManager::install('graph')"
RUN R -e "install.packages(c('reshape2','gridExtra','ggtern','ggalt','igraph'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
COPY . /home/CELLOR
WORKDIR /home/CELLOR
CMD ["R"]
