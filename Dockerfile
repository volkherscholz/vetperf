FROM rocker/rstudio

# Install R packages
RUN install2.r --error \
    ggplot2 \
    data.table \
    lme4 \
    knitr \
    glmnet \
    roxygen2 \
    devtools
