FROM alpine:latest
COPY . .
COPY --from=staphb/kma /usr/local/bin/* /usr/local/bin/
COPY --from=bschiffthaler/diamond /usr/local/bin/diamond /usr/local/bin/diamond

RUN apk update

ENV PYTHONUNBUFFERED=1
RUN apk add --update --no-cache python3 && ln -sf python3 /usr/bin/python

RUN apk add --update R R-dev wget curl alpine-sdk libc6-compat gcompat

RUN R -e "install.packages('parallel',repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager',repos = 'http://cran.us.r-project.org')"
# Use the Bioconductor release compatible with R, requiring at least 3.22.
RUN R -e "if (BiocManager::version() < numeric_version('3.22')) stop('Bioconductor >= 3.22 is required; upgrade R'); BiocManager::install(c('pwalign', 'Biostrings'), ask=FALSE)"

LABEL version="0.3"
LABEL maintainer="Cailean Carter"
LABEL email="cailean.carter@quadram.ac.uk"
