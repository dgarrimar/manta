#!/bin/bash

# pre-computed checksum
MD5SUM="85ec4dd714957756413a834a7a3842f6" # Obtained with R vs 3.4.1 (same as in Travis)

# run mlm2 on sample dataset
R CMD INSTALL --preclean --no-multiarch --with-keep.source ../mlm
R -e "library(mlm); fit <- mlm2(biomarkers ~ .^2, data = patients); save(fit, file = 'ci/fit.rda')"
[[ $(md5sum ci/fit.rda | awk '{$0=$1}1') == $MD5SUM ]]
