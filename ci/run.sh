#!/bin/bash

# pre-computed checksum
MD5SUM="642b1da659136e9e8df4052fa594ffce" # Obtained with R vs 3.4.1 (same as in Travis)

# run mlm2 on sample dataset
R CMD INSTALL --preclean --no-multiarch --with-keep.source ../mlm
R -e "library(mlm); fit <- mlm2(biomarkers ~ .^2, data = patients); save(fit, file = 'ci/smm.rda')"
[[ $(md5sum ci/smm.rda | awk '{$0=$1}1') == $MD5SUM ]]
