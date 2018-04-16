#!/bin/bash

# pre-computed checksum
MD5SUM="3cb8404a12d70d63b7e614ec0b02c339" # Obtained with R vs 3.4.1 (same as in Travis)

# run mlm on sample dataset
R CMD INSTALL --preclean --no-multiarch --with-keep.source ../mlm
R -e "library(mlm); res <- mlm(biomarkers ~ .^2, data = patients); save(res, file = 'ci/res.rda')"
[[ $(md5sum ci/res.rda | awk '{$0=$1}1') == $MD5SUM ]]
