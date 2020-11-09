#!/bin/bash

# Pre-computed checksum obtained using Docker and R vs 4.0.3 (as in Travis)

MD5SUM="b78a1116c530ae5a05e5218ba4e1ce84"

# Install and run mlm on sample dataset

COMM1='R CMD INSTALL --preclean --no-multiarch --with-keep.source ../mlm'
COMM2='R --vanilla -e "library(mlm); res <- mlm(biomarkers ~ .^2, data = patients); save(res, file = \"ci/res.rda\")"'

docker run --rm --entrypoint /bin/bash -w $PWD -v $PWD:$PWD dgarrimar/mlm:r4.0.3 -c "$COMM1;$COMM2"

[[ $(md5sum ci/res.rda | awk '{$0=$1}1') == $MD5SUM ]]
