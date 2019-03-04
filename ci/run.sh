#!/bin/bash

# Pre-computed checksum obtained using Docker and R vs 3.5.2 (as in Travis)

MD5SUM="12a83357254c789197b9076cfdbc8b63"

# Install and run mlm on sample dataset

COMM1='R CMD INSTALL --preclean --no-multiarch --with-keep.source ../mlm --library=ci/R/lib'
COMM2='R --vanilla -e "library(mlm, lib = \"ci/R/lib\"); res <- mlm(biomarkers ~ .^2, data = patients); save(res, file = \"ci/res.rda\")"'

docker run --rm --entrypoint /bin/bash -w $PWD -v $PWD:$PWD dgarrimar/mlm:r3.5.2 -c "$COMM1;$COMM2"

[[ $(md5sum ci/res.rda | awk '{$0=$1}1') == $MD5SUM ]]
