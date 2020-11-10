#!/bin/bash

Rscript --vanilla -e "library(mlm); res <- mlm(biomarkers ~ .^2, data = patients); print.default(res)" > ci/res

if [[ $(diff ci/res ci/res0) ]]; then exit 1; fi
