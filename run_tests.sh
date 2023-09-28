#!/bin/bash

#
# Runs all insilicoSV tests.
#

set -eux -o pipefail

for TEST_GROUP in processing simulate known_svs
do
    PYTHONPATH=$PWD python test/test_${TEST_GROUP}.py
done
