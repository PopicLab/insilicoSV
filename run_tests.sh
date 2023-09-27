#!/bin/bash

#
# Runs all insilicoSV tests.
#

set -eux -o pipefail

for TEST_GROUP in processing simulate known_svs
do
    python -m unittest test.test_${TEST_GROUP}
done
