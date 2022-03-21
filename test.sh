#!/bin/bash

TARGET=icdf_tests

# setup script
set -e
cd "$(dirname "$0")"

# build and run tests in debug mode (to catch issues with address sanitizer etc)
./build.sh ${TARGET} DEBUG
cmake-build-debug/src/${TARGET} $@

# build and run tests in debug mode (to catch potential deviations due to DNDEBUG etc)
./build.sh ${TARGET} RELEASE
cmake-build-release/src/${TARGET} $@
