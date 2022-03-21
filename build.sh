#!/bin/bash

# Setup script
source .env
set -e
cd "$(dirname "$0")"

# Parse arguments
TARGET=${1:-"icdf_tests"}
BUILD_TYPE=${2:-"RELEASE"}
BUILD_DIR="cmake-build-$(echo "${BUILD_TYPE}" | awk '{print tolower($0)}')/"

# Make sure repo is setup properly
./setup.sh ${BUILD_TYPE}

# Link compile_commands.json
ln -fs ${BUILD_DIR}compile_commands.json compile_commands.json

# Build
cmake \
  --build ${BUILD_DIR} \
  --target ${TARGET} \
  -j
