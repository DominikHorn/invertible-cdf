#!/bin/bash

TARGET_NAME="icdf_books_experiment"
MODE="DEBUG"
LOWERCASE_MODE=$(echo "${MODE}" | awk '{print tolower($0)}')

./build.sh ${TARGET_NAME} ${MODE}
cmake-build-${LOWERCASE_MODE}/src/${TARGET_NAME}
