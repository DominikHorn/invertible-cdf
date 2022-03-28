#!/bin/bash

TARGET_NAME="icdf_books_experiment"
./build.sh ${TARGET_NAME} RELEASE
cmake-build-release/src/${TARGET_NAME}
