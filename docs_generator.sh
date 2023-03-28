#!/bin/bash

MODULE_PATHS=src/DrillPyze
SOURCE_DIR=docs/source
BUILD_DIR=docs/build


# Generate source files from module path
sphinx-apidoc -f -o $SOURCE_DIR $MODULE_PATH

# Generate html documentation from source files
sphinx-build -b html $SOURCE_DIR $BUILD_DIR