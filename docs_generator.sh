#!/bin/bash

MODULE_PATHS=src/DrillPyze
SOURCE_DIR=docs/source
BUILD_DIR=docs/build


# Generate source files from module path
sphinx-apidoc -f -o $SOURCE_DIR $MODULE_PATH


# Replace delimiter
#SOURCE_FILES=${MODULE_PATHS/ /\\n   }

# Insert source file name
#sed "s/Contents:/Contents:\n\n   $SOURCE_FILES/g" $SOURCE_DIR/index.tmpl.rst > $SOURCE_DIR/index.rst

# Generate html documentation from source files
sphinx-build -b html $SOURCE_DIR $BUILD_DIR