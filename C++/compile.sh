#!/bin/bash

PYTHON_BIN_PATH="/Users/linigoga/anaconda3/bin/python"
PYTHON_CONFIG_PATH="/Users/linigoga/anaconda3/bin/python3.10-config"
PYTHON_LIB_DIR=$($PYTHON_BIN_PATH -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
#PYBIND_INCLUDE_DIR="/Users/linigoga/anaconda3/lib/python3.10/site-packages/pybind11/include/pybind11"
PYBIND_INCLUDE_DIR="/Users/linigoga/anaconda3/lib/python3.10/site-packages/pybind11/include"


g++ -O0 -g -Wall -shared -std=c++20 -fPIC \
-I$PYBIND_INCLUDE_DIR \
-L$PYTHON_LIB_DIR \
`$PYTHON_CONFIG_PATH --includes` \
-Wl,-rpath,$PYTHON_LIB_DIR \
-undefined dynamic_lookup \
binding.cpp QEDProcesses.cpp\
 -o QEDProcesses`$PYTHON_CONFIG_PATH --extension-suffix`


