#!/bin/bash

# Python Configuration (Adapt for your actual Python version if needed)
PYTHON_BIN_PATH="/Users/lucasinigogamiz/anaconda3/bin/python"  
PYTHON_CONFIG_PATH="/Users/lucasinigogamiz/anaconda3/bin/python3.11-config" 

# Project Directories
PROJECT_DIR="/Users/lucasinigogamiz/Documents/GitHub/QED_rates/C++"
BUILD_DIR="$PROJECT_DIR/build"

# Create build directory if it doesn't exist
mkdir -p $BUILD_DIR

# Change to the build directory
cd $BUILD_DIR

# Configure CMake (adjust -DPYTHON_EXECUTABLE if your Python is elsewhere)
cmake -DPYTHON_EXECUTABLE=$PYTHON_BIN_PATH ..

# Build with Make
make