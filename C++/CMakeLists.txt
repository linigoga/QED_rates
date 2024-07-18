cmake_minimum_required(VERSION 3.4...3.24)
project(QEDProcesses)

# Find Python (using the specified interpreter)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
set(Python_EXECUTABLE "/Users/lucasinigogamiz/anaconda3/bin/python") # Match your compile.sh

# Add pybind11
find_package(pybind11 REQUIRED)

# Additional compiler and linker flags from your .sh script
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -Wall -std=c++20 -fPIC")

# Get Python include directories and libraries (similar to your .sh script)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
    OUTPUT_VARIABLE PYTHON_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND "${Python_CONFIG_PATH}" --includes
    OUTPUT_VARIABLE PYTHON_INCLUDES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Include the Python libraries
link_directories("${PYTHON_LIB_DIR}")

# Add your source files and link with pybind11 and Python
add_library(QEDProcesses MODULE QEDProcesses.cpp binding.cpp)
target_include_directories(QEDProcesses PRIVATE ${pybind11_INCLUDE_DIRS} ${PYTHON_INCLUDES})
target_link_libraries(QEDProcesses PRIVATE pybind11::module ${PYTHON_LIBRARIES})

# Set the output directory and module name (similar to before)
set_target_properties(QEDProcesses PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"  # Put module in build dir
    PREFIX ""  # Remove any default prefix from the module name
)
