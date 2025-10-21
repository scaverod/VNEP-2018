#!/bin/bash

# CPLEX Configuration
# Update these paths to match your CPLEX installation
CPLEX_ROOT="/opt/ibm/ILOG/CPLEX_Studio2211"
CPLEX_INCLUDE="${CPLEX_ROOT}/cplex/include"
CONCERT_INCLUDE="${CPLEX_ROOT}/concert/include"
CPLEX_LIB="${CPLEX_ROOT}/cplex/lib/x86-64_linux/static_pic/libcplex.a"
ILOCPLEX_LIB="${CPLEX_ROOT}/cplex/lib/x86-64_linux/static_pic/libilocplex.a"
CONCERT_LIB="${CPLEX_ROOT}/concert/lib/x86-64_linux/static_pic/libconcert.a"

# Clean old build
rm -rf build

# Configure with CMake
cmake -S . -B build \
  -DCPLEX_INCLUDE_DIR="${CPLEX_INCLUDE}" \
  -DCONCERT_INCLUDE_DIR="${CONCERT_INCLUDE}" \
  -DCPLEX_LIBRARY="${CPLEX_LIB}" \
  -DILOCPLEX_LIBRARY="${ILOCPLEX_LIB}" \
  -DCONCERT_LIBRARY="${CONCERT_LIB}"

# Build
if [ $? -eq 0 ]; then
    echo "Configuration successful, starting build..."
    cmake --build build -j$(nproc)
else
    echo "Configuration failed!"
    exit 1
fi
