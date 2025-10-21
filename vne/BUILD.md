# Building VNE with CPLEX

## Prerequisites

1. CPLEX Studio installed (tested with CPLEX_Studio2211)
2. Boost libraries (version 1.54.0 or higher)
3. CMake (version 2.7 or higher)
4. C++ compiler with C++11 support

## Build Instructions

### Option 1: Using the build script (Recommended)

1. Edit `build_with_cplex.sh` and update the CPLEX paths to match your installation
2. Make the script executable:
   ```bash
   chmod +x build_with_cplex.sh
   ```
3. Run the script:
   ```bash
   ./build_with_cplex.sh
   ```

### Option 2: Manual CMake configuration

```bash
cmake -S . -B build \
  -DCPLEX_INCLUDE_DIR="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include" \
  -DCONCERT_INCLUDE_DIR="/opt/ibm/ILOG/CPLEX_Studio2211/concert/include" \
  -DCPLEX_LIBRARY="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/libcplex.a" \
  -DILOCPLEX_LIBRARY="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/libilocplex.a" \
  -DCONCERT_LIBRARY="/opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic/libconcert.a"

cmake --build build -j
```

**Note:** Adjust the paths according to your CPLEX installation location and architecture.

## Common Issues

### Issue: `ilcplex/ilocplex.h: No such file or directory`

This means the CPLEX include directories are not being found. Make sure:
- The CPLEX_INCLUDE_DIR points to the directory containing the `ilcplex` folder
- The CONCERT_INCLUDE_DIR points to the directory containing the `ilconcert` folder
- Check the CMake output for the detected paths

### Issue: `No rule to make target 'libilocplex.a'`

The `libilocplex.a` file might not exist in some CPLEX installations. The CMake configuration has been updated to handle this case - it will use only `libcplex.a` if `libilocplex.a` is not found.

### Issue: `std::iota` or `std::accumulate` not found

This has been fixed by adding `#include <numeric>` to the required source files.

## Files Modified

The following fixes have been applied to make the code compile with modern C++ compilers and CPLEX:

1. **CMake configuration** (`cmake/Findcplex.cmake`): Made ILOCPLEX library optional
2. **Include headers**: Added `#include <numeric>` to:
   - `src/ac.cpp`
   - `src/bprice.cpp`
   - `src/vnedefs.cpp`
   - `include/backtrack.h`
3. **CMakeLists.txt**: Added debug output and improved include directory handling
