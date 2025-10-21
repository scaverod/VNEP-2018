# VNEP-2018 Build Fixes Summary

## Issues Fixed

### 1. CMake not respecting command-line CPLEX variables
**Problem:** When you passed `-DCPLEX_INCLUDE_DIR=...` to cmake, it was being ignored because the `Findcplex.cmake` script was unconditionally searching for paths.

**Solution:** Modified `vne/cmake/Findcplex.cmake` to only search for CPLEX paths if they haven't already been provided via command-line arguments. Now all `FIND_PATH` and `FIND_LIBRARY` commands are wrapped in `IF(NOT VARIABLE_NAME)` checks.

### 2. Missing ILOCPLEX library
**Problem:** The file `libilocplex.a` doesn't exist in some CPLEX installations.

**Solution:** Modified `vne/cmake/Findcplex.cmake` to make `ILOCPLEX_LIBRARY` optional - the build will work with just `libcplex.a` if `libilocplex.a` is not available.

### 3. Missing `#include <numeric>` header
**Problem:** Compilation errors for `std::iota` and `std::accumulate` - these require the `<numeric>` header in C++11.

**Solution:** Added `#include <numeric>` to:
- `vne/src/ac.cpp`
- `vne/src/bprice.cpp`
- `vne/src/vnedefs.cpp`
- `vne/include/backtrack.h`

### 4. CMakeLists.txt improvements
**Problem:** No visibility into what CPLEX paths were being detected.

**Solution:** 
- Added debug output messages to show detected CPLEX paths
- Improved include directory handling to use `CPLEX_INCLUDE_DIRS` when available

## How to Build

### Quick Start (Using the build script)

```bash
cd vne
chmod +x build_with_cplex.sh
# Edit build_with_cplex.sh to set your CPLEX paths
./build_with_cplex.sh
```

### Manual Build

```bash
cd vne

# Clean previous build
rm -rf build

# Configure (adjust paths for your system)
cmake -S . -B build \
  -DCPLEX_INCLUDE_DIR="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include" \
  -DCONCERT_INCLUDE_DIR="/opt/ibm/ILOG/CPLEX_Studio2211/concert/include" \
  -DCPLEX_LIBRARY="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/libcplex.a" \
  -DCONCERT_LIBRARY="/opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic/libconcert.a"

# Note: ILOCPLEX_LIBRARY is now optional - you can omit it if the file doesn't exist

# Build
cmake --build build -j
```

### What to Look for in CMake Output

After running cmake configure, you should see:

```
-- Found CPLEX: /opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/libcplex.a
-- CPLEX_INCLUDE_DIR: /opt/ibm/ILOG/CPLEX_Studio2211/cplex/include
-- CONCERT_INCLUDE_DIR: /opt/ibm/ILOG/CPLEX_Studio2211/concert/include
-- CPLEX_INCLUDE_DIRS: /opt/ibm/ILOG/CPLEX_Studio2211/cplex/include;/opt/ibm/ILOG/CPLEX_Studio2211/concert/include;
```

If these paths are empty or wrong, the build will fail with "ilcplex/ilocplex.h: No such file or directory".

## Testing the Fixes

After applying these changes and rebuilding from scratch (deleting the `build` directory first), the compilation should succeed without the previous errors.

## Files Modified

1. `vne/cmake/Findcplex.cmake` - Made all FIND_PATH/FIND_LIBRARY conditional, made ILOCPLEX optional
2. `vne/CMakeLists.txt` - Added debug output and improved include handling
3. `vne/src/ac.cpp` - Added `#include <numeric>`
4. `vne/src/bprice.cpp` - Added `#include <numeric>`
5. `vne/src/vnedefs.cpp` - Added `#include <numeric>`
6. `vne/include/backtrack.h` - Added `#include <numeric>`

## Files Created

1. `vne/build_with_cplex.sh` - Convenient build script
2. `vne/BUILD.md` - Build instructions
3. `vne/FIXES.md` - This file

## Next Steps

1. Delete your existing `build` directory: `rm -rf vne/build`
2. Run the cmake configure command with the correct CPLEX paths
3. Check that the CMake output shows the correct paths
4. Run the build: `cmake --build build -j`

The build should now complete successfully!
