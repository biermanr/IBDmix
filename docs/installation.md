# Installation

## Prerequisites

IBDmix requires:
- CMake 3.14 or higher
- C++ compiler with C++11 support
- Git

## Building from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/IBDmix.git
   cd IBDmix
   ```

2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

3. Configure with CMake:
   ```bash
   cmake ..
   ```

4. Build:
   ```bash
   cmake --build .
   ```

5. The executables will be available in `build/src/`:
   - `ibdmix` - Main IBD detection tool
   - `generate_gt` - Genotype file generator

## Testing the Installation

You can verify the installation by running the unit tests:

```bash
cd build
ctest
```

## Dependencies

The following dependencies are automatically handled by CMake:
- Standard C++ libraries
- Any additional dependencies will be listed here
