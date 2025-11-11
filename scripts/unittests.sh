#!/bin/bash
# Run the same unit tests that are executed in the GitHub Actions workflow
# This script replicates the steps from .github/workflows/tests.yml

set -e  # Exit on error

# Get the project root directory (parent of scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build"

echo "==> Project root: ${PROJECT_ROOT}"
echo "==> Build directory: ${BUILD_DIR}"
echo ""

# Step 1: Configure CMake
echo "==> Configuring CMake with BUILD_TESTS=ON"
cmake -B "${BUILD_DIR}" -DBUILD_TESTS=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --version
echo ""

# Step 2: Build
echo "==> Building project"
cmake --build "${BUILD_DIR}"
echo ""

# Step 3: Run tests with ctest
echo "==> Running tests with ctest"
cd "${BUILD_DIR}"
ctest --output-on-failure
echo ""

# Step 4: Run valgrind (Linux only)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "==> Checking for valgrind"
    if ! command -v valgrind &> /dev/null; then
        echo "valgrind not found. Install with: sudo apt-get install -y valgrind"
        echo "Skipping valgrind tests."
    else
        echo "==> Running valgrind memory tests"
        ls tests/*_test | \
            xargs -I{} valgrind --error-exitcode=1 --leak-check=yes ./{} --gtest_shuffle --gtest_repeat=10 > /dev/null
        echo "Valgrind tests passed!"
    fi
    echo ""
else
    echo "==> Skipping valgrind (only runs on Linux)"
    echo ""
fi

echo "==> All tests completed successfully!"
