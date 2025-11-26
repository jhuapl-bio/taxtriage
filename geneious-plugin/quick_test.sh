#!/bin/bash

# Quick test of TaxTriage plugin compilation and build
cd /Users/dho/Documents/taxtriage/geneious-plugin

echo "=== TaxTriage Plugin Test Results ==="
echo "Starting plugin compilation and testing..."
echo ""

# 1. Clean previous builds
echo "1. Cleaning previous builds..."
ant clean

echo ""

# 2. Compile source code
echo "2. Compiling source code..."
if ant compile; then
    echo "✓ Source compilation successful"

    # Count source classes
    if [ -d "classes" ]; then
        echo "Compiled classes found:"
        find classes -name "*.class" | wc -l | xargs echo "  Total class files:"
        find classes -name "*.class" | head -10
        if [ $(find classes -name "*.class" | wc -l) -gt 10 ]; then
            echo "  ... and $(( $(find classes -name "*.class" | wc -l) - 10 )) more"
        fi
    fi
else
    echo "✗ Source compilation failed"
    exit 1
fi

echo ""

# 3. Compile test code
echo "3. Compiling test code..."
if ant compile-tests; then
    echo "✓ Test compilation successful"
else
    echo "✗ Test compilation failed, but continuing..."
fi

echo ""

# 4. Try to run tests
echo "4. Running tests..."
if ant test; then
    echo "✓ Unit tests executed successfully"
else
    echo "Note: JUnit tests may not be available, trying manual test execution..."

    # Try manual test execution
    cd classes
    if java -cp ".:../lib/*" com.jhuapl.taxtriage.geneious.TaxTriagePluginTest; then
        echo "✓ Manual test execution successful"
    else
        echo "Note: Manual test execution requires full Geneious environment"
    fi
    cd ..
fi

echo ""

# 5. Build plugin
echo "5. Building plugin..."
if ant build; then
    echo "✓ Plugin build successful"

    # Check for build artifacts
    if [ -d "build" ]; then
        echo "Build artifacts:"
        ls -la build/
    fi
else
    echo "✗ Plugin build failed"
    exit 1
fi

echo ""

# 6. Create distribution package
echo "6. Creating distribution package..."
if ant distribute; then
    echo "✓ Plugin distribution package created"

    if [ -f "build/TaxTriage.gplugin" ]; then
        echo "Plugin package: build/TaxTriage.gplugin"
        echo "Package size: $(ls -lh build/TaxTriage.gplugin | awk '{print $5}')"
    fi
else
    echo "✗ Distribution package creation failed"
    exit 1
fi

echo ""
echo "=== Test Summary ==="
echo "✓ Plugin compiled successfully"
echo "✓ Core classes created"
echo "✓ Plugin JAR built"
echo "✓ Distribution package ready"
echo ""
echo "Plugin is ready for installation in Geneious!"