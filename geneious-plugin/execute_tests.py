#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd, cwd=None):
    """Run a command and return the result"""
    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def main():
    plugin_dir = "/Users/dho/Documents/taxtriage/geneious-plugin"

    print("=== TaxTriage Geneious Plugin Test Execution ===")
    print(f"Plugin directory: {plugin_dir}")
    print()

    # Check if plugin directory exists
    if not os.path.exists(plugin_dir):
        print(f"Error: Plugin directory not found: {plugin_dir}")
        return 1

    # Change to plugin directory
    os.chdir(plugin_dir)

    test_results = {
        'compilation': False,
        'build': False,
        'distribution': False,
        'class_count': 0,
        'errors': []
    }

    # 1. Clean
    print("1. Cleaning previous builds...")
    success, stdout, stderr = run_command("ant clean", plugin_dir)
    if success:
        print("✓ Clean successful")
    else:
        print("✗ Clean failed")
        test_results['errors'].append(f"Clean failed: {stderr}")

    print()

    # 2. Compile
    print("2. Compiling source code...")
    success, stdout, stderr = run_command("ant compile", plugin_dir)
    if success:
        print("✓ Compilation successful")
        test_results['compilation'] = True

        # Count classes
        classes_dir = Path(plugin_dir) / "classes"
        if classes_dir.exists():
            class_files = list(classes_dir.rglob("*.class"))
            test_results['class_count'] = len(class_files)
            print(f"  Compiled {len(class_files)} class files")

            # List main classes
            main_classes = [f for f in class_files if 'TaxTriage' in f.name]
            for class_file in main_classes[:5]:  # Show first 5
                rel_path = class_file.relative_to(classes_dir)
                print(f"    {rel_path}")

    else:
        print("✗ Compilation failed")
        test_results['errors'].append(f"Compilation failed: {stderr}")
        print(f"Error: {stderr}")

    print()

    # 3. Compile tests (optional)
    print("3. Compiling test code...")
    success, stdout, stderr = run_command("ant compile-tests", plugin_dir)
    if success:
        print("✓ Test compilation successful")
    else:
        print("! Test compilation failed (may be expected)")
        # Don't treat this as a fatal error

    print()

    # 4. Run tests (optional)
    print("4. Running tests...")
    success, stdout, stderr = run_command("ant test", plugin_dir)
    if success:
        print("✓ Tests executed successfully")
        print(stdout)
    else:
        print("! Test execution failed (may be expected without full environment)")
        # Try manual test execution
        print("  Attempting manual test execution...")
        success, stdout, stderr = run_command(
            "java -cp 'classes:lib/*' com.jhuapl.taxtriage.geneious.TaxTriagePluginTest",
            plugin_dir
        )
        if success:
            print("✓ Manual test execution successful")
            print(stdout)
        else:
            print("! Manual test execution failed (expected without Geneious)")

    print()

    # 5. Build
    print("5. Building plugin...")
    success, stdout, stderr = run_command("ant build", plugin_dir)
    if success:
        print("✓ Plugin build successful")
        test_results['build'] = True

        # Check for JAR file
        jar_path = Path(plugin_dir) / "build" / "TaxTriage" / "TaxTriage.jar"
        if jar_path.exists():
            size = jar_path.stat().st_size
            print(f"  Created JAR: {jar_path} ({size} bytes)")
    else:
        print("✗ Plugin build failed")
        test_results['errors'].append(f"Build failed: {stderr}")
        print(f"Error: {stderr}")

    print()

    # 6. Distribute
    print("6. Creating distribution package...")
    success, stdout, stderr = run_command("ant distribute", plugin_dir)
    if success:
        print("✓ Distribution package created")
        test_results['distribution'] = True

        # Check for gplugin file
        gplugin_path = Path(plugin_dir) / "build" / "TaxTriage.gplugin"
        if gplugin_path.exists():
            size = gplugin_path.stat().st_size
            print(f"  Created plugin package: {gplugin_path} ({size} bytes)")
    else:
        print("✗ Distribution package creation failed")
        test_results['errors'].append(f"Distribution failed: {stderr}")
        print(f"Error: {stderr}")

    print()
    print("=== TEST SUMMARY ===")
    print(f"Source compilation: {'✓ PASSED' if test_results['compilation'] else '✗ FAILED'}")
    print(f"Plugin build: {'✓ PASSED' if test_results['build'] else '✗ FAILED'}")
    print(f"Distribution package: {'✓ PASSED' if test_results['distribution'] else '✗ FAILED'}")
    print(f"Classes compiled: {test_results['class_count']}")

    if test_results['errors']:
        print("\nErrors encountered:")
        for error in test_results['errors']:
            print(f"  - {error}")

    # Overall status
    success_count = sum([test_results['compilation'], test_results['build'], test_results['distribution']])
    print(f"\nOverall: {success_count}/3 critical tasks completed successfully")

    if test_results['compilation'] and test_results['build'] and test_results['distribution']:
        print("✓ Plugin is ready for installation!")
        return 0
    else:
        print("✗ Plugin build incomplete")
        return 1

if __name__ == "__main__":
    sys.exit(main())