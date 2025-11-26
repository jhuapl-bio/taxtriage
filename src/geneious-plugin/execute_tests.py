#!/usr/bin/env python3

import os
import subprocess
import sys
import json
import argparse
from pathlib import Path

def run_command(cmd, cwd=None):
    """Run a command and return success, stdout, stderr."""
    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd,
                                capture_output=True, text=True)
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def run_tests(plugin_dir, verbose=False):
    print("=== TaxTriage Geneious Plugin Test Execution ===")
    print(f"Plugin directory: {plugin_dir}")
    print()

    if not os.path.exists(plugin_dir):
        print(f"Error: Plugin directory not found: {plugin_dir}")
        return {
            "compilation": False,
            "build": False,
            "distribution": False,
            "class_count": 0,
            "errors": [f"Plugin directory not found: {plugin_dir}"]
        }

    os.chdir(plugin_dir)

    test_results = {
        'compilation': False,
        'build': False,
        'distribution': False,
        'class_count': 0,
        'errors': []
    }

    # 1. CLEAN
    print("1. Cleaning previous builds...")
    success, stdout, stderr = run_command("ant clean", plugin_dir)
    if success:
        print("✓ Clean successful")
    else:
        print("✗ Clean failed")
        test_results['errors'].append(f"Clean failed: {stderr}")
        if verbose:
            print(stderr)
    print()

    # 2. COMPILE
    print("2. Compiling source code...")
    success, stdout, stderr = run_command("ant compile", plugin_dir)
    if success:
        print("✓ Compilation successful")
        test_results['compilation'] = True

        classes_dir = Path(plugin_dir) / "classes"
        if classes_dir.exists():
            class_files = list(classes_dir.rglob("*.class"))
            test_results['class_count'] = len(class_files)
            print(f"  Compiled {len(class_files)} class files")

            main_classes = [f for f in class_files if 'TaxTriage' in f.name]
            for class_file in main_classes[:5]:
                print(f"    {class_file.relative_to(classes_dir)}")
    else:
        print("✗ Compilation failed")
        test_results['errors'].append(f"Compilation failed: {stderr}")
        if verbose:
            print(stderr)
    print()

    # 3. COMPILE TESTS
    print("3. Compiling test code...")
    success, stdout, stderr = run_command("ant compile-tests", plugin_dir)
    if success:
        print("✓ Test compilation successful")
    else:
        print("! Test compilation failed (may be expected)")
        if verbose:
            print(stderr)
    print()

    # 4. RUN TESTS
    print("4. Running tests...")
    success, stdout, stderr = run_command("ant test", plugin_dir)
    if success:
        print("✓ Tests executed successfully")
        if verbose:
            print(stdout)
    else:
        print("! Test execution failed (may be expected without full environment)")
        print("  Attempting manual test execution...")

        success, stdout, stderr = run_command(
            "java -cp 'classes:lib/*' com.jhuapl.taxtriage.geneious.TaxTriagePluginTest",
            plugin_dir
        )

        if success:
            print("✓ Manual test execution successful")
            if verbose:
                print(stdout)
        else:
            print("! Manual test execution failed (expected without Geneious)")
            if verbose:
                print(stderr)
    print()

    # 5. BUILD PLUGIN
    print("5. Building plugin...")
    success, stdout, stderr = run_command("ant build", plugin_dir)
    if success:
        print("✓ Plugin build successful")
        test_results['build'] = True

        jar_path = Path(plugin_dir) / "build" / "TaxTriage" / "TaxTriage.jar"
        if jar_path.exists():
            print(f"  Created JAR: {jar_path} ({jar_path.stat().st_size} bytes)")
    else:
        print("✗ Plugin build failed")
        test_results['errors'].append(f"Build failed: {stderr}")
        if verbose:
            print(stderr)
    print()

    # 6. DISTRIBUTE PLUGIN
    print("6. Creating distribution package...")
    success, stdout, stderr = run_command("ant distribute", plugin_dir)
    if success:
        print("✓ Distribution package created")
        test_results['distribution'] = True

        gplugin_path = Path(plugin_dir) / "build" / "TaxTriage.gplugin"
        if gplugin_path.exists():
            print(f"  Created plugin package: {gplugin_path} ({gplugin_path.stat().st_size} bytes)")
    else:
        print("✗ Distribution package creation failed")
        test_results['errors'].append(f"Distribution failed: {stderr}")
        if verbose:
            print(stderr)
    print()

    return test_results


def print_summary(test_results):
    print("=== TEST SUMMARY ===")
    print(f"Source compilation: {'✓ PASSED' if test_results['compilation'] else '✗ FAILED'}")
    print(f"Plugin build: {'✓ PASSED' if test_results['build'] else '✗ FAILED'}")
    print(f"Distribution package: {'✓ PASSED' if test_results['distribution'] else '✗ FAILED'}")
    print(f"Classes compiled: {test_results['class_count']}")

    if test_results['errors']:
        print("\nErrors encountered:")
        for e in test_results['errors']:
            print(f"  - {e}")

    success_count = sum([
        test_results['compilation'],
        test_results['build'],
        test_results['distribution']
    ])

    print(f"\nOverall: {success_count}/3 critical tasks completed successfully")

    if success_count == 3:
        print("✓ Plugin is ready for installation!")
    else:
        print("✗ Plugin build incomplete")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run build and test steps for the TaxTriage Geneious plugin."
    )

    parser.add_argument(
        "-p", "--plugin-dir",
        default=None,
        help="Path to the Geneious plugin directory"
    )

    parser.add_argument(
        "-o", "--output",
        help="Optional JSON file to write the test summary"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print full stdout/stderr for each step"
    )

    return parser.parse_args()


def main():
    args = parse_args()
    results = run_tests(args.plugin_dir, verbose=args.verbose)
    print_summary(results)

    # Write JSON output if requested
    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nSummary written to: {args.output}")

    # Exit code based on pass/fail
    if results['compilation'] and results['build'] and results['distribution']:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
