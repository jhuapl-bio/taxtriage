#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path
import json
from datetime import datetime

def run_command(cmd, cwd=None, capture_output=True):
    """Run a command and return success, stdout, stderr"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=capture_output,
            text=True,
            timeout=300  # 5 minute timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)

def count_files(directory, pattern="*.class"):
    """Count files matching pattern in directory"""
    try:
        path = Path(directory)
        if not path.exists():
            return 0
        if pattern == "*.class":
            return len(list(path.rglob("*.class")))
        elif pattern == "*.java":
            return len(list(path.rglob("*.java")))
        else:
            return len(list(path.rglob(pattern)))
    except Exception:
        return 0

def get_file_size(filepath):
    """Get file size in human readable format"""
    try:
        size = os.path.getsize(filepath)
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size < 1024.0:
                return f"{size:.1f} {unit}"
            size /= 1024.0
        return f"{size:.1f} TB"
    except Exception:
        return "Unknown"

def main():
    plugin_dir = "/Users/dho/Documents/taxtriage/geneious-plugin"
    test_report = {
        "timestamp": datetime.now().isoformat(),
        "plugin_directory": plugin_dir,
        "results": {},
        "errors": [],
        "summary": {},
        "files_created": []
    }

    print("=" * 60)
    print("TaxTriage Geneious Plugin Comprehensive Test Suite")
    print("=" * 60)
    print(f"Plugin directory: {plugin_dir}")
    print(f"Test started: {test_report['timestamp']}")
    print()

    if not os.path.exists(plugin_dir):
        print(f"❌ Plugin directory not found: {plugin_dir}")
        return 1

    os.chdir(plugin_dir)

    # Test 1: Source code analysis
    print("📊 PHASE 1: Source Code Analysis")
    print("-" * 40)

    src_java_count = count_files("src", "*.java")
    test_java_count = count_files("test", "*.java")

    print(f"Source Java files: {src_java_count}")
    print(f"Test Java files: {test_java_count}")

    test_report["summary"]["source_files"] = src_java_count
    test_report["summary"]["test_files"] = test_java_count

    # List main source files
    src_path = Path("src")
    if src_path.exists():
        java_files = list(src_path.rglob("*.java"))
        print("Main source files:")
        for java_file in java_files[:10]:  # Show first 10
            rel_path = java_file.relative_to(src_path)
            print(f"  ✓ {rel_path}")
        if len(java_files) > 10:
            print(f"  ... and {len(java_files) - 10} more")

    print()

    # Test 2: Dependencies check
    print("📋 PHASE 2: Dependencies Check")
    print("-" * 40)

    required_jars = [
        "lib/GeneiousPublicAPI.jar",
        "lib/jdom.jar",
        "lib/jebl.jar"
    ]

    dependencies_ok = True
    for jar in required_jars:
        if os.path.exists(jar):
            size = get_file_size(jar)
            print(f"  ✓ {jar} ({size})")
        else:
            print(f"  ❌ {jar} - Missing")
            dependencies_ok = False
            test_report["errors"].append(f"Missing dependency: {jar}")

    test_report["results"]["dependencies"] = dependencies_ok
    print()

    # Test 3: Clean build
    print("🧹 PHASE 3: Clean Build")
    print("-" * 40)

    success, stdout, stderr = run_command("ant clean")
    if success:
        print("  ✓ Clean successful")
        test_report["results"]["clean"] = True
    else:
        print(f"  ❌ Clean failed: {stderr}")
        test_report["results"]["clean"] = False
        test_report["errors"].append(f"Clean failed: {stderr}")

    print()

    # Test 4: Compilation
    print("🔨 PHASE 4: Source Compilation")
    print("-" * 40)

    success, stdout, stderr = run_command("ant compile")
    if success:
        print("  ✓ Compilation successful")
        test_report["results"]["compilation"] = True

        # Count compiled classes
        class_count = count_files("classes", "*.class")
        print(f"  📦 Compiled {class_count} class files")
        test_report["summary"]["compiled_classes"] = class_count

        # List some compiled classes
        classes_path = Path("classes")
        if classes_path.exists():
            class_files = list(classes_path.rglob("*.class"))
            print("  Key compiled classes:")
            for class_file in sorted(class_files)[:8]:
                rel_path = class_file.relative_to(classes_path)
                class_name = str(rel_path).replace("/", ".").replace(".class", "")
                print(f"    ✓ {class_name}")

            test_report["files_created"] = [str(f) for f in class_files]

    else:
        print(f"  ❌ Compilation failed")
        print(f"  Error: {stderr}")
        test_report["results"]["compilation"] = False
        test_report["errors"].append(f"Compilation failed: {stderr}")

    print()

    # Test 5: Test compilation (optional)
    print("🧪 PHASE 5: Test Compilation")
    print("-" * 40)

    success, stdout, stderr = run_command("ant compile-tests")
    if success:
        print("  ✓ Test compilation successful")
        test_report["results"]["test_compilation"] = True
    else:
        print("  ⚠️ Test compilation failed (may be expected)")
        print(f"  Info: {stderr}")
        test_report["results"]["test_compilation"] = False

    print()

    # Test 6: Basic plugin validation
    print("🔍 PHASE 6: Plugin Validation")
    print("-" * 40)

    # Check if main plugin class exists
    main_class_path = "classes/com/jhuapl/taxtriage/geneious/TaxTriagePlugin.class"
    if os.path.exists(main_class_path):
        print("  ✓ Main plugin class compiled")
        test_report["results"]["main_class"] = True
    else:
        print("  ❌ Main plugin class missing")
        test_report["results"]["main_class"] = False
        test_report["errors"].append("Main plugin class not found")

    # Check if operation class exists
    operation_class_path = "classes/com/jhuapl/taxtriage/geneious/TaxTriageOperation.class"
    if os.path.exists(operation_class_path):
        print("  ✓ Operation class compiled")
        test_report["results"]["operation_class"] = True
    else:
        print("  ❌ Operation class missing")
        test_report["results"]["operation_class"] = False
        test_report["errors"].append("Operation class not found")

    # Check if options class exists
    options_class_path = "classes/com/jhuapl/taxtriage/geneious/TaxTriageOptions.class"
    if os.path.exists(options_class_path):
        print("  ✓ Options class compiled")
        test_report["results"]["options_class"] = True
    else:
        print("  ❌ Options class missing")
        test_report["results"]["options_class"] = False
        test_report["errors"].append("Options class not found")

    print()

    # Test 7: Plugin build
    print("📦 PHASE 7: Plugin Build")
    print("-" * 40)

    if test_report["results"].get("compilation", False):
        success, stdout, stderr = run_command("ant build")
        if success:
            print("  ✓ Plugin build successful")
            test_report["results"]["build"] = True

            # Check for JAR file
            jar_path = "build/TaxTriage/TaxTriage.jar"
            if os.path.exists(jar_path):
                size = get_file_size(jar_path)
                print(f"  📦 Created JAR: {jar_path} ({size})")
                test_report["summary"]["jar_size"] = size
            else:
                print("  ⚠️ JAR file not found at expected location")

        else:
            print(f"  ❌ Plugin build failed: {stderr}")
            test_report["results"]["build"] = False
            test_report["errors"].append(f"Build failed: {stderr}")
    else:
        print("  ⏭️ Skipping build (compilation failed)")
        test_report["results"]["build"] = False

    print()

    # Test 8: Distribution package
    print("📤 PHASE 8: Distribution Package")
    print("-" * 40)

    if test_report["results"].get("build", False):
        success, stdout, stderr = run_command("ant distribute")
        if success:
            print("  ✓ Distribution package created")
            test_report["results"]["distribution"] = True

            # Check for gplugin file
            gplugin_path = "build/TaxTriage.gplugin"
            if os.path.exists(gplugin_path):
                size = get_file_size(gplugin_path)
                print(f"  📦 Created plugin: {gplugin_path} ({size})")
                test_report["summary"]["plugin_size"] = size
                test_report["summary"]["plugin_ready"] = True
            else:
                print("  ⚠️ Plugin file not found")
                test_report["summary"]["plugin_ready"] = False

        else:
            print(f"  ❌ Distribution failed: {stderr}")
            test_report["results"]["distribution"] = False
            test_report["errors"].append(f"Distribution failed: {stderr}")
    else:
        print("  ⏭️ Skipping distribution (build failed)")
        test_report["results"]["distribution"] = False

    print()

    # Test 9: Manual test execution (if possible)
    print("🧪 PHASE 9: Manual Test Execution")
    print("-" * 40)

    if test_report["results"].get("compilation", False):
        # Try to run the manual test class
        test_class = "com.jhuapl.taxtriage.geneious.TaxTriagePluginTest"
        classpath = "classes:lib/*"

        print(f"  Attempting to run: {test_class}")
        success, stdout, stderr = run_command(f"java -cp '{classpath}' {test_class}")

        if success:
            print("  ✓ Manual test execution successful")
            print("  Test output:")
            for line in stdout.split('\n'):
                if line.strip():
                    print(f"    {line}")
            test_report["results"]["manual_test"] = True
        else:
            print("  ⚠️ Manual test execution failed (expected without full Geneious)")
            print(f"  Info: {stderr}")
            test_report["results"]["manual_test"] = False
    else:
        print("  ⏭️ Skipping manual test (compilation failed)")
        test_report["results"]["manual_test"] = False

    print()

    # Generate summary
    print("📋 FINAL SUMMARY")
    print("=" * 60)

    passed_tests = sum(1 for result in test_report["results"].values() if result)
    total_tests = len(test_report["results"])

    print(f"Test Results: {passed_tests}/{total_tests} passed")
    print()

    print("✅ PASSED:")
    for test_name, result in test_report["results"].items():
        if result:
            print(f"  ✓ {test_name.replace('_', ' ').title()}")

    print()
    failed_tests = [name for name, result in test_report["results"].items() if not result]
    if failed_tests:
        print("❌ FAILED:")
        for test_name in failed_tests:
            print(f"  ❌ {test_name.replace('_', ' ').title()}")
        print()

    if test_report["errors"]:
        print("⚠️ ERRORS:")
        for error in test_report["errors"]:
            print(f"  - {error}")
        print()

    # Key metrics
    print("📊 KEY METRICS:")
    print(f"  Source files: {test_report['summary'].get('source_files', 0)}")
    print(f"  Test files: {test_report['summary'].get('test_files', 0)}")
    print(f"  Compiled classes: {test_report['summary'].get('compiled_classes', 0)}")

    if 'jar_size' in test_report['summary']:
        print(f"  JAR size: {test_report['summary']['jar_size']}")
    if 'plugin_size' in test_report['summary']:
        print(f"  Plugin size: {test_report['summary']['plugin_size']}")

    print()

    # Overall status
    critical_tests = ['compilation', 'build', 'distribution']
    critical_passed = sum(1 for test in critical_tests if test_report["results"].get(test, False))

    if critical_passed == len(critical_tests):
        print("🎉 OVERALL STATUS: SUCCESS")
        print("   Plugin is ready for installation in Geneious!")
        overall_status = 0
    elif critical_passed >= 2:
        print("⚠️ OVERALL STATUS: PARTIAL SUCCESS")
        print("   Plugin compiled but has issues with packaging")
        overall_status = 1
    else:
        print("❌ OVERALL STATUS: FAILURE")
        print("   Plugin has critical compilation or build issues")
        overall_status = 2

    print()
    print("📁 NEXT STEPS:")
    if test_report["summary"].get("plugin_ready", False):
        print("  1. Install the plugin: build/TaxTriage.gplugin")
        print("  2. Copy to Geneious plugins directory")
        print("  3. Restart Geneious to load the plugin")
    else:
        print("  1. Fix compilation errors")
        print("  2. Ensure all dependencies are available")
        print("  3. Re-run the build process")

    # Save detailed report
    report_file = "/Users/dho/Documents/taxtriage/detailed_test_report.json"
    try:
        with open(report_file, 'w') as f:
            json.dump(test_report, f, indent=2)
        print(f"\n📄 Detailed report saved: {report_file}")
    except Exception as e:
        print(f"\n⚠️ Could not save detailed report: {e}")

    print()
    print("=" * 60)

    return overall_status

if __name__ == "__main__":
    sys.exit(main())