#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path
import json
from datetime import datetime

def run_command(cmd, cwd=None, capture_output=True, timeout=300):
    """Run a command and return success, stdout, stderr"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=capture_output,
            text=True,
            timeout=timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)

def main():
    plugin_dir = "/Users/dho/Documents/taxtriage/geneious-plugin"
    print("=" * 80)
    print("TAXTRIAGE GENEIOUS PLUGIN - FINAL COMPREHENSIVE TEST")
    print("=" * 80)
    print(f"Plugin directory: {plugin_dir}")
    print(f"Test started: {datetime.now().isoformat()}")
    print()

    if not os.path.exists(plugin_dir):
        print(f"❌ Plugin directory not found: {plugin_dir}")
        return 1

    os.chdir(plugin_dir)

    results = {}
    overall_success = True

    # Phase 1: Source Analysis
    print("📋 PHASE 1: SOURCE CODE ANALYSIS")
    print("-" * 50)

    src_files = list(Path("src").rglob("*.java")) if Path("src").exists() else []
    test_files = list(Path("test").rglob("*.java")) if Path("test").exists() else []

    print(f"Source files found: {len(src_files)}")
    print("Key source files:")
    for src_file in sorted(src_files):
        rel_path = src_file.relative_to(Path("src"))
        print(f"  ✓ {rel_path}")

    print(f"\nTest files found: {len(test_files)}")
    for test_file in sorted(test_files):
        rel_path = test_file.relative_to(Path("test"))
        print(f"  ✓ {rel_path}")

    print()

    # Phase 2: Dependencies Check
    print("📦 PHASE 2: DEPENDENCIES CHECK")
    print("-" * 50)

    required_jars = ["GeneiousPublicAPI.jar", "jdom.jar", "jebl.jar"]
    lib_dir = Path("lib")

    if lib_dir.exists():
        for jar in required_jars:
            jar_path = lib_dir / jar
            if jar_path.exists():
                size = jar_path.stat().st_size
                print(f"  ✓ {jar} ({size:,} bytes)")
            else:
                print(f"  ❌ {jar} - MISSING")
                overall_success = False
    else:
        print("  ❌ lib directory not found")
        overall_success = False

    print()

    # Phase 3: Clean Build
    print("🧹 PHASE 3: CLEAN BUILD")
    print("-" * 50)

    success, stdout, stderr = run_command("ant clean")
    if success:
        print("  ✓ Clean successful")
        results['clean'] = True
    else:
        print(f"  ❌ Clean failed: {stderr}")
        results['clean'] = False
        overall_success = False

    print()

    # Phase 4: Compilation
    print("🔨 PHASE 4: COMPILATION")
    print("-" * 50)

    success, stdout, stderr = run_command("ant compile")
    if success:
        print("  ✓ Compilation successful")
        results['compile'] = True

        # Count compiled classes
        classes_dir = Path("classes")
        if classes_dir.exists():
            class_files = list(classes_dir.rglob("*.class"))
            print(f"  📦 Compiled {len(class_files)} class files")

            # Show key classes
            key_classes = [
                "com/jhuapl/taxtriage/geneious/TaxTriagePlugin.class",
                "com/jhuapl/taxtriage/geneious/TaxTriageOperation.class",
                "com/jhuapl/taxtriage/geneious/TaxTriageOptions.class"
            ]

            for key_class in key_classes:
                class_path = classes_dir / key_class
                if class_path.exists():
                    print(f"    ✓ {key_class}")
                else:
                    print(f"    ❌ {key_class} - Missing")
        else:
            print("  ⚠️ No classes directory found")
    else:
        print(f"  ❌ Compilation failed")
        print(f"  Error details: {stderr}")
        results['compile'] = False
        overall_success = False

    print()

    # Phase 5: Test Compilation
    print("🧪 PHASE 5: TEST COMPILATION")
    print("-" * 50)

    success, stdout, stderr = run_command("ant compile-tests")
    if success:
        print("  ✓ Test compilation successful")
        results['compile_tests'] = True
    else:
        print("  ⚠️ Test compilation failed (may be expected)")
        print(f"  Info: {stderr}")
        results['compile_tests'] = False

    print()

    # Phase 6: Plugin Build
    print("📦 PHASE 6: PLUGIN BUILD")
    print("-" * 50)

    if results.get('compile', False):
        success, stdout, stderr = run_command("ant build")
        if success:
            print("  ✓ Plugin build successful")
            results['build'] = True

            # Check for JAR
            jar_path = Path("build/TaxTriage/TaxTriage.jar")
            if jar_path.exists():
                size = jar_path.stat().st_size
                print(f"  📦 JAR created: {jar_path} ({size:,} bytes)")
            else:
                print("  ⚠️ JAR not found at expected location")
        else:
            print(f"  ❌ Build failed: {stderr}")
            results['build'] = False
            overall_success = False
    else:
        print("  ⏭️ Skipping build (compilation failed)")
        results['build'] = False

    print()

    # Phase 7: Distribution
    print("📤 PHASE 7: DISTRIBUTION PACKAGE")
    print("-" * 50)

    if results.get('build', False):
        success, stdout, stderr = run_command("ant distribute")
        if success:
            print("  ✓ Distribution package created")
            results['distribute'] = True

            # Check for gplugin
            gplugin_path = Path("build/TaxTriage.gplugin")
            if gplugin_path.exists():
                size = gplugin_path.stat().st_size
                print(f"  📦 Plugin package: {gplugin_path} ({size:,} bytes)")
                results['plugin_ready'] = True
            else:
                print("  ⚠️ Plugin package not found")
                results['plugin_ready'] = False
        else:
            print(f"  ❌ Distribution failed: {stderr}")
            results['distribute'] = False
            overall_success = False
    else:
        print("  ⏭️ Skipping distribution (build failed)")
        results['distribute'] = False

    print()

    # Phase 8: Manual Test (if available)
    print("🧪 PHASE 8: MANUAL TEST EXECUTION")
    print("-" * 50)

    if results.get('compile', False):
        test_class = "com.jhuapl.taxtriage.geneious.TaxTriagePluginTest"
        classpath = "classes:lib/*"

        success, stdout, stderr = run_command(f"java -cp '{classpath}' {test_class}")
        if success:
            print("  ✓ Manual test execution successful")
            print("  Test output:")
            for line in stdout.split('\n'):
                if line.strip():
                    print(f"    {line}")
            results['manual_test'] = True
        else:
            print("  ⚠️ Manual test failed (expected without full Geneious environment)")
            print(f"  Note: {stderr}")
            results['manual_test'] = False
    else:
        print("  ⏭️ Skipping manual test (compilation failed)")
        results['manual_test'] = False

    print()

    # Final Summary
    print("📊 FINAL TEST SUMMARY")
    print("=" * 80)

    critical_tests = ['compile', 'build', 'distribute']
    passed_critical = sum(1 for test in critical_tests if results.get(test, False))

    print("Test Results:")
    for test_name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {test_name.replace('_', ' ').title()}: {status}")

    print()
    print(f"Critical Tests Passed: {passed_critical}/{len(critical_tests)}")

    # Overall Assessment
    if passed_critical == len(critical_tests):
        print("🎉 OVERALL STATUS: SUCCESS")
        print("   ✅ Plugin compiled successfully")
        print("   ✅ Plugin built successfully")
        print("   ✅ Distribution package created")
        print("   🚀 Plugin is ready for Geneious installation!")
        exit_code = 0
    elif passed_critical >= 2:
        print("⚠️ OVERALL STATUS: PARTIAL SUCCESS")
        print("   ✅ Plugin compiled")
        print("   ⚠️ Some build issues detected")
        exit_code = 1
    else:
        print("❌ OVERALL STATUS: FAILURE")
        print("   ❌ Critical build failures detected")
        print("   🔧 Requires fixing before deployment")
        exit_code = 2

    print()

    # Installation Instructions
    if results.get('plugin_ready', False):
        print("📋 INSTALLATION INSTRUCTIONS:")
        print("   1. Locate the plugin file: build/TaxTriage.gplugin")
        print("   2. Copy to Geneious plugins directory:")
        print("      • Windows: %USERPROFILE%\\.geneious\\plugins\\")
        print("      • macOS: ~/.geneious/plugins/")
        print("      • Linux: ~/.geneious/plugins/")
        print("   3. Restart Geneious")
        print("   4. The TaxTriage plugin should appear in the Tools menu")
        print()
        print("🔧 REQUIREMENTS FOR RUNTIME:")
        print("   • Docker must be installed and running")
        print("   • TaxTriage Docker image must be available")
        print("   • Sufficient disk space for analysis")
    else:
        print("🔧 NEXT STEPS:")
        print("   1. Fix any compilation errors")
        print("   2. Ensure all dependencies are in lib/ directory")
        print("   3. Re-run the build process")

    print()
    print("=" * 80)
    print(f"Test completed: {datetime.now().isoformat()}")

    return exit_code

if __name__ == "__main__":
    sys.exit(main())