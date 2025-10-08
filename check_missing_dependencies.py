#!/usr/bin/env python3

import os
from pathlib import Path

def check_plugin_structure():
    """Check the plugin structure and identify missing dependencies"""
    plugin_dir = "/Users/dho/Documents/taxtriage/geneious-plugin"

    print("TaxTriage Plugin Dependency Check")
    print("=" * 50)
    print(f"Plugin directory: {plugin_dir}")
    print()

    if not os.path.exists(plugin_dir):
        print("❌ Plugin directory not found!")
        return False

    # Check essential files
    essential_files = {
        "build.xml": "Ant build script",
        "plugin.properties": "Plugin configuration",
        "src/com/jhuapl/taxtriage/geneious/TaxTriagePlugin.java": "Main plugin class",
        "src/com/jhuapl/taxtriage/geneious/TaxTriageOperation.java": "Plugin operation",
        "src/com/jhuapl/taxtriage/geneious/TaxTriageOptions.java": "Plugin options"
    }

    print("📋 Essential Files Check:")
    all_essential_present = True
    for file_path, description in essential_files.items():
        full_path = os.path.join(plugin_dir, file_path)
        if os.path.exists(full_path):
            print(f"  ✓ {file_path} - {description}")
        else:
            print(f"  ❌ {file_path} - {description} (MISSING)")
            all_essential_present = False

    print()

    # Check dependencies
    lib_dir = os.path.join(plugin_dir, "lib")
    print("📚 Dependencies Check:")

    if not os.path.exists(lib_dir):
        print("  ❌ lib directory not found!")
        return False

    required_jars = [
        "GeneiousPublicAPI.jar",
        "jdom.jar",
        "jebl.jar"
    ]

    optional_jars = [
        "junit-jupiter-api-5.8.2.jar",
        "junit-jupiter-engine-5.8.2.jar",
        "junit-platform-launcher-1.8.2.jar",
        "mockito-core-4.6.1.jar",
        "mockito-junit-jupiter-4.6.1.jar",
        "byte-buddy-1.12.10.jar",
        "objenesis-3.2.jar"
    ]

    deps_ok = True
    for jar in required_jars:
        jar_path = os.path.join(lib_dir, jar)
        if os.path.exists(jar_path):
            size = os.path.getsize(jar_path)
            print(f"  ✓ {jar} ({size:,} bytes)")
        else:
            print(f"  ❌ {jar} (REQUIRED - MISSING)")
            deps_ok = False

    print("\n  Optional (for testing):")
    for jar in optional_jars:
        jar_path = os.path.join(lib_dir, jar)
        if os.path.exists(jar_path):
            size = os.path.getsize(jar_path)
            print(f"  ✓ {jar} ({size:,} bytes)")
        else:
            print(f"  - {jar} (optional)")

    print()

    # List all JAR files actually present
    print("📦 JAR Files Present:")
    jar_files = list(Path(lib_dir).glob("*.jar"))
    if jar_files:
        for jar_file in sorted(jar_files):
            size = jar_file.stat().st_size
            print(f"  • {jar_file.name} ({size:,} bytes)")
    else:
        print("  No JAR files found in lib directory")

    print()

    # Check source structure
    src_dir = os.path.join(plugin_dir, "src")
    print("📁 Source Structure:")

    if os.path.exists(src_dir):
        java_files = list(Path(src_dir).rglob("*.java"))
        print(f"  Total Java files: {len(java_files)}")

        # Group by package
        packages = {}
        for java_file in java_files:
            rel_path = java_file.relative_to(Path(src_dir))
            package = str(rel_path.parent).replace("/", ".")
            if package not in packages:
                packages[package] = []
            packages[package].append(java_file.name)

        for package, files in sorted(packages.items()):
            print(f"  📂 {package}:")
            for file in sorted(files):
                print(f"    • {file}")
    else:
        print("  ❌ src directory not found!")
        all_essential_present = False

    print()

    # Missing dependencies check
    print("🔍 Missing Dependencies Analysis:")

    # Check for classes that might need creation
    missing_classes = []

    # Look for imports in existing source files that might reference missing classes
    if os.path.exists(src_dir):
        java_files = list(Path(src_dir).rglob("*.java"))

        for java_file in java_files:
            try:
                with open(java_file, 'r') as f:
                    content = f.read()

                # Look for imports that might be missing
                if 'com.jhuapl.taxtriage.geneious.config' in content:
                    config_classes = ['ConfigGenerator', 'SampleSheetBuilder']
                    for cls in config_classes:
                        if cls in content:
                            config_path = src_dir + f"/com/jhuapl/taxtriage/geneious/config/{cls}.java"
                            if not os.path.exists(config_path):
                                missing_classes.append(f"config.{cls}")

                if 'com.jhuapl.taxtriage.geneious.docker' in content:
                    docker_classes = ['DockerManager', 'DockerException', 'ExecutionResult']
                    for cls in docker_classes:
                        if cls in content:
                            docker_path = src_dir + f"/com/jhuapl/taxtriage/geneious/docker/{cls}.java"
                            if not os.path.exists(docker_path):
                                missing_classes.append(f"docker.{cls}")
            except Exception as e:
                print(f"  Warning: Could not read {java_file}: {e}")

    if missing_classes:
        print("  ❌ Potentially missing classes:")
        for cls in sorted(set(missing_classes)):
            print(f"    - {cls}")
    else:
        print("  ✓ No obviously missing classes detected")

    print()

    # Overall assessment
    print("📊 Overall Assessment:")
    if all_essential_present and deps_ok:
        print("  ✅ Plugin structure looks complete")
        print("  ✅ All required dependencies present")
        print("  🎯 Ready for compilation")
        return True
    else:
        print("  ❌ Plugin has missing components")
        if not all_essential_present:
            print("  ❌ Essential files are missing")
        if not deps_ok:
            print("  ❌ Required dependencies are missing")
        print("  🔧 Requires fixes before compilation")
        return False

if __name__ == "__main__":
    success = check_plugin_structure()
    exit(0 if success else 1)