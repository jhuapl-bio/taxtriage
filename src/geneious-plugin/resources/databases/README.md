# Bundling Databases with TaxTriage Plugin

## Overview

The TaxTriage plugin supports bundling reference databases to avoid repeated downloads. This directory is where bundled databases should be placed before building the plugin.

## Database Structure

Each database should be placed in its own subdirectory:

```
resources/databases/
├── viral/          # Viral database files
│   ├── hash.k2d    # Kraken2 hash file
│   ├── opts.k2d    # Kraken2 options file
│   ├── taxo.k2d    # Kraken2 taxonomy file
│   └── seqid2taxid.map  # Sequence ID mapping
├── standard/       # Standard database (optional)
└── minikraken/     # MiniKraken database (optional)
```

## Adding Viral Database

1. **Download the viral database:**

   ```bash
   # Example using Kraken2
   kraken2-build --download-taxonomy --db viral_db
   kraken2-build --download-library viral --db viral_db
   kraken2-build --build --db viral_db --threads 8
   ```

2. **Copy database files to this directory:**

   ```bash
   cp viral_db/*.k2d resources/databases/viral/
   cp viral_db/*.map resources/databases/viral/
   ```

3. **Rebuild the plugin:**
   ```bash
   ant clean distribute
   ```

## Database Size Considerations

- **Viral Database**: ~200-500MB (recommended for bundling)
- **MiniKraken**: ~5-8GB (consider caching instead)
- **Standard**: ~50-100GB (too large for bundling)

## Database Versioning

The DatabaseManager automatically tracks database versions. When bundling:

- Viral database version: Set in DatabaseType.VIRAL.getBundledVersion()
- Update version string when updating bundled databases

## Notes

- Only include essential databases to keep plugin size manageable
- Larger databases will be downloaded and cached on first use
- The plugin will prefer cached databases over bundled ones if newer
