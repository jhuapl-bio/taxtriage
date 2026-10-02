# Offline report libraries

Local copies of everything the interactive report (`all.odr.html`) otherwise
loads from the internet when it is opened. Used to build a report that opens
with **no network access**:

| File(s)                                                                        | Library                                                      | License                      |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------ | ---------------------------- |
| `d3.min.js`                                                                    | D3 7.8.5                                                     | ISC                          |
| `xlsx.full.min.js`                                                             | SheetJS Community Edition 0.18.5                             | Apache-2.0                   |
| `jspdf.umd.min.js`                                                             | jsPDF 2.5.1                                                  | MIT                          |
| `leaflet.js`, `leaflet.css`, `layers*.png`, `marker-icon.png`                  | Leaflet 1.9.4                                                | BSD-2-Clause                 |
| `leaflet.markercluster.js`, `MarkerCluster.css`                                | Leaflet.markercluster 1.5.3                                  | MIT                          |
| `all.min.css`, `fa-*.woff2`                                                    | Font Awesome Free 6.5.0                                      | CSS: MIT, fonts: SIL OFL 1.1 |
| `ne_110m_admin_0_countries.geojson`, `ne_50m_admin_1_states_provinces.geojson` | Natural Earth boundaries (slimmed: names + rounded geometry) | Public domain                |

`manifest.json` maps each CDN URL in `assets/heatmap.html` to its file here, so
the build embeds exactly the version the template asks for and warns when this
folder is out of date.

How it is used:

- `--offline_report` embeds these files and downloads only what is missing.
- `--offline_report_files <dir>` embeds a folder like this one with no network
  at all (pass `assets/offline_report_libs` to use this copy).

Refresh after bumping a library version in `assets/heatmap.html` (needs internet):

```bash
python scripts/fetch_offline_report_libs.py
```
