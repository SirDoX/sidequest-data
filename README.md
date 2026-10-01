# sidequest-data

Prepares the walkable street network for the Sidequest app, so the app
doesn't have to ask the public Overpass servers tile by tile.

Every Sunday a GitHub Actions job (`.github/workflows/build.yml`)

1. downloads the DACH extract from Geofabrik,
2. keeps only highways (`osmium tags-filter`),
3. writes every walkable way into zoom-13 tiles (~3 × 3 km) with
   `pipeline/build_packs.py`, using the same filters as the app,
4. publishes the result as a release: `manifest.json`, an index and the
   data parts. The app reads the latest release and fetches single tiles
   with HTTP range requests.

Run it locally, e.g. for Berlin:

```bash
pip install -r pipeline/requirements.txt
curl -LO https://download.geofabrik.de/europe/germany/berlin-latest.osm.pbf
curl -LO https://download.geofabrik.de/europe/germany/berlin.poly
python pipeline/build_packs.py berlin-latest.osm.pbf berlin.poly out berlin test
```

Data © OpenStreetMap contributors, available under the Open Database
License (ODbL).
