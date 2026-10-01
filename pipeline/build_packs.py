"""Builds Sidequest street-network packs from an OpenStreetMap extract.

The app used to load the walkable street network tile by tile from the
public Overpass API, which is slow and often overloaded. This script does
the same filtering once, offline, and writes the result as one archive the
app reads with HTTP range requests:

* ``<name>-<n>.bin``   concatenated gzip chunks, one per zoom-13 tile
                       (~3 x 3 km); split into parts below 1.9 GB, the
                       GitHub release asset limit being 2 GB
* ``<name>-index.bin.gz``  where each chunk is (see INDEX_ENTRY)
* ``manifest.json``    format, zoom, parts, creation date

A chunk is gzip-compressed JSON ``{"w": [[id, highway, name, [node ids],
polyline6], ...]}`` with every way that touches the tile; the app splits
them into segments exactly as it does with Overpass data. The filters
below must stay in sync with the app's OverpassClient.

Usage: build_packs.py EXTRACT.osm.pbf BOUNDARY.poly OUT_DIR NAME TAG
"""

import gzip
import io
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

import osmium
from shapely.geometry import Polygon, box
from shapely.prepared import prep

ZOOM = 13
PART_LIMIT = 1_900_000_000
FORMAT = 2

# x, y, part, flags, offset, length, street metres. Flags: 1 = the tile
# crosses the boundary of the extract (data beyond it is missing), 2 = no
# ways. Street metres: the streets that count for exploring (no sidewalks,
# crossings or cycleway connectors), each piece between two nodes counted
# in the tile of its middle, so the tiles of an area add up to its total.
INDEX_ENTRY = struct.Struct("<HHBBIII")
CONNECTORS = {"sidewalk", "crossing", "cycleway_connector"}
FLAG_PARTIAL = 1
FLAG_EMPTY = 2

STREET_HIGHWAYS = {
    "footway", "path", "pedestrian", "living_street", "residential",
    "unclassified", "tertiary", "tertiary_link", "secondary",
    "secondary_link", "primary", "primary_link", "service", "track", "steps",
}
NO_ACCESS = {"private", "no"}
NO_SERVICE = {"parking_aisle", "driveway", "drive-through"}


def classify(tags):
    """The app's highway type for a way, or None if not walkable.

    Mirrors OverpassClient.streetsIn / connectorQuery and OsmData._highway.
    """
    highway = tags.get("highway")
    if highway is None or tags.get("access") in NO_ACCESS:
        return None
    foot = tags.get("foot")
    if highway == "cycleway":
        if foot in ("yes", "designated"):
            return "cycleway"
        if foot in ("no", "private"):
            return None
        return "cycleway_connector"
    if highway not in STREET_HIGHWAYS or foot in ("no", "private"):
        return None
    footway = tags.get("footway")
    if highway == "footway" and footway in ("sidewalk", "crossing"):
        return footway
    if tags.get("service") in NO_SERVICE or tags.get("area") == "yes":
        return None
    return highway


def tile_of(lat, lon):
    n = 1 << ZOOM
    x = int((lon + 180) / 360 * n)
    lat_rad = math.radians(lat)
    y = int((1 - math.log(math.tan(lat_rad) + 1 / math.cos(lat_rad)) / math.pi) / 2 * n)
    return min(max(x, 0), n - 1), min(max(y, 0), n - 1)


def tile_bounds(x, y):
    n = 1 << ZOOM

    def lat(yy):
        return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * yy / n))))

    return x / n * 360 - 180, lat(y + 1), (x + 1) / n * 360 - 180, lat(y)


def encode_polyline6(coords):
    """Google polyline with 6 decimals, (lat, lon) order like Valhalla's."""
    out = io.StringIO()
    last_lat = last_lon = 0
    for lon, lat in coords:
        ilat, ilon = round(lat * 1e6), round(lon * 1e6)
        for delta in (ilat - last_lat, ilon - last_lon):
            value = ~(delta << 1) if delta < 0 else delta << 1
            while value >= 0x20:
                out.write(chr((0x20 | (value & 0x1F)) + 63))
                value >>= 5
            out.write(chr(value + 63))
        last_lat, last_lon = ilat, ilon
    return out.getvalue()


def tiles_along(coords):
    """Zoom-13 tiles a line touches: those of its nodes, plus points every
    ~100 m in between so long straight stretches don't skip a tile."""
    tiles = set()
    for i, (lon, lat) in enumerate(coords):
        tiles.add(tile_of(lat, lon))
        if i == 0:
            continue
        plon, plat = coords[i - 1]
        meters = math.hypot(
            (lon - plon) * 111195 * math.cos(math.radians(lat)),
            (lat - plat) * 111195,
        )
        steps = int(meters // 100)
        for s in range(1, steps + 1):
            t = s / (steps + 1)
            tiles.add(tile_of(plat + (lat - plat) * t, plon + (lon - plon) * t))
    return tiles


def read_poly(path):
    """Geofabrik .poly file -> shapely polygon. Sections whose name starts
    with "!" are holes."""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()][1:]
    outer, holes = [], []
    i = 0
    while i < len(lines) and lines[i] != "END":
        hole = lines[i].startswith("!")
        ring = []
        i += 1
        while lines[i] != "END":
            lon, lat = map(float, lines[i].split()[:2])
            ring.append((lon, lat))
            i += 1
        i += 1
        if len(ring) >= 3:
            (holes if hole else outer).append(Polygon(ring))
    shape = outer[0]
    for r in outer[1:]:
        shape = shape.union(r)
    for h in holes:
        shape = shape.difference(h)
    return shape


def piece_meters(a, b):
    (lon1, lat1), (lon2, lat2) = a, b
    x = (lon2 - lon1) * math.cos(math.radians((lat1 + lat2) / 2))
    return math.hypot(x, lat2 - lat1) * 111195


def write_ways(pbf, out, street_meters):
    """Pass over the extract: every walkable way, once per tile it touches,
    as a line ``x,y<TAB>json`` (sorted afterwards). Adds the street length
    per tile to [street_meters]."""
    count = 0
    processor = (
        osmium.FileProcessor(pbf)
        .with_locations()
        .with_filter(osmium.filter.KeyFilter("highway"))
    )
    for obj in processor:
        if not obj.is_way():
            continue
        highway = classify(obj.tags)
        if highway is None:
            continue
        try:
            coords = [(n.lon, n.lat) for n in obj.nodes]
        except osmium.InvalidLocationError:
            continue
        if len(coords) < 2:
            continue
        record = json.dumps(
            [
                obj.id,
                highway,
                obj.tags.get("name"),
                [n.ref for n in obj.nodes],
                encode_polyline6(coords),
            ],
            ensure_ascii=False,
            separators=(",", ":"),
        )
        for x, y in tiles_along(coords):
            out.write(f"{x:05d},{y:05d}\t{record}\n")
        if highway not in CONNECTORS:
            for a, b in zip(coords, coords[1:]):
                key = tile_of((a[1] + b[1]) / 2, (a[0] + b[0]) / 2)
                street_meters[key] = street_meters.get(key, 0) + piece_meters(a, b)
        count += 1
    return count


def build(pbf, poly_path, out_dir, name, tag):
    os.makedirs(out_dir, exist_ok=True)
    boundary = read_poly(poly_path)
    inside = prep(boundary)
    edge = prep(boundary.boundary)

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        # Compressed throughout: for a whole country the lines take several
        # times the disk space of the extract.
        raw = os.path.join(tmp, "ways.tsv.gz")
        with gzip.open(raw, "wt", encoding="utf-8", compresslevel=1) as f:
            street_meters = {}
            ways = write_ways(pbf, f, street_meters)
        print(f"{ways} walkable ways", flush=True)
        sort = subprocess.Popen(
            f"gzip -dc '{raw}' | sort -t \"$(printf '\\t')\" -k1,1 "
            f"-S {os.environ.get('SORT_MEMORY', '1G')} "
            f"--compress-program=gzip -T '{tmp}'",
            shell=True,
            stdout=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        sorted_lines = io.TextIOWrapper(sort.stdout, encoding="utf-8")

        entries = {}
        parts = []
        part_file = None
        offset = 0

        def open_part():
            nonlocal part_file, offset
            if part_file:
                part_file.close()
            parts.append(f"{name}-{len(parts) + 1}.bin")
            part_file = open(os.path.join(out_dir, parts[-1]), "wb")
            offset = 0

        def flush(key, lines):
            nonlocal offset
            if key is None:
                return
            x, y = map(int, key.split(","))
            body = '{"w":[' + ",".join(lines) + "]}"
            data = gzip.compress(body.encode("utf-8"), compresslevel=9, mtime=0)
            if part_file is None or offset + len(data) > PART_LIMIT:
                open_part()
            part_file.write(data)
            entries[(x, y)] = (len(parts) - 1, 0, offset, len(data))
            offset += len(data)

        key, lines = None, []
        for line in sorted_lines:
            k, record = line.rstrip("\n").split("\t", 1)
            if k != key:
                flush(key, lines)
                key, lines = k, []
            lines.append(record)
        flush(key, lines)
        if sort.wait() != 0:
            sys.exit("sort failed")
        if part_file:
            part_file.close()

    # Every tile inside the extract gets an entry, empty ones too, so the
    # app knows there is nothing to load; tiles on the border are flagged.
    west, south, east, north = boundary.bounds
    x0, y0 = tile_of(north, west)
    x1, y1 = tile_of(south, east)
    # Ways reaching beyond the extract put data into tiles outside it.
    for key, (part, _, off, length) in list(entries.items()):
        entries[key] = (part, FLAG_PARTIAL, off, length)
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            b = box(*tile_bounds(x, y))
            if not inside.intersects(b):
                continue
            flags = FLAG_PARTIAL if edge.intersects(b) else 0
            part, _, off, length = entries.get((x, y), (0, FLAG_EMPTY, 0, 0))
            if length == 0:
                flags |= FLAG_EMPTY
            entries[(x, y)] = (part, flags, off, length)

    index_name = f"{name}-index.bin.gz"
    index = b"".join(
        INDEX_ENTRY.pack(
            x, y, *entries[(x, y)], round(street_meters.get((x, y), 0))
        )
        for (x, y) in sorted(entries)
    )
    with open(os.path.join(out_dir, index_name), "wb") as f:
        f.write(gzip.compress(index, compresslevel=9, mtime=0))

    manifest = {
        "format": FORMAT,
        "zoom": ZOOM,
        "tag": tag,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "index": index_name,
        "parts": parts,
        "bbox": [west, south, east, north],
    }
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    size = sum(os.path.getsize(os.path.join(out_dir, p)) for p in parts)
    print(f"{len(entries)} tiles, {len(parts)} part(s), {size / 1e6:.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    build(*sys.argv[1:])
