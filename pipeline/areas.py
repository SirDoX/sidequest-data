"""Administrative areas for the Sidequest statistics.

From an extract with boundaries and place nodes this writes

* ``<name>-areas.bin``           one gzip JSON record per area:
  ``{"id", "l": admin level, "v": is a village, "n": names, "m": street
  metres inside, "o": [polyline6 rings]}``
* ``<name>-areas-index.bin.gz``  id (int64), offset, length (uint32) per area
* ``<name>-chunk-areas.bin.gz``  per zoom-13 chunk: x, y, count (uint16) and
  the ids (int64) of the areas touching it

Villages: a part of a municipality (level 9/10) is a village when a town,
village or hamlet of the same name lies inside. Where a municipality has
several villages without boundaries of their own, it is split into cells
around their place nodes (Voronoi); those areas get the negative node id.
Cities are never split.
"""

import gzip
import json
import struct
from array import array

import numpy as np
import osmium
import shapely
import shapely.wkb
from shapely.geometry import MultiPoint, Point, box
from shapely.prepared import prep

AREA_INDEX = struct.Struct("<qII")
CHUNK_HEAD = struct.Struct("<HHH")
SETTLEMENTS = {"town", "village", "hamlet"}

# Simplification in degrees: ~5 m; ~20 m for states (Berlin and Hamburg
# among them, which are also shown as places); ~100 m for countries.
SIMPLIFY = 0.00005
SIMPLIFY_STATE = 0.0002
SIMPLIFY_COUNTRY = 0.001


class Pieces:
    """Street pieces (between two nodes) per chunk: middle and length."""

    def __init__(self):
        self.by_chunk = {}

    def add(self, chunk, lon, lat, meters):
        arrays = self.by_chunk.get(chunk)
        if arrays is None:
            arrays = self.by_chunk[chunk] = (array("f"), array("f"), array("f"))
        arrays[0].append(lon)
        arrays[1].append(lat)
        arrays[2].append(meters)


def _names(tags):
    return {t.k: t.v for t in tags if t.k == "name" or t.k.startswith("name:")}


def read_admin(pbf):
    """Boundaries (relations with admin_level 2-10) and settlement nodes."""
    factory = osmium.geom.WKBFactory()
    areas, places = [], []
    processor = (
        osmium.FileProcessor(pbf)
        .with_areas(osmium.filter.KeyFilter("boundary"))
        .with_filter(osmium.filter.KeyFilter("boundary", "place"))
    )
    for obj in processor:
        if obj.is_area():
            level = obj.tags.get("admin_level", "")
            if (
                obj.from_way()
                or obj.tags.get("boundary") != "administrative"
                or not level.isdigit()
                or not 2 <= int(level) <= 10
            ):
                continue
            try:
                geom = shapely.wkb.loads(factory.create_multipolygon(obj), hex=True)
            except Exception:
                continue
            if not geom.is_valid:
                # Self-touching rings happen in OSM; repair rather than lose
                # the area.
                geom = shapely.make_valid(geom)
            if geom.is_empty:
                continue
            areas.append(
                {
                    "id": obj.orig_id(),
                    "l": int(level),
                    "v": False,
                    "n": _names(obj.tags),
                    "g": geom,
                }
            )
        elif obj.is_node() and obj.tags.get("place") in SETTLEMENTS | {"city"}:
            if "name" not in obj.tags:
                continue
            places.append(
                {
                    "id": obj.id,
                    "kind": obj.tags["place"],
                    "n": _names(obj.tags),
                    "p": Point(obj.location.lon, obj.location.lat),
                }
            )
    return areas, places


def mark_villages(areas, places):
    """Parts of a municipality with a settlement of the same name inside."""
    tree = shapely.STRtree([p["p"] for p in places])
    for area in areas:
        if area["l"] < 9:
            continue
        name = area["n"].get("name")
        for i in tree.query(area["g"], predicate="contains"):
            place = places[i]
            if place["kind"] in SETTLEMENTS and place["n"].get("name") == name:
                area["v"] = True
                break


def synthesize_villages(areas, places):
    """Cells around the villages of municipalities that have several but
    no village boundaries."""
    tree = shapely.STRtree([p["p"] for p in places])
    village_parts = [a for a in areas if a["v"]]
    parts_tree = shapely.STRtree([a["g"] for a in village_parts]) if village_parts else None
    created = []
    for town in [a for a in areas if a["l"] == 8]:
        geom = town["g"]
        inside = [places[i] for i in tree.query(geom, predicate="contains")]
        # A city stays whole, however many villages it swallowed.
        if any(p["kind"] == "city" for p in inside):
            continue
        settlements = [p for p in inside if p["kind"] in SETTLEMENTS]
        if len(settlements) < 2:
            continue
        if parts_tree is not None and len(parts_tree.query(geom, predicate="intersects")):
            # Has village boundaries of its own.
            continue
        cells = shapely.voronoi_polygons(
            MultiPoint([p["p"] for p in settlements]), extend_to=geom
        )
        for cell in cells.geoms:
            place = next((p for p in settlements if cell.contains(p["p"])), None)
            if place is None:
                continue
            shape = cell.intersection(geom)
            if shape.is_empty:
                continue
            created.append(
                {"id": -place["id"], "l": 10, "v": True, "n": place["n"], "g": shape}
            )
    return created


def chunks_touching(geom, tile_of, tile_bounds):
    """Zoom-13 chunks overlapping the geometry, with whether they lie
    wholly inside."""
    west, south, east, north = geom.bounds
    x0, y0 = tile_of(north, west)
    x1, y1 = tile_of(south, east)
    shape = prep(geom)
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            b = box(*tile_bounds(x, y))
            if not shape.intersects(b):
                continue
            yield (x, y), shape.contains(b)


def street_meters(geom, pieces, chunk_meters, tile_of, tile_bounds):
    """Length of the street pieces whose middle lies inside."""
    total = 0.0
    for chunk, whole in chunks_touching(geom, tile_of, tile_bounds):
        if whole:
            total += chunk_meters.get(chunk, 0)
            continue
        arrays = pieces.by_chunk.get(chunk)
        if arrays is None:
            continue
        lons = np.frombuffer(arrays[0], dtype=np.float32).astype(np.float64)
        lats = np.frombuffer(arrays[1], dtype=np.float32).astype(np.float64)
        meters = np.frombuffer(arrays[2], dtype=np.float32)
        total += float(meters[shapely.contains_xy(geom, lons, lats)].sum())
    return total


def rings(geom):
    polygons = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
    for polygon in polygons:
        if polygon.geom_type != "Polygon":
            continue
        yield list(polygon.exterior.coords)
        for interior in polygon.interiors:
            yield list(interior.coords)


def write_areas(admin_pbf, pieces, chunk_meters, out_dir, name, helpers):
    """Writes the three area files; returns their names for the manifest."""
    tile_of, tile_bounds, encode_polyline6 = helpers
    areas, places = read_admin(admin_pbf)
    mark_villages(areas, places)
    areas += synthesize_villages(areas, places)
    print(f"{len(areas)} areas, {sum(a['v'] for a in areas)} villages", flush=True)

    data_name = f"{name}-areas.bin"
    index_name = f"{name}-areas-index.bin.gz"
    chunks_name = f"{name}-chunk-areas.bin.gz"
    by_chunk = {}
    index = []
    offset = 0
    with open(f"{out_dir}/{data_name}", "wb") as data:
        for area in sorted(areas, key=lambda a: a["id"]):
            geom = area["g"]
            for chunk, _ in chunks_touching(geom, tile_of, tile_bounds):
                by_chunk.setdefault(chunk, []).append(area["id"])
            tolerance = (
                SIMPLIFY_COUNTRY
                if area["l"] <= 2
                else SIMPLIFY_STATE if area["l"] <= 4 else SIMPLIFY
            )
            simple = geom.simplify(tolerance, preserve_topology=True)
            record = {
                "id": area["id"],
                "l": area["l"],
                "v": area["v"],
                "n": area["n"],
                "m": round(street_meters(geom, pieces, chunk_meters, tile_of, tile_bounds)),
                "o": [encode_polyline6(r) for r in rings(simple)],
            }
            blob = gzip.compress(
                json.dumps(record, ensure_ascii=False, separators=(",", ":")).encode(),
                compresslevel=9,
                mtime=0,
            )
            data.write(blob)
            index.append(AREA_INDEX.pack(area["id"], offset, len(blob)))
            offset += len(blob)

    with open(f"{out_dir}/{index_name}", "wb") as f:
        f.write(gzip.compress(b"".join(index), compresslevel=9, mtime=0))
    chunk_blob = bytearray()
    for (x, y) in sorted(by_chunk):
        ids = by_chunk[(x, y)]
        chunk_blob += CHUNK_HEAD.pack(x, y, len(ids))
        chunk_blob += struct.pack(f"<{len(ids)}q", *ids)
    with open(f"{out_dir}/{chunks_name}", "wb") as f:
        f.write(gzip.compress(bytes(chunk_blob), compresslevel=9, mtime=0))
    print(f"areas: {offset / 1e6:.1f} MB, {len(by_chunk)} chunks", flush=True)
    return {"data": data_name, "index": index_name, "chunks": chunks_name}
