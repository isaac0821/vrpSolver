import geojson
import math
import shapely
import networkx as nx

from .geometry import *

def createRoadNetworkFromGeoJSON(
    geoJSONPath: str,
    boundaryLatLon: poly|None = None,
    projType: str = 'LatLon',
    buildingIncludedFlag: bool = False
    ) -> dict:

    """Given a path to geoJSON file, returns a road network dictionary with roads (and buildings)


    Parameters
    ----------
    
    geoJSONPath: string, required
        The path to a geoJSON file (".geojson")
    boundaryLatLon: list[pt], optionan, default as None
        Filter out the area that are not within the boundary, if given
    projType: string, optional, default as 'LatLon'
        Determine the axises of exported data, options are ['LatLon', 'Mercator']
    buildingIncludedFlag: bool, optional, default as False
        True if buildings are included

    Returns
    -------

    dict
        A road network dictionary, in the following formatt:
            >>> road = {
            ...     'boundary': boundary,
            ...     'road': road,
            ...     'building': building
            ... }
    """

    # Open the geojson file ===================================================
    with open(geoJSONPath, encoding = 'utf-8') as f:
        gj = geojson.load(f)

    # Create roads from geoJSON ===============================================
    road = {}
    roadID = 0

    # Create buildings from geoJSON ===========================================
    building = {}
    buildingID = 0

    # Loop through data
    # FIXME: Need to truncate the roads on the edge of polygon
    for i in range(len(gj["features"])):
        if (gj["features"][i]["geometry"]["type"] == "LineString" 
            and "properties" in gj["features"][i] 
            and "highway" in gj["features"][i]["properties"]):

            roadEntireIncludedFlag = False   # All nodes should be included
            roadPartialIncludedFlag = False  # Some (but not all) nodes should be included

            # The shape of the road is projected using Mercator projection
            # FIXME: truncate the roads that partially inside polygon
            line = []
            for p in gj["features"][i]["geometry"]["coordinates"]:
                if (boundaryLatLon == None or isPtOnPoly([p[1], p[0]], boundaryLatLon)):
                    roadEntireIncludedFlag = True
                if (projType == 'Mercator'):
                    line.append(list(ptLatLon2XYMercator([p[1], p[0]])))
                elif (projType == 'LatLon'):
                    line.append([p[1], p[0]])

            if (roadEntireIncludedFlag):
                road[roadID] = {}
                road[roadID]['shape'] = line

                # Name of road
                if ("name" in gj["features"][i]["properties"]):
                    road[roadID]['name'] = gj["features"][i]["properties"]["name"]

                # Speed of road
                if ("maxspeed" in gj["features"][i]["properties"]):
                    road[roadID]['maxspeed'] = gj["features"][i]["properties"]["maxspeed"]
                else:
                    road[roadID]['maxspeed'] = 30 # [mph]

                # Is it one-way
                if ("oneway" in gj["features"][i]["properties"]):
                    road[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]    
                else:
                    road[roadID]['oneway'] = False

                # Road class
                road[roadID]['class'] = gj["features"][i]["properties"]["highway"]
                
                roadID += 1

        if (buildingIncludedFlag):
            if (gj["features"][i]["geometry"]["type"] in ["MultiPolygon"]
                and "properties" in gj["features"][i]
                and "building" in gj["features"][i]["properties"]):

                # The shape of the road is projected using Mercator projection
                poly = []
                if (projType == 'Mercator'):
                    for p in gj["features"][i]["geometry"]["coordinates"][0][0]:
                        poly.append(list(ptLatLon2XYMercator([p[1], p[0]])))
                elif (projType == 'LatLon'):
                    for p in gj["features"][i]["geometry"]["coordinates"][0][0]:
                        poly.append([p[1], p[0]])
                building[buildingID] = {}
                building[buildingID]['shape'] = poly

                # Building type
                if (gj["features"][i]["properties"]["building"] == 'yes'):
                    building[buildingID]['type'] = "building"
                else:
                    building[buildingID]['type'] = gj["features"][i]["properties"]["building"]

                buildingID += 1

    # Fix the position ========================================================
    minX = None
    maxX = None
    minY = None
    maxY = None
    for r in road:
        for p in road[r]['shape']:
            if (minX == None or p[0] < minX):
                minX = p[0]
            if (maxX == None or p[0] > maxX):
                maxX = p[0]
            if (minY == None or p[1] < minY):
                minY = p[1]
            if (maxY == None or p[1] > maxY):
                maxY = p[1]
    if (projType == 'Mercator'):
        for r in road:
            for p in road[r]['shape']:
                p[0] -= minX
                p[1] -= minY

    # Define boundary =========================================================
    boundary = []
    if (projType == 'Mercator'):
        boundary = [(0, 0), (maxX - minX, 0), (maxX - minX, maxY - minY), (0, maxY - minY)]
    elif (projType == 'LatLon'):
        boundary = [(minX, minY), (maxX, minY), (maxX, maxY), (minX, maxY)]

    return {
        'boundary': boundary,
        'projType': projType,
        'road': road,
        'building': building
    }

def clipRoadsByPoly(
    roads: dict,
    poly: poly) -> dict:

    # FIXME: Currently using a stupid method, since it is a one-time function    
    # Roads ===================================================================
    clip = {}
    maxRoadID = max(roads.keys()) + 1
    for r in roads:
        # 最笨的办法，一条一条路处理
        seqInt = intSeq2Poly(roads[r]['shape'], poly)

        # 如果是完整的路径
        if (type(seqInt) != list):
            if (seqInt['status'] == 'Cross' and seqInt['intersectType'] == 'Segment'):
                clip[r] = {}
                for k in roads[r]:
                    clip[r][k] = roads[r][k]
                clip[r]['shape'] = seqInt['intersect']

        # 如果路径被切割开
        else:
            for parti in seqInt:
                if (parti['status'] == 'Cross' and parti['intersectType'] == 'Segment'):
                    clip[maxRoadID] = {}
                    for k in roads[r]:
                        clip[maxRoadID][k] = roads[r][k]
                    clip[maxRoadID]['shape'] = parti['intersect']
                    maxRoadID += 1

    return clip

def clipRoadsByMultiPoly(
    roads: dict,
    multiPoly: polys) -> dict:

    # FIXME: Currently using a stupid method, since it is a one-time function    
    # Roads ===================================================================
    clip = {}
    roadID = 0
    for poly in multiPoly:
        for r in roads:
            # 最笨的办法，一条一条路处理
            seqInt = intSeq2Poly(roads[r]['shape'], poly)

            # 如果是完整的路径
            if (type(seqInt) != list):
                if (seqInt['status'] == 'Cross' and seqInt['intersectType'] == 'Segment'):
                    clip[roadID] = {}
                    for k in roads[r]:
                        clip[roadID][k] = roads[r][k]
                    clip[roadID]['shape'] = seqInt['intersect']
                    roadID += 1

            # 如果路径被切割开
            else:
                for parti in seqInt:
                    if (parti['status'] == 'Cross' and parti['intersectType'] == 'Segment'):
                        clip[roadID] = {}
                        for k in roads[r]:
                            clip[roadID][k] = roads[r][k]
                        clip[roadID]['shape'] = parti['intersect']
                        roadID += 1

    return clip