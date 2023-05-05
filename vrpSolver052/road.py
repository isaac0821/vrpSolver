import geojson
import math

from .const import *
from .geometry import *

def createRoadNetworkFromGeoJSON(
    geoJSONPath:    "Path to geojson file" = None,
    boundaryLatLon: "Filter out the area that are not within the boundary, if given" = None,
    exportType:     "1) String, 'LatLon', exported as in lat/lon, or\
                     2) String, 'Mercator', exported as projected in Mercator" = 'LatLon',
    buildingIncludedFlag: "True if buildings are included" = False
    ) -> "Given a geoJSON file, returns a road network dictionaries which only includes road info":

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
                if (exportType == 'Mercator'):
                    line.append(list(ptLatLon2XYMercator([p[1], p[0]])))
                elif (exportType == 'LatLon'):
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

                # Is it one-way
                if ("oneway" in gj["features"][i]["properties"]):
                    road[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]    

                # Road class
                if (gj["features"][i]["properties"]["highway"] in ['motorway', 'motorway_link']):
                    road[roadID]['class'] = 'motorway'
                elif (gj["features"][i]["properties"]["highway"] in ['truck', 'truck_link']):
                    road[roadID]['class'] = 'truck'
                elif (gj["features"][i]["properties"]["highway"] in ['primary', 'primary_link']):
                    road[roadID]['class'] = 'primary'
                elif (gj["features"][i]["properties"]["highway"] in ['secondary', 'secondary_link']):
                    road[roadID]['class'] = 'secondary'
                elif (gj["features"][i]["properties"]["highway"] in ['tertiary', 'tertiary_link']):
                    road[roadID]['class'] = 'tertiary'
                elif (gj["features"][i]["properties"]["highway"] in ['residential']):
                    road[roadID]['class'] = 'residential'
                else:
                    road[roadID]['class'] = gj["features"][i]["properties"]["highway"]
                
                roadID += 1

        if (buildingIncludedFlag):
            if (gj["features"][i]["geometry"]["type"] in ["MultiPolygon"]
                and "properties" in gj["features"][i]
                and "building" in gj["features"][i]["properties"]):

                # The shape of the road is projected using Mercator projection
                poly = []
                if (exportType == 'Mercator'):
                    for p in gj["features"][i]["geometry"]["coordinates"][0][0]:
                        poly.append(list(ptLatLon2XYMercator([p[1], p[0]])))
                elif (exportType == 'LatLon'):
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

    if (exportType == 'Mercator'):
        for r in road:
            for p in road[r]['shape']:
                p[0] -= minX
                p[1] -= minY

    # Define boundary =========================================================
    if (exportType == 'Mercator'):
        boundary = [(0, 0), (maxX - minX, 0), (maxX - minX, maxY - minY), (0, maxY - minY)]
    elif (exportType == 'LatLon'):
        boundary = [(minX, minY), (maxX, minY), (maxX, maxY), (minX, maxY)]

    return {
        'boundary': boundary,
        'road': road,
        'building': building
    }
