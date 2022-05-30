import geojson
import math

from .const import *

def createRoadNetworkFromGeoJSON(
    geoJSONPath:    "Path to geojson file" = None
    ) -> "Given a geoJSON file, returns a list of network dictionaries":

    # Open the geojson file ===================================================
    with open(geoJSONPath, encoding = 'utf-8') as f:
        gj = geojson.load(f)

    # Create roads from geoJSONN ==============================================
    road = {}
    roadID = 0
    # Loop through data
    for i in range(len(gj["features"])):
        if (gj["features"][i]["geometry"]["type"] == "LineString" 
            and "properties" in gj["features"][i] 
            and "highway" in gj["features"][i]["properties"]):
            # The shape of the road is projected using Mercator projection
            line = []
            for p in gj["features"][i]["geometry"]["coordinates"]:
                line.append(list(_ptLatLon2XYMercator([p[1], p[0]])))
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
            elif (gj["features"][i]["properties"]["highway"] in [
                'truck', 
                'truck_link', 
                'primary', 
                'primary_link', 
                'secondary',
                'secondary_link',
                'tertiary',
                'tertiary_link',
                'unclassified']):
                road[roadID]['class'] = 'country'
            elif (gj["features"][i]["properties"]["highway"] in [
                'residential']):
                road[roadID]['class'] = 'residential'
            else:
                road[roadID]['class'] = 'unknown'
            roadID += 1
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
    for r in road:
        for p in road[r]['shape']:
            p[0] -= minX
            p[1] -= minY

    # Define boundary =========================================================
    boundary = [(0, 0), (maxX - minX, 0), (maxX - minX, maxY - minY), (0, maxY - minY)]

    return {
        'boundary': boundary,
        'road': road        
    }

def roadNetwork2Arcs(
    road:       "Road network dictionary" = None
    ) -> "Given a road network, convert it into a list of arcs and a list of nodes":

    return {
        'node': node,
        'arc': arc
    }