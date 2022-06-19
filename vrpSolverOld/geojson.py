import geojson

def readHighwaysFromGeoJSON(
    geoJSONPath:    "Path to geojson file" = None
    ) -> "Given a geoJSON file, returns a list of network dictionaries":

    # Open the geojson file ===================================================
    with open(geoJSONPath, encoding = 'utf-8') as f:
        gj = geojson.load(f)

    # Categorize highways =====================================================
    motorway = {}
    country = {}
    residential = {}
    unknown = {}

    roadID = 0
    for i in range(len(gj["features"])):
        if (gj["features"][i]["geometry"]["type"] == "LineString" 
            and "properties" in gj["features"][i] 
            and "highway" in gj["features"][i]["properties"]):
            if (gj["features"][i]["properties"]["highway"] in ['motorway', 'motorway_link']):
                line = []
                for p in gj["features"][i]["geometry"]["coordinates"]:
                    line.append([p[1], p[0]])
                motorway[roadID] = {}
                motorway[roadID]['line'] = line
                if ("name" in gj["features"][i]["properties"]):
                    motorway[roadID]['name'] = gj["features"][i]["properties"]["name"]
                if ("maxspeed" in gj["features"][i]["properties"]):
                    motorway[roadID]['maxspeed'] = gj["features"][i]["properties"]["maxspeed"]
                if ("oneway" in gj["features"][i]["properties"]):
                    motorway[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]
                roadID += 1
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
                line = []
                for p in gj["features"][i]["geometry"]["coordinates"]:
                    line.append([p[1], p[0]])
                country[roadID] = {}
                country[roadID]['line'] = line
                if ("name" in gj["features"][i]["properties"]):
                    country[roadID]['name'] = gj["features"][i]["properties"]["name"]
                if ("maxspeed" in gj["features"][i]["properties"]):
                    country[roadID]['maxspeed'] = gj["features"][i]["properties"]["maxspeed"]
                if ("oneway" in gj["features"][i]["properties"]):
                    country[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]
                roadID += 1
            elif (gj["features"][i]["properties"]["highway"] in [
                'residential']):
                line = []
                for p in gj["features"][i]["geometry"]["coordinates"]:
                    line.append([p[1], p[0]])
                residential[roadID] = {}
                residential[roadID]['line'] = line
                if ("name" in gj["features"][i]["properties"]):
                    residential[roadID]['name'] = gj["features"][i]["properties"]["name"]
                if ("maxspeed" in gj["features"][i]["properties"]):
                    residential[roadID]['maxspeed'] = gj["features"][i]["properties"]["maxspeed"]
                if ("oneway" in gj["features"][i]["properties"]):
                    residential[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]
                roadID += 1
            else:
                line = []
                for p in gj["features"][i]["geometry"]["coordinates"]:
                    line.append([p[1], p[0]])
                unknown[roadID] = {}
                unknown[roadID]['line'] = line
                if ("name" in gj["features"][i]["properties"]):
                    unknown[roadID]['name'] = gj["features"][i]["properties"]["name"]
                if ("maxspeed" in gj["features"][i]["properties"]):
                    unknown[roadID]['maxspeed'] = gj["features"][i]["properties"]["maxspeed"]
                if ("oneway" in gj["features"][i]["properties"]):
                    unknown[roadID]['oneway'] = gj["features"][i]["properties"]["oneway"]
                roadID += 1
    return {
        'motorway': motorway,
        'country': country,
        'residential': residential,
        'unknown': unknown
    }