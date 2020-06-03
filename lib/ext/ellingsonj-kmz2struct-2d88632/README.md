# kmz2struct
This function will convert kml and kmz files to a matlab structure. If converting a kmz file it will extract it to a directory called '.kml2struct' in your home directory. This directory will be deleted when the function exits.

The output of this function should be similar to 'shaperead' except that it will add another field, "Folder", for the kml folder where the shape was file.

This function will only handle kml/kmz files with Point, LineString, and Polygon geometries. If you try to run this on a kml/kmz with different elements those elements will be omitted from the result
