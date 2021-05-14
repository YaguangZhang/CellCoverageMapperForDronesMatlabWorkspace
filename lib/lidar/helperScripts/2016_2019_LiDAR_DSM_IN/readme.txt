Note: to include buildings and trees, we need the point cloud version of the LiDAR data.

From Dr. Sam Noel:
"DEM is generally used for ground elevations and applications that would prefer that"
"DTM (digital terrain model) is more general and is not necessarily just the ground"

Here we have the 2016-2019 LiDAR data set.

-----------------------------------
About Loading the Data into PostGIS
-----------------------------------
SRID (Spatial Reference IDentifier)

Install GDAL (e.g., via OSGeo4W64) and run:

    `gdalinfo pathToTifFile`

The projection information obtained this way can be used to get the SRID:

    NAD_1983_HARN_StatePlane_Indiana_East_FIPS_1301: 2792
    NAD_1983_HARN_StatePlane_Indiana_West_FIPS_1302: 2793
    NAD_1983_HARN_StatePlane_Indiana_East_FIPS_1301_Feet: 2967
    NAD_1983_HARN_StatePlane_Indiana_West_FIPS_1302_Feet: 2968

Then run

    `raster2pgsql -s 2968 -t 128x128 "D:\One Drive - Purdue\OneDrive - purdue.edu\OATS\CellCoverageMapper\Lidar_2019\IN\DSM_Demo_ACRE\Raw\in2018_29701895_12_dsm.tif" -I -C -M lidar_z | psql "sslmode=allow host=localhost port=5433 dbname=lidar_in user=postgres"`

Refs:
    https://www.compose.com/articles/geofile-postgis-and-raster-data/
    https://developers.arcgis.com/rest/services-reference/projected-coordinate-systems.htm
    https://github.com/lcalisto/workshop-postgis-raster