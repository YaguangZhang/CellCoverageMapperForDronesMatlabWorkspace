#!/bin/bash
# Bulk download the LiDAR data from OpenTopo for selected counties.
#
# Here we use the counties within Wabash Heartland Innovation Network.

mkdir -p '/mnt/d/One Drive - Purdue/OneDrive - purdue.edu/OATS/CellCoverageMapper/Lidar/Tipp_Extended/IN_2011_2013_W'
cd '/mnt/d/One Drive - Purdue/OneDrive - purdue.edu/OATS/CellCoverageMapper/Lidar/Tipp_Extended/IN_2011_2013_W'
pwd
read -n 1 -s -r -p "Press any key to continue"
echo " "
echo " "
echo "====== Downloading ======"
counties=(
    Pulaski
    White
    Cass
    Benton
    Carroll
    Tippecanoe
    Warren
    Fountain
    Montgomery
    Clinton 
)
for c in "${counties[@]}"; do
    mkdir -p "./${c}"
    wget --no-parent -r -nd -N -P "./${c}" "https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/IN_2011_2013_W/IN_2011_2013_W_be/${c}/"
done
echo " "
echo "====== Done ======"