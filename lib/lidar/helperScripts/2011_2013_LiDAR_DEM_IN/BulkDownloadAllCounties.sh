#!/bin/bash
# Bulk download the LiDAR data from OpenTopo for all IN counties.

mkdir -p '/mnt/d/One Drive - Purdue/OneDrive - purdue.edu/OATS/CellCoverageMapper/Lidar/IN'
cd '/mnt/d/One Drive - Purdue/OneDrive - purdue.edu/OATS/CellCoverageMapper/Lidar/IN'
pwd
read -n 1 -s -r -p "Press any key to continue"
echo " "
echo " "
echo "====== Downloading ======"
datasets=(
    IN_2011_2013_W
    IN_2011_2013_E
)
for d in "${datasets[@]}"; do
    mkdir -p "./${d}"
    wget --no-parent -r -nd -N -P "./${d}" "https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/${d}/"
done
echo " "
echo "====== Done ======"