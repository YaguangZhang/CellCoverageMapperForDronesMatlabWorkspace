#!/bin/bash
# Please see BulkDownloadAllCountiesDsm.sh for more details.
mkdir -p '/mnt/y/data/CellCoverageMapper/Lidar_2019/IN/'
cd '/mnt/y/data/CellCoverageMapper/Lidar_2019/IN/'
pwd
read -n 1 -s -r -p "Press any key to continue"
echo " "
echo " "
echo "====== Downloading Cloud Point Data ======"
mkdir -p "./CloudPoint/QL2_3DEP_LiDAR_IN_2017_2019_laz"
wget --no-parent -e robots=off -U mozilla -A '*.laz' -r -nd -N -P "./CloudPoint/QL2_3DEP_LiDAR_IN_2017_2019_laz" "https://lidar.jinha.org/QL2_3DEP_LiDAR_IN_2017_2019_laz/"
echo " "
echo "====== Done ======"