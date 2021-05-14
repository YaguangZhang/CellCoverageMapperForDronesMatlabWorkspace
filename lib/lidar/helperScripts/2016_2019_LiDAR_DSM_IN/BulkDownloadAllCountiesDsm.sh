#!/bin/bash
# Bulk download the LiDAR data from OpenTopo for all IN counties.
#
# Ref dir: '/mnt/d/One Drive - Purdue/OneDrive - purdue.edu/OATS/CellCoverageMapper/Lidar_2019/'
# Ref dir for Data Depot: '/mnt/y/data/CellCoverageMapper/Lidar_2019/IN/'
#    One would need to mount the network drive (e.g., Y) first via something like:
#       sudo mkdir -p /mnt/y
#       sudo mount -t drvfs Y: /mnt/y
# Notes:
#    (1) To save all files under one directory without preserving the online folder structures
#           wget --no-parent -e robots=off -U mozilla -A '*_dsm.tif' -r -nd -N -P "dir" "url"
#        e.g.,
#           wget --no-parent -e robots=off -U mozilla -A '*_dsm.tif' -r -nd -N -P "./DSM/QL2_3DEP_LiDAR_IN_2017_2019_l2" "https://lidar.jinha.org/QL2_3DEP_LiDAR_IN_2017_2019_l2/"
#    (2) To save files following the online folder structures
#           wget --no-parent -e robots=off -U mozilla -A '*_dsm.tif' -r -x -nH -N -P "dir" "url"
#        e.g.,
#           wget --no-parent -e robots=off -U mozilla -A '*_dsm.tif' -r -x -nH -N -P "./DSM/" "https://lidar.jinha.org/QL2_3DEP_LiDAR_IN_2017_2019_l2/"
mkdir -p '/mnt/y/data/CellCoverageMapper/Lidar_2019/IN/'
cd '/mnt/y/data/CellCoverageMapper/Lidar_2019/IN/'
pwd
read -n 1 -s -r -p "Press any key to continue"
echo " "
echo " "
echo "====== Downloading DSM Data ======"
mkdir -p "./DSM/QL2_3DEP_LiDAR_IN_2017_2019_l2"
wget --no-parent -e robots=off -U mozilla -A '*_dsm.tif' -r -nd -N -P "./DSM/QL2_3DEP_LiDAR_IN_2017_2019_l2" "https://lidar.jinha.org/QL2_3DEP_LiDAR_IN_2017_2019_l2/"
echo " "
echo "====== Done ======"
# echo " "
# echo " "
# echo "====== Downloading Cloud Point Data ======"
# mkdir -p "./CloudPoint/QL2_3DEP_LiDAR_IN_2017_2019_laz"
# wget --no-parent -e robots=off -U mozilla -A '*.laz' -r -nd -N -P "./CloudPoint/QL2_3DEP_LiDAR_IN_2017_2019_laz" "https://lidar.jinha.org/QL2_3DEP_LiDAR_IN_2017_2019_laz/"
# echo " "
# echo "====== Done ======"