#!/bin/bash
# Remove the write permission for all users in the LiDAR folder.

find /home/coverage/CellCoverageMapper/Lidar -type f -exec chmod a-w '{}' \;