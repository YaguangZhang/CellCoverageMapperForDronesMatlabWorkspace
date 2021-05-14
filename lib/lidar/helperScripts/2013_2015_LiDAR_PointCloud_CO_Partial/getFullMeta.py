# Download the meta files for the CO LiDAR data.
#
# Developed in Python 3.
#
# Yaguang Zhang, 2021/05/12

import pip
import os

# We will use a modified local copy of directory_downloader
from directory_downloader import DDownloader

def safeImport(packageStr):
    try:
        return __import__(packageStr)
    except ImportError:
        pip.main(['install', packageStr])
        return __import__(packageStr)

asyncio = safeImport('asyncio')
safeImport('lxml')

def main():
    numOfWorkers = 10

    curDirPath = os.path.dirname(os.path.realpath(__file__))
    pathToSaveMetaFiles = os.path.join(curDirPath, '..', 'meta')
    if not os.path.exists(pathToSaveMetaFiles):
        os.makedirs(pathToSaveMetaFiles)

    urlToMetaFiles = r"https://thor-f5.er.usgs.gov/ngtoc/metadata/waf/elevation/lidar_point_cloud/las/USGS_LPC_CO_SoPlatteRiver_Lot2a_2013_LAS_2015/"
    asyncio.run(paraDownload(urlToMetaFiles, pathToSaveMetaFiles,
                ['xml'], numOfWorkers))

async def paraDownload(url, destDir, extStr, numOfWorker=5):

    downloader = DDownloader(url, workers=numOfWorker, directory=destDir)
    await downloader.fetch_file_links(extensions=extStr)
    await downloader.download_files(flagNoDir=True)

if __name__ == '__main__':
    main()