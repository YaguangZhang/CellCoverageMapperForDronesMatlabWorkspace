import asyncio
from directory_downloader import DDownloader


async def main():
    url = r"https://ftp.mozilla.org/pub/artwork/"
    downloader = DDownloader(url)
    await downloader.fetch_file_links(extensions=['.zip'])  # returns set of downloadable file urls
    await downloader.download_files()  # download all files to current directory


if __name__ == '__main__':
    asyncio.run(main())