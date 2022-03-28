# For fetching URLs with Python 3.
#
# Yaguang Zhang, Purdue, 09/25/2019

import threading
import urllib.request
import ssl
from socket import timeout
from urllib.error import HTTPError, URLError
import logging

# Ignore SSL certificate errors
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def fetch_url(url, htmls, idx):
    TIMEOUT_IN_S = 10

    try:
        htmls[idx] = urllib.request.urlopen(url,
            timeout=TIMEOUT_IN_S, context=ctx).read()
    except (HTTPError, URLError) as error:
        logging.warning('Data not retrieved because %s\nURL: %s', error, url)
    except timeout:
        logging.warning('Socket timed out - URL %s', url)
    else:
        logging.info('Access successful.')

def concurrent_fetch_urls(urls):
    numOfUrls = len(urls)

    threads = [None] * numOfUrls
    htmls = [None] * numOfUrls

    for idx in range(numOfUrls):
        threads[idx] = threading.Thread(target=fetch_url, \
            args=(urls[idx], htmls, idx))
        threads[idx].start()

    for thread in threads:
        thread.join()

    return htmls

if __name__ == '__main__':
    # For testing.
    urls = ["https://www.google.com", "http://www.yahoo.com", "http://www.baidu.com", \
        "https://nationalmap.gov/epqs/pqs.php?x=-80&y=42&units=Meters&output=json"]
    htmls = concurrent_fetch_urls(urls)
    print('Successfully fetched', sum(h is not None for h in htmls), 'out of', len(urls), 'URLs!')
    for i, val in enumerate(htmls):
        if val == None: print ('    Failed URL: ', urls[i])