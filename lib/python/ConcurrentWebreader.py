# For fetching URLs with Python 3.
#
# Yaguang Zhang, Purdue, 09/25/2019

import threading
import urllib.request
import queue
import ssl

# Ignore SSL certificate errors
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def fetch_url(url):
    return urllib.request.urlopen(url, context=ctx).read()

def concurrent_fetch_urls(urls):
    numOfUrls = len(urls)

    threads = [None] * numOfUrls
    htmls = [None] * numOfUrls

    results = queue.Queue(numOfUrls)

    for idx in range(numOfUrls):
        threads[idx] = threading.Thread(target=lambda q, fctArgs: q.put((idx, fetch_url(fctArgs))), \
            args=(results, urls[idx]))
        threads[idx].start()

    for thread in threads:
        thread.join()

    for idx in range(results.qsize()):
        (idx, html) = results.get()
        htmls[idx] = html 

    return htmls

if __name__ == '__main__':
    # For testing.
    urls = ["https://www.google.com", "http://www.yahoo.com", "http://www.baidu.com", \
        "https://nationalmap.gov/epqs/pqs.php?x=-80&y=42&units=Meters&output=json"]
    htmls = concurrent_fetch_urls(urls)
    print('Successfully fetched', sum(h is not None for h in htmls), 'out of', len(urls), 'URLs!')
