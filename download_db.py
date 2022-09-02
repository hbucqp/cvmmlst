#!/usr/bin/python3

import orjson
import urllib.request
import pandas as pd
import requests

# request = "https://rest.pubmlst.org/db"

# with urllib.request.urlopen(request) as response:
#     response_content = response.read()


# data = orjson.loads(response_content)

# for resource in data:
#     if resource['databases']:
#         for db in resource['databases']:
#             print(str(db) + '\n')

r = requests.get(
    'http://rest.pubmlst.org/db/pubmlst_blicheniformis_seqdef/schemes/14').json()
print(r)
