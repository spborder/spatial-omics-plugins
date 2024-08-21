"""

Plugin implementation to perform selective aggregation of properties to other intersecting structures

"""

import os
import sys

import girder_client
import json

import requests
import geopandas as gpd
from shapely.geometry import Polygon
import pandas as pd

from io import BytesIO
from typing_extensions import Union

from tqdm import tqdm


















































