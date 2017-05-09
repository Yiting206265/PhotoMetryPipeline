#!/usr/bin/python
# -*- coding: utf-8 -*-
# File: gkastro.py
# Created: 2016-03-31 by gks 
"""
Description: Helpful for astronomy
"""

from __future__ import print_function

import everest
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import math

from astropy.visualization import LogStretch, SqrtStretch, AsinhStretch, HistEqStretch,ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize

norm_mean_sub = lambda x: x - np.nanmean(x)
norm_mean     = lambda x: x/np.nanmean(x)
norm_median   = lambda x: x/np.median(x)
compactString = lambda string: string.replace(' ', '').replace('-', '').lower()
cosd = lambda x : np.cos(np.deg2rad(x))
sind = lambda x : np.sin(np.deg2rad(x))

def stretch_data(data,method="HistEqStretch"):
    """
    methods = 
    LogStretch,
    SqrtStretch,
    AsinhStretch,
    HistEqStretch
    """
    if method=="LogStretch":
        norm = ImageNormalize(stretch=LogStretch(data))
    elif method=="SqrtStretch":
        norm = ImageNormalize(stretch=SqrtStretch(data))
    elif method=="AsinhStretch":
        norm = ImageNormalize(stretch=AsinhStretch(data))
    elif method=="HistEqStretch":
        norm = ImageNormalize(stretch=HistEqStretch(data))
    else:
        norm = data
    return norm

