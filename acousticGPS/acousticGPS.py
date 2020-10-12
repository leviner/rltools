"""
GPS from raw files

For details on use and alternative, see acousticGPS.ipynb

## pyEcholab

Using the [`pyEcholab`](https://github.com/CI-CMG/pyEcholab.git) toolkit developed by NOAA for reading echosounder data, to pull nmea data from EK60 raw files, you can install directly from the repository:

```python
!pip install git+https://github.com/CI-CMG/pyEcholab.git
```

This is setup to use pyEcholab rather than [`echoPype`](https://echopype.readthedocs.io/en/latest/index.html) which uses an intermediate step of forming the netCDF structure. Pulling the lat/lon with `echoPype` will be included as an example at the bottom of the notebook at some point.

The downside of the setup below is that it doesn't write out the file until it's done, but that is the trade off for maintaining iterations intervals of < 2s. File building within the loop exponentially increases iteration time due to loading the document and essentially having to read the dataframe twice each step.

"""

from glob import glob
from echolab2.instruments import EK60
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import dask

def acousticGPS(rawFileDir,save=0, outFile='out.csv') # Input is raw file path, outFile is csv to be saved to if save = 1
    rawfiles = sorted(glob(rawFileDir+'*.raw'))
    df_list= []
    for i in trange(2000,len(rawfiles)):
        try:
            ek60 = EK60.EK60()
            ek60.read_raw(rawfiles[i])
            rawd = ek60.get_raw_data(channel_number=1)
            sv = rawd.get_Sv()
            positions = ek60.nmea_data.interpolate(sv, 'position')
            dfGPS = dask.delayed(pd.DataFrame)(positions)
            df_list.append(dfGPS)
        except:
            print(rawfiles[i])
    df = dask.delayed(pd.concat)(df_list).compute()
    if save == 1:
        print('Saving to '+outFile)
        df.to_csv(outFile)
    return df