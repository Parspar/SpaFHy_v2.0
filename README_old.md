# SpaFHy v2

*SpaFHy version for gridded simulation of forests and peatlands*

 - Canopy and bucket mostly as in SpaFHy v1
 - Option for physically-based 2D groundwater flow and groundwater storage
 - Option for TOPMODEL as in SpaFHy v1

### Example for running model and plotting some results
Data for example simulation in folder testcase_inputs (1000 nodes)
```
from model_driver import driver
from iotools import read_results
import matplotlib.pyplot as plt

# runs model
outputfile = driver(create_ncf=True)

# reads results from .nc-file
results = read_results(outputfile)

# plots ground water level for first ten nodes
plt.figure()
results['soil_ground_water_level'][:,0,:10].plot.line(x='date')


