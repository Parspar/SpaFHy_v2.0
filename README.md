# SpaFHy-Peat
*SpaFHy version for gridded simulation of drained peatland forests.*

 - Canopy and mosslayer as in SpaFHy
 - Soil water storage based on equilibrium state and Hooghoudt's drainage equation.

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
```
**Edits:**
5-Jul-2019: Modified interpolation function calls to speed up computation
26-Aug-2019: Rootzone moisture restriction on transpiration from sompa branch
28-Aug-2019: Reading forcing files adjusted
