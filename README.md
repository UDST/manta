# manta

Microsimulation Analysis for Network Traffic Assignment. MANTA employs a highly parallelized GPU implementation that is fast enough to run simulations on large-scale demand and networks within a few minutes - metropolitan and regional scale with hundreds of thousands of nodes and edges and millions of trips. We test our platform to simulate the entire Bay Area metropolitan region over the course of the morning using half-second time steps. The runtime for the nine-county Bay Area simulation is just over four minutes, not including routing and initialization. This computational performance significantly improves state of the art in large-scale traffic microsimulation and offers new capacity for analyzing the detailed travel patterns and travel choices of individuals for infrastructure planning and emergency management.

![](https://github.com/UDST/manta/blob/main/bay_bridge_trips.png)

## Dependencies

 - Boost 1.59
 - OpenCV (used versions: 3.2.0 in Ubuntu)
 - CUDA (used versions: 9.0 in Ubuntu)
 - g++ (used versions: 6.4.0 in Ubuntu)
 - Qt5 (used versions: 5.9.5 in Ubuntu)
 - qmake (used versions: 3.1 in Ubuntu)

## Installation & Compilation

Once the necessary dependencies are installed, you can use the following lines to make sure the
correct versions of each one are used:
```bash
export PATH=PATH=/usr/local/cuda-9.0/bin:$PATH
export LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LIBRARY_PATH 
export LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/{YOUR_USERNAME}/pandana/src
```

You can also add the `export` lines at the end of your user's `~/.bashrc` to
avoid re-entering them in each session.

Clone the repo in your home directory with:
```bash
git clone git@github.com:udst/manta.git ~/manta && cd ~/manta
```

Clone the [Pandana repository](https://github.com/UDST/pandana) to your home directory. This is necessary since MANTA now uses a fast contraction hierarchies framework for shortest path routing. Previously implemented shortest path frameworks include Johnson's all pairs shortest path and a parallelized Dijkstra's priority queue.

If necessary, you can checkout a different branch than master (`edge_speeds_over_time` for instance):
```bash
git checkout edge_speeds_over_time
```

Create `Makefile` and compile with:
```bash
sudo qmake LivingCity/LivingCity.pro && sudo make -j
```

## Running

If you wish to edit the microsimulation configuration, modify `command_line_options.ini`, which contains the following:

```[General]
GUI=false
USE_CPU=false
NETWORK_PATH=berkeley_2018/new_full_network/
USE_JOHNSON_ROUTING=false
USE_SP_ROUTING=true
USE_PREV_PATHS=true
LIMIT_NUM_PEOPLE=256000
ADD_RANDOM_PEOPLE=false
NUM_PASSES=1
TIME_STEP=0.5
START_HR=5
END_HR=12
```

Here, you can modify the:

1. `GUI` - deprecated. Do not touch.
2. `USE_CPU` - deprecated. Do not touch.
3. `NETWORK_PATH` - specific path to the network files. Default is `berkeley_2018/new_full_network/`.
4. `USE_JOHNSON_ROUTING` - uses Johnson's all pairs shortest path routing. This should always be set to `false`.
5. `USE_SP_ROUTING` - uses new SP routing framework. This should always be set to `true`.
6. `USE_PREV_PATHS` - uses paths already produced and saved to file. Set to `false` if running for the first time. Set to `true` if the simulation was already run and it was saved to file. 
7. `LIMIT_NUM_PEOPLE` - deprecated. Do not touch.
8. `ADD_RANDOM_PEOPLE` - deprecated. Do not touch.
9. `NUM_PASSES` - the number of times the simulation is run. Set to 1.
10. `TIME_STEP` - timestep. Default is .5 seconds.
11. `START_HR` - start hour of the simulation. Default is 5am.
12. `END_HR` - end hour of the simulation. Default is 12pm.

Run with:
```bash
cd LivingCity
./LivingCity
```

## Development

Should you wish to make any changes, please create a new branch. In addition, once the original Makefile is created, you can simply run `sudo make -j` from the `manta/` directory to compile any new changes.

## Acknowledgments

This repository and code have been developed and maintained by Pavan Yedavalli, Ignacio Garcia Dorado, Krishna Kumar, and Paul Waddell. This work heavily derives from Ignacio Garcia Dorado's [Automatic Urban Modeling project](http://www.ignaciogarciadorado.com/p/2014_EG/2014_EG.html).

If this code is used in any shape or form for your project, please cite this paper accordingly:

P. Yedavalli, K. Kumar, and P. Waddell, “Microsimulation Analysis for Network Traffic Assignment (MANTA) at Metropolitan-Scale for Agile Transportation Planning,” arXiv:2007.03614 [physics], Jul. 2020, Available: http://arxiv.org/abs/2007.03614.

Thank you!




