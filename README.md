# manta
Microsimulation Analysis for Network Traffic Assignment
## Dependencies

 - Boost 1.59
 - OpenCV (used versions: 3.2.0 in Ubuntu)
 - CUDA (used versions: 9.0 in Ubuntu)
 - g++ (used versions: 6.4.0 in Ubuntu)
 - Qt5 (used versions: 5.9.5 in Ubuntu)
 - qmake (used verions: 3.1 in Ubuntu)

## Installation & Compilation

Once the necessary dependencies are installed, you can use the following lines to make sure the
correct versions of each one are used:
```bash
export PATH=PATH=/usr/local/cuda-9.0/bin:$PATH
export LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LIBRARY_PATH 
export LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LD_LIBRARY_PATH 
```

You can also add the `export` lines at the end of your user's `~/.bashrc` to
avoid re-entering them in each session.

Clone the repo in your home directory with:
```bash
git clone git@github.com:udst/manta.git ~/manta && cd ~/manta
```

If necessary, you can checkout a different branch than master (`edge_speeds_over_time` for instance):
```bash
git checkout edge_speeds_over_time
```

Create `Makefile` and compile with:
```bash
sudo qmake LivingCity/LivingCity.pro && sudo make -j
```

## Running

If you wish to edit the microsimulation configuration, modify `command_line_options.ini`.

Run with:
```bash
cd LivingCity
./LivingCity
```
In `command_line_options.ini`, you can modify the timestep (default dt = .5 seconds), simulation time range (default 5am-12pm), specific path to the network and OD demand, whether you would like to compute the shortest paths or use a previous version, and which type of routing.

