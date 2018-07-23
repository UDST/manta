# traffic-simulator-extended
Traffic Microsimulator with Extensions

## Dependencies

 - Boost 1.59
 - OpenCV (used versions: 2.4.12 in Windows; 3.2.0 in Ubuntu)
 - CUDA (used versions: 8.0.61 in Windows; 9.0 in Ubuntu)
 - g++ (used verions: 6.4.0 in Ubuntu)
 - Qt5 (used verions: 5.9.5 in Ubuntu)
 - qmake (used verions: 3.1 in Ubuntu)

Once the necessaries dependencies are installed you can use the following lines to make sure the
correct verions of each one are used:
```bash
export PATH=PATH=/usr/local/cuda-9.0/bin:$PATH
export LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LIBRARY_PATH 
```
Note that I used `LIBRARY_PATH` instead of `LD_LIBRARY_PATH` because the last one was not working 
in the server with my user (view this 
[SO thread](https://stackoverflow.com/questions/13292261/ld-library-path-doesnt-seem-to-work)).

In Linux you can also add those lines at the end of your user's `~/.bashrc` to avoid re-entering
them on each session.

## Installation

Clone in your home directory with:
```bash
git clone git@github.com:ual/traffic-simulator-extended.git ~/traffic-simulator-extended && cd traffic-simulator-extended
```
Create `Makefile` and compile with:
```bash
sudo qmake LivingCity/LivingCity.pro && sudo make -j
```

If you wish to edit configuration, modify `command_line_options.ini`.

## Running

Run with:
```bash
cd LivingCity
./LivingCity
```

