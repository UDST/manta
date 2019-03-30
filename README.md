# traffic-simulator-extended
Traffic Microsimulator with Extensions

Traffic flow is modeled using the
[Intelligent Driver Model (IDM)](https://en.wikipedia.org/wiki/Intelligent_driver_model)
.

## Dependencies

 - Boost 1.59
 - OpenCV (used versions: 2.4.12 in Windows; 3.2.0 in Ubuntu)
 - CUDA (used versions: 8.0.61 in Windows; 9.0 in Ubuntu)
 - g++ (used versions: 6.4.0 in Ubuntu)
 - Qt5 (used versions: 5.9.5 in Ubuntu)
 - qmake (used versions: 3.1 in Ubuntu)

## Installation & Compilation

Once the necessaries dependencies are installed you can use the following lines to make sure the
correct versions of each one are used:
```bash
export PATH=PATH=/usr/local/cuda-9.0/bin:$PATH
export LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LIBRARY_PATH 
export LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LD_LIBRARY_PATH 
```

In Linux you can also add those `export`statements lines at the end of your user's `~/.bashrc` to
avoid re-entering them on each session.

Clone in your home directory with:
```bash
git clone git@github.com:ual/traffic-simulator-extended.git ~/traffic-simulator-extended && cd ~/traffic-simulator-extended
```

If necessary checkout a different branch than master (`maintenance` for instance):
```bash
git checkout maintenance
```

Create `Makefile` and compile with:
```bash
sudo qmake LivingCity/LivingCity.pro && sudo make -j
```

## Running

If you wish to edit configuration, modify `command_line_options.ini`.

Run with:
```bash
cd LivingCity
./LivingCity
```

## Functions called when running
1) Once ./LivingCity is invoked without the GUI, LC_main() is called. 
2) Within LC_main(), runB18Simulation() is called, which begins the traffic microsimulation.
	- runB18Simulation(), located in traffic/b18CommandLineVersion.cpp.
		- loadB2018RoadGraph() loads the road graph from Geoff's street network file. specifically, it loads from the csv files that contain the full bay area street network (bay_area_full_strongly_nodes.csv), the edges (full_edges_speed_capacity.csv), and the origin destination pairs (od.csv). Each file is then parsed to create a network based on lat-long.
		- initSimulator() initializes the newly loaded street network.
		- createB2018People() generates a synthetic population with certain attributes. specifically, it calls resetTrafficPersonJob() and loadB18TrafficPeople().
			- resetTrafficPersonJob() makes every active person in the synthetic population to be 0 (resets each person)
			- loadB18TrafficPeople() creates a random seed from boost using normal_distribution() with mean 7.5 and variance .75. Then, using that random seed, it calls randomPerson().
				- randomPerson(), which takes as arguments the index of the person, the person pointer, the origin lat/long, destination lat/long, and start time. It then sets the person's acceleration, braking, time headway, the number of steps, the number of lanes in the edge, and the state of its lane changes. This is then called by an overarching randomPerson() function that finds the closest nodes to the person based on their house position and job position. all of these persons are now in a persons vector (trafficPersonVec)
		- simulateInCPU_MultiPass() if using CPU or simulateInGPU() if using GPU (we do use the GPU)
			- createLaneMap() is called, which creates edge data and finds the required sizes, checks the distribution of the street lengths and widths, generates intersection info (which determines the number of in/out edges for the vertex. now every attribute about the edge is stored in the edgesData vector. it also creates an intersection vector, which contains the state of the intersection, the number of in/out edges of the intersection, and what the next 'event' is. it also determines the edge angle from the vertex itself. each intersection also has "edge," which describes all the edges to/from that intersection.
			- generateRoutes() actually does traffic assignment using Johnson's shortest path method. Fede is suggesting using a different shortest path algorithm here. This takes in the road graph, the trafficPersonVec, indexPathVec, lane map mapping, weight, and peoplePathSampling.
			- generateRoutesMulti() uses Djikstra's algorithm, where it creates a counter to keep track where to put more edges, and then you update the weight edges (which could be travel time or cost), use the average speed of the edge sampled from a former step. 
				- calculateSeveralPeopleRoute() sees whether the person is at their destination; if not, it uses djikstra to figure out each person's path to their destination. it seems like indexPathVec gives the indices of each edge that is traversed by a person.
			- simulateInCPU()
				- simulateOneIntersectionCPU() takes the current time, the intersections vector pointer, and the traffic lights pointer. it first checks if current time is after the next event for the intersection and if the intersection has real edges. if so, it writes the old traffic lights to 0 in memory. then it loops through the intersection's edges and changes the state of the intersection to the next light, and then it sets that edge's characteristics (such as numLines), and then it sets the traffic light to green, and finally sets the intersection's nextEvent to the current time plus the delta in time to the next event.
				- simulateOnePersonCPU() takes the delta time, current time, map for read and write shifting, trafficPersonVec, indexPathVec, edgesData, laneMap, intersections, and trafficlights. it first checks if person should wait (if time_departure > currentTime) or start (time_departure < currentTime). In the latter, it then finds the first edge and then updates that person's edge data (i.e. creates the trafficPersonVec with the specific edges for that particular person's route). It then also creates the lane map for that person's first edge, determining the exact lanes the person will pass on that edge (it does this by first counting all the empty cells on the edge to see if there is space for the car). then, it checks whether the car is ready for the next edge, so it will change it so that it is in the intersection. now the person/car starts moving, trying to go from the current edge to the next edge. to do so, first it finds the car in front of it in the same lane. the next thing is reaching the traffic light, and it first sees if the car is allowed to do move before and at the intersection. then, it moves into the new edge after the traffic light by updating the current edge and the lane in the edge. now that it's in the new lane, it starts checking for the speed of the car in front. if it's slower than us, you make a discretionary lane change. in this entire process, it is also calculating emissions and gas consumption. if the route requires a mandatory lane change, it checks to see if there are any open cells and the speeds of the adjacent cars, and then waits x amount of time to make the lane change (this is part of calculateGapsLC(). this process repeats along the full route assigned to that person.
				- sampleTraffic() counts the number of people in the trafficPersonVec, then it loops through the number of people, and if each person's trafficPersonVec is active, it updates the accSpeedPerLinePerTimeInterval for each person's car with "offset" as the increment, which is the numLines * samplingNumber.

Overall, when the simulator is called, it first loads a street network from file, which has edges and nodes. After this, the actual simulator is run. Since the goal is to simulate a synthetic population's movement in traffic, the simulator first creates a synthetic population, where each person is an element in the vector trafficPersonVec. It creates this population randomly, and it determines each person's origin and destination, their start time, and then also microscopic details like their accelerations, number of lane changes, etc.

Once we have this information about each person, it starts simulating the actual route and the actual edges it will take by using either Johnson's path algorithm or Djikstra's shortest path algorithm. Once the edges are known, then it begins more of the real-time movement. Specifically, it will create a lanemap of each edge that shows which cells are taken and which ones are empty. Based on this and discretionary/mandatory lane changing logic, the person's position will move, given that we know what its next edge needs to be. Then, it will reach an intersection (node), where there is a traffic light (or some level of traffic control). Based on this traffic control, it determines whether it can move at the intersection (based again on cells being populated as a proxy for traffic.) If it notices it can move onto the new edge, the new edge now becomes the current edge. 

This is run in a loop across all trafficPersonVecs (people) and dynamically changes according to various person positions.
