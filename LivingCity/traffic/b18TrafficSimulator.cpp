#include "B18TrafficSimulator.h"

#include "../global.h"
#include "../LC_GLWidget3D.h"
#include "../LC_UrbanMain.h"
#include <thread>

#include "b18TrafficJohnson.h"
#include "b18CUDA_trafficSimulator.h"

#define DEBUG_TRAFFIC 0
#define DEBUG_SIMULATOR 0
#define DEBUG_T_LIGHT 0

#ifdef __linux__ 
#include <unistd.h>
void printPercentageMemoryUsed() {
  // TODO
}
#elif _WIN32
#include <windows.h>
void printPercentageMemoryUsed() {
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  printf("** Memory usage %ld%% --> Free %I64d MB\n", status.dwMemoryLoad, status.ullAvailPhys / (1024*1024));
}
#endif
namespace LC {

//LCUrbanMain *clientMain;

const float s_0 = 1.5f * 4.12f; // ALWAYS USED
const float intersectionClearance = 7.8f;
const bool calculatePollution = false;

B18TrafficSimulator::B18TrafficSimulator() {
  initialized = false;
  simRoadGraph = 0;
  threadNumber = 0;
}//

B18TrafficSimulator::~B18TrafficSimulator() {
}//

void B18TrafficSimulator::initSimulator(
  float _deltaTime,
  RoadGraph *originalRoadGraph,
  LCUrbanMain *urbanMain) {
  deltaTime = _deltaTime;

  if (simRoadGraph == 0) { //!!!!! NOT COPY IF IT IS THERE
    simRoadGraph = new RoadGraph(*originalRoadGraph);
  }

  clientMain = urbanMain;
  initialized = true;
}

void B18TrafficSimulator::createRandomPeople(float startTime, float endTime,
    int numberPeople,
    PeopleJobInfoLayers &peopleJobInfoLayers) {
  if (initialized == false) {
    printf("Error: initSimulator was not called\n");
    return;
  }

  b18TrafficOD.resetTrafficPersonJob(trafficPersonVec);
  b18TrafficOD.createRandomPeople(numberPeople, startTime, endTime,
                                  trafficPersonVec,
                                  peopleJobInfoLayers, simRoadGraph->myRoadGraph_BI);
}//

void B18TrafficSimulator::createB2018People(float startTime, float endTime, int limitNumPeople) {
  if (initialized == false) {
    printf("Error: initSimulator was not called\n");
    return;
  }

  b18TrafficOD.resetTrafficPersonJob(trafficPersonVec);
  b18TrafficOD.loadB18TrafficPeople(startTime, endTime, trafficPersonVec, simRoadGraph->myRoadGraph_BI, limitNumPeople);

}


void B18TrafficSimulator::resetPeopleJobANDintersections() {
  if (initialized == false) {
    printf("Error: initSimulator was not called\n");
    return;
  }

  //printf("resetPeopleJobANDintersections\n");
  b18TrafficOD.resetTrafficPersonJob(trafficPersonVec);
  b18TrafficLaneMap.resetIntersections(intersections, trafficLights);
}//


void B18TrafficSimulator::createLaneMap() { //
  if (initialized == false) {
    printf("Error: initSimulator was not called\n");
    return;
  }

  b18TrafficLaneMap.createLaneMap(*simRoadGraph, laneMap, edgesData,
                                   intersections, trafficLights, laneMapNumToEdgeDesc, edgeDescToLaneMapNum);
}//

void B18TrafficSimulator::generateCarPaths(bool useJohnsonRouting) { //
  if (initialized == false) {
    printf("Error (generateCarPaths): initSimulator was not called\n");
    return;
  }
  if (useJohnsonRouting) {
    printf("***generateCarPaths Start generateRoute Johnson\n");
    B18TrafficJohnson::generateRoutes(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, 0);
  } else {
    printf("***generateCarPaths Start generateRoutesMulti Disktra\n");
    cudaTrafficPersonShortestPath.generateRoutesMulti(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, 0);
  }

}//

void B18TrafficSimulator::calculateAndDisplayTrafficDensity() {
  if (initialized == false) {
    printf("Error: initSimulator was not called\n");
    return;
  }

  calculateAndDisplayTrafficDensity(accSpeedPerLinePerTimeInterval,
                                    numVehPerLinePerTimeInterval, edgeDescToLaneMapNum, laneMapNumToEdgeDesc,
                                    trafficLights.size());//num lanes
}//




//////////////////////////////////////////////////
// GPU
//////////////////////////////////////////////////
void B18TrafficSimulator::simulateInGPU(float startTimeH, float endTimeH, bool useJohnsonRouting) {

  //////////////////////////////////////////
  printf("Simulate GPU  createLaneMap();");
  createLaneMap();

  ////////////////////////////
  printf("Simulate GPU  routes();");
  QTime pathTimer;
  pathTimer.start();

  int weigthMode = 0;
  int nP = 0;
  float peoplePathSampling[] = {1.0f, 1.0f, 0.5f, 0.25f, 0.12f, 0.67f};

  if (useJohnsonRouting) {
    printf("***Start generateRoute Johnson\n");
    B18TrafficJohnson::generateRoutes(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, weigthMode, peoplePathSampling[nP]);
  } else {
    printf("***Start generateRoutesMulti Disktra\n");
    cudaTrafficPersonShortestPath.generateRoutesMulti(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, weigthMode, peoplePathSampling[nP]);
  }
  printf("***Path time %f\n", pathTimer.elapsed());

  /////////////////////////////////////
	printf("  >>Start Simulation\n");

  float startTime = startTimeH*3600.0f;//7.0f
  float endTime = endTimeH*3600.0f;//8.0f//10.0f

	float currentTime=23.99f*3600.0f;

	for(int p=0;p<trafficPersonVec.size() ;p++){
		if(currentTime>trafficPersonVec[p].time_departure){
			currentTime=trafficPersonVec[p].time_departure;
		}
	}
	int numInt=currentTime/deltaTime;//floor
	currentTime=numInt*deltaTime;
	uint steps=0;
	steps=(currentTime-startTime)/deltaTime;
  // start as early as starting time
  if (currentTime < startTime) { //as early as the starting time
    currentTime = startTime;
  }
	QTime timer;
	G::global()["cuda_render_displaylist_staticRoadsBuildings"]=1;//display list
	timer.start();
	// 1. Init Cuda
  b18InitCUDA(trafficPersonVec, edgesData, laneMap, trafficLights, intersections);
	// 2. Excute
	printf("First time_departure %f\n",currentTime);
	while(currentTime<endTime){//24.0f){
    b18SimulateTrafficCUDA(currentTime, trafficPersonVec.size(), intersections.size());

    if (clientMain->ui.b18RenderSimulationCheckBox->isChecked()) {
			steps++;
      if (steps%clientMain->ui.b18RenderStepSpinBox->value() == 0) {//each "steps" steps, render
        b18GetDataCUDA(trafficPersonVec, trafficLights);
				QString timeT;

				int timeH=currentTime/3600.0f;
				int timeM=(currentTime-timeH*3600.0f)/60.0f;
				timeT.sprintf("%d:%02d",timeH,timeM);
        clientMain->ui.b18TimeLCD->display(timeT);

        clientMain->glWidget3D->updateGL();
				QApplication::processEvents();
			}
		}
		currentTime+=deltaTime;
	}
	// 3. Finish
  b18GetDataCUDA(trafficPersonVec, trafficLights);
	/////
	uint totalNumSteps=0;
		float totalGas=0;
		for(int p=1;p<trafficPersonVec.size();p++){
			totalNumSteps+=trafficPersonVec[p].num_steps;
			totalGas+=trafficPersonVec[p].gas;
		}
		avgTravelTime=(totalNumSteps)/(trafficPersonVec.size()*60.0f);//in min
		printf("Total num steps %u Avg %f min Avg Gas %f. Calculated in %d ms\n",totalNumSteps,avgTravelTime,totalGas/trafficPersonVec.size(),timer.elapsed());

	///
  b18FinishCUDA();
	printf("  <<End Simulation TIME: %d ms.\n",timer.elapsed());
  G::global()["cuda_render_displaylist_staticRoadsBuildings"] = 3;//kill display list
  clientMain->glWidget3D->updateGL();
	printf(" <<simulate\n");
}//


//////////////////////////////////////////////////
// CPU
//////////////////////////////////////////////////

// calculate if the car will turn left center or right (or the same)
void calculateLaneCarShouldBe(
  uint curEdgeLane,
  uint nextEdge,
  std::vector<B18IntersectionData> &intersections,
  uint edgeNextInters,
  ushort edgeNumLanes,
  ushort &initOKLanes,
  ushort &endOKLanes) {
  initOKLanes = 0;
  endOKLanes = edgeNumLanes;

  if (DEBUG_TRAFFIC == 1) {
    printf("curEdgeLane %05x nextEdge %05x\n", curEdgeLane, nextEdge);
  }

  if (DEBUG_TRAFFIC == 1) {
    for (int eN = 0; eN < intersections[edgeNextInters].totalInOutEdges; eN++) {
      printf("* procEdge %05x\n", intersections[edgeNextInters].edge[eN]);
    }
  }

  bool currentEdgeFound = false;
  bool exitFound = false;
  ushort numExitToTake = 0;
  ushort numExists = 0;

  for (int eN = intersections[edgeNextInters].totalInOutEdges - 1; eN >= 0; eN--) {  // clockwise
    uint procEdge = intersections[edgeNextInters].edge[eN];

    if ((procEdge & kMaskLaneMap) == curEdgeLane) { //current edge 0xFFFFF
      if (DEBUG_TRAFFIC == 1) {
        printf("CE procEdge %05x\n", procEdge);
      }

      currentEdgeFound = true;

      if (exitFound == false) {
        numExitToTake = 0;
      }

      continue;
    }


    if ((procEdge & kMaskInEdge) == 0x0) { //out edge 0x800000
      if (DEBUG_TRAFFIC == 1) {
        printf("   procEdge %05x\n", procEdge);
      }

      numExists++;

      if (currentEdgeFound == true) {
        numExitToTake++;
      }

      if (currentEdgeFound == false && exitFound == false) {
        numExitToTake++;
      }
    }

    if ((procEdge & kMaskInEdge) == nextEdge) {
      exitFound = true;
      currentEdgeFound = false;

      if (DEBUG_TRAFFIC == 1) {
        printf("NE procEdge %05x\n", procEdge);
      }
    }
  }

  if (DEBUG_TRAFFIC == 1) {
    printf("Num extis %u Num exit to take %u%\n", numExists, numExitToTake);
  }

  if (edgeNumLanes == 0) {
    printf("ERRRROR\n");
  }

  switch (edgeNumLanes) {
  /// ONE LANE
  case 1:
    initOKLanes = 0;
    endOKLanes = 1;
    break;

  /// TWO LANE
  case 2:
    switch (numExists) {
    case 1:
    case 2://all okay
      initOKLanes = 0;
      endOKLanes = 2;
      break;

    case 3:
      if (numExitToTake > 2) { //left
        initOKLanes = 0;
        endOKLanes = 1;
        break;
      }

      initOKLanes = 1;
      endOKLanes = 2;
      break;

    default:
      if (DEBUG_TRAFFIC == 1) {
        printf("2 defalt\n");
      }

      if (numExitToTake >= numExists - 1) {
        initOKLanes = 0;
        endOKLanes = 1;
        break;
      }

      initOKLanes = 1;
      endOKLanes = 2;
      break;
    }

    break;

  /// THREE LANE
  case 3:
    switch (numExists) {
    case 1:
    case 2://all okay
      initOKLanes = 0;
      endOKLanes = 3;
      break;

    case 3:
      if (numExitToTake > 2) { //left
        initOKLanes = 0;
        endOKLanes = 1;
        break;
      }

      initOKLanes = 1;
      endOKLanes = 3;
      break;

    default:
      if (numExitToTake >= numExists - 1) {
        initOKLanes = 0;
        endOKLanes = 1;
        break;
      }

      initOKLanes = 1;
      endOKLanes = 2;
      break;
    }

    break;

  case 4:
    switch (numExists) {
    case 1:
    case 2://all okay
      initOKLanes = 0;
      endOKLanes = 4;
      break;

    case 3:
      if (numExitToTake == 1) { //right
        initOKLanes = 3;
        endOKLanes = 4;
      }

      if (numExitToTake > 3) { //left
        initOKLanes = 0;
        endOKLanes = 1;
        break;
      }

      initOKLanes = 1;
      endOKLanes = 4;
      break;

    default:
      if (numExitToTake == 1) { //right
        initOKLanes = edgeNumLanes - 1;
        endOKLanes = edgeNumLanes;
      }

      if (numExitToTake >= numExists - 2) {
        initOKLanes = 0;
        endOKLanes = 2;
        break;
      }

      initOKLanes = 1; //also lane 2
      endOKLanes = edgeNumLanes;
    }

    break;

  default:
    switch (numExists) {
    case 1:
    case 2://all okay
      initOKLanes = 0;
      endOKLanes = edgeNumLanes;
      break;

    case 3:
      if (numExitToTake == 1) { //right
        initOKLanes = edgeNumLanes - 1;
        endOKLanes = edgeNumLanes;
      }

      if (numExitToTake > edgeNumLanes - 2) { //left
        initOKLanes = 0;
        endOKLanes = 2;
        break;
      }

      initOKLanes = 1;
      endOKLanes = edgeNumLanes;
      break;

    default:
      if (numExitToTake < 2) { //right
        initOKLanes = edgeNumLanes - 2;
        endOKLanes = edgeNumLanes;
      }

      if (numExitToTake >= numExists - 2) {
        initOKLanes = 0;
        endOKLanes = 2;
        break;
      }

      initOKLanes = 1; //also lane 2
      endOKLanes = edgeNumLanes - 1;
    }

    break;
  }
}//

void calculateGapsLC(
  uint mapToReadShift,
  std::vector<uchar> &laneMap,
  uchar trafficLightState,
  uint laneToCheck,
  ushort numLinesEdge,
  float posInMToCheck,
  float length,
  uchar &v_a,
  uchar &v_b,
  float &gap_a,
  float &gap_b) {
  ushort numOfCells = ceil(length);
  ushort initShift = ceil(posInMToCheck);
  uchar laneChar;
  bool found = false;

  // CHECK FORWARD
  //printf("initShift %u numOfCells %u\n",initShift,numOfCells);
  for (ushort b = initShift - 1; (b < numOfCells) &&  (found == false); b++) { //NOTE -1 to make sure there is none in at the same level
    // laneChar = laneMap[mapToReadShift + maxWidth * (laneToCheck) + b];
    const uint posToSample = mapToReadShift + kMaxMapWidthM * (laneToCheck + (((int) (b / kMaxMapWidthM)) * numLinesEdge)) + b % kMaxMapWidthM;
    laneChar = laneMap[posToSample];

    if (laneChar != 0xFF) {
      gap_a = ((float)b - initShift); //m
      v_a = laneChar; //laneChar is in 3*ms (to save space in array)
      found = true;
      break;
    }
  }

  if (found == false) {
    if (trafficLightState == 0x00) { //red
      //found=true;
      gap_a = gap_b = 1000.0f; //force to change to the line without vehicle
      v_a = v_b = 0xFF;
      return;
    }
  }

  if (found == false) {
    gap_a = 1000.0f;
  }

  // CHECK BACKWARDS
  found = false;

  //printf("2initShift %u numOfCells %u\n",initShift,numOfCells);
  for (int b = initShift + 1; (b >= 0) && (found == false); b--) {  // NOTE +1 to make sure there is none in at the same level
    //laneChar = laneMap[mapToReadShift + maxWidth * (laneToCheck) + b];
    const uint posToSample = mapToReadShift + kMaxMapWidthM * (laneToCheck + (((int) (b / kMaxMapWidthM)) * numLinesEdge)) + b % kMaxMapWidthM;
    laneChar = laneMap[posToSample];
    if (laneChar != 0xFF) {
      gap_b = ((float)initShift - b); //m
      v_b = laneChar; //laneChar is in 3*ms (to save space in array)
      found = true;
      break;
    }
  }

  //printf("3initShift %u numOfCells %u\n",initShift,numOfCells);
  if (found == false) {
    gap_b = 1000.0f;
  }

}//

void simulateOnePersonCPU(
  uint p,
  float deltaTime,
  float currentTime,
  uint mapToReadShift,
  uint mapToWriteShift,
  std::vector<LC::CUDATrafficPerson> &trafficPersonVec,
  std::vector<LC::B18EdgeData> &edgesData,
  std::vector<uchar> &laneMap,
  std::vector<B18IntersectionData> &intersections,
  std::vector<uchar> &trafficLights) {
  //if(DEBUG_TRAFFIC==1)printf("currentTime %f   0 Person: %d State %d Time Dep %f\n",currentTime,p,trafficPersonVec[p].active, trafficPersonVec[p].time_departure);
  ///////////////////////////////
  //2.0. check if finished
  if (trafficPersonVec[p].active == 2) {
    return;
  }

  ///////////////////////////////
  //2.1. check if person should still wait or should start
  if (trafficPersonVec[p].active == 0) {

    //printf("  1. Person: %d active==0\n",p);
    if (trafficPersonVec[p].time_departure > currentTime) { //wait
      //1.1 just continue waiting
      //printf("   1.1 Person: %d wait\n",p);
      return;
    } else { //start
      //1.2 find first edge
      trafficPersonVec[p].currPathEdge = 0;
      uint firstEdge = trafficPersonVec[p].personPath[0];

      if (firstEdge == -1) {
        trafficPersonVec[p].active = 2;
        //printf("0xFFFF\n");
        return;
      }

      if (DEBUG_TRAFFIC == 1) {
        printf("   1.2 Person: %d TRY put in first edge\n", p, firstEdge);
      }

      //1.3 update person edgeData
      //if(DEBUG_TRAFFIC==1)printf("   1.3 Person: %d put in first edge %u\n",p,firstEdge);
      //printf("edgesData %d\n",edgesData);

      // COPY DATA FROM EDGE TO PERSON
      trafficPersonVec[p].edgeNumLanes = edgesData[firstEdge].numLines;
      trafficPersonVec[p].edgeNextInters = edgesData[firstEdge].nextInters;

      trafficPersonVec[p].length = edgesData[firstEdge].length;
      trafficPersonVec[p].maxSpeedMperSec = edgesData[firstEdge].maxSpeedMperSec;
      //printf("edgesData %.10f\n",edgesData[firstEdge].maxSpeedCellsPerDeltaTime);
      //1.4 try to place it in middle of edge
      ushort numOfCells = ceil(trafficPersonVec[p].length);
      ushort initShift = (ushort)(0.5f * numOfCells); //number of cells it should be placed (half of road)

      uchar laneChar;
      bool placed = false;

      ushort numCellsEmptyToBePlaced = s_0;
      ushort countEmptyCells = 0;

      for (ushort b = initShift; (b < numOfCells) && (placed == false); b++) {
        ushort lN = trafficPersonVec[p].edgeNumLanes - 1; //just right LANE !!!!!!!
        laneChar = laneMap[mapToReadShift + kMaxMapWidthM * (firstEdge + lN) +
                                          b]; //get byte of edge (proper line)

        if (laneChar != 0xFF) {
          countEmptyCells = 0;
          continue;
        }

        countEmptyCells++;// ensure there is enough room to place the car

        if (countEmptyCells < numCellsEmptyToBePlaced) {
          continue;
        }

        trafficPersonVec[p].numOfLaneInEdge = lN;
        trafficPersonVec[p].posInLaneM = b; //m
        uchar vInMpS = (uchar)(trafficPersonVec[p].v *
                               3); //speed in m/s *3 (to keep more precision
        laneMap[mapToWriteShift + kMaxMapWidthM * (firstEdge + lN) + b] = vInMpS;
        placed = true;
        //printf("Placed\n");
        break;
        //}
      }

      if (placed == false) { //not posible to start now
        return;
      }

      trafficPersonVec[p].v = 0; 
      trafficPersonVec[p].LC_stateofLaneChanging = 0;

      //1.5 active car
      if (DEBUG_TRAFFIC == 1) {
        printf("   1.2 Person: %d PUT in first edge %u Pos %f of %f\n", p, firstEdge,
               trafficPersonVec[p].posInLaneM, trafficPersonVec[p].length);
      }

      trafficPersonVec[p].active = 1;
      trafficPersonVec[p].isInIntersection = 0;
      trafficPersonVec[p].num_steps = 1;
      trafficPersonVec[p].gas = 0;
      //trafficPersonVec[p].nextPathEdge++;//incremet so it continues in next edge
      // set up next edge info
      uint nextEdge = trafficPersonVec[p].personPath[1];

      //trafficPersonVec[p].nextEdge=nextEdge;
      if (nextEdge != -1) {
        trafficPersonVec[p].nextEdgemaxSpeedMperSec =
          edgesData[nextEdge].maxSpeedMperSec;
        trafficPersonVec[p].nextEdgeNumLanes = edgesData[nextEdge].numLines;
        trafficPersonVec[p].nextEdgeNextInters = edgesData[nextEdge].nextInters;
        trafficPersonVec[p].nextEdgeLength = edgesData[nextEdge].length;
        //trafficPersonVec[p].nextPathEdge++;
        trafficPersonVec[p].LC_initOKLanes = 0xFF;
        trafficPersonVec[p].LC_endOKLanes = 0xFF;
      }
      return;
    }
  }

  if (DEBUG_TRAFFIC == 1) {
    printf("    2. Person: %d moving\n", p);
  }

  ///////////////////////////////
  //2. it is moving
  trafficPersonVec[p].num_steps++;
  //2.1 try to move
  float numMToMove;
  bool getToNextEdge = false;
  bool nextVehicleIsATrafficLight = false;
  uint currentEdge = trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge];
  uint nextEdge = trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge + 1];
  if (DEBUG_TRAFFIC == 1) {
    printf("    2. Person: %d Try to move current Edge %u Next %u\n", p, currentEdge, nextEdge);
  }
  // www.vwi.tu-dresden.de/~treiber/MicroApplet/IDM.html
  // IDM
  float thirdTerm = 0;
  ///////////////////////////////////////////////////
  // 2.1.1 Find front car
  int numCellsCheck = std::max<float>(30.0f, trafficPersonVec[p].v * deltaTime * 2); //30 or double of the speed*time

  // a) SAME LINE (BEFORE SIGNALING)
  bool found = false;
  bool noFirstInLaneBeforeSign = false; //use for stop control (just let 1st to pass)
  bool noFirstInLaneAfterSign = false; //use for stop control (just let 1st to pass)
  float s;
  float delta_v;
  uchar laneChar;
  ushort byteInLine = (ushort)floor(trafficPersonVec[p].posInLaneM);
  ushort numOfCells = ceil((trafficPersonVec[p].length - intersectionClearance));

  for (ushort b = byteInLine + 2; (b < numOfCells) && (found == false) && (numCellsCheck > 0); b++, numCellsCheck--) {
    // laneChar = laneMap[mapToReadShift + maxWidth * (trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge] + trafficPersonVec[p].numOfLaneInEdge) + b];
    // ShiftRead + WIDTH * (width number * # edges + # laneInEdge) + b 
    const uint posToSample = mapToReadShift + kMaxMapWidthM * (trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge] + (((int) (byteInLine / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + b % kMaxMapWidthM;
    laneChar = laneMap[posToSample];

    if (laneChar != 0xFF) {
      s = ((float)(b - byteInLine)); //m
      delta_v = trafficPersonVec[p].v - (laneChar / 3.0f); //laneChar is in 3*ms (to save space in array)
      found = true;
      noFirstInLaneBeforeSign = true;
      if (DEBUG_TRAFFIC == 1) {
        printf("    2. Person: %d BEFORE Intersection: first obstacle -> car\n", p);
      }
      break;
    }
  }

  // b) TRAFFIC LIGHT
  if (byteInLine < numOfCells && found == false && numCellsCheck > 0) { //before traffic signaling (and not cell limited)
    if (trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] == 0x00) { //red
      s = ((float)(numOfCells - byteInLine)); //m
      delta_v = trafficPersonVec[p].v - 0; //it should be treated as an obstacle
      nextVehicleIsATrafficLight = true;
      if (DEBUG_TRAFFIC == 1) {
        printf("    2. Person: %d first obstacle -> traffic light\n", p);
      }
      //printf("\nFOUND TL\n",s,delta_v);
      found = true;
    }
  }

  // c) SAME LINE (AFTER SIGNALING)
  for (ushort b = byteInLine + 2; (b < numOfCells) && (found == false) && (numCellsCheck > 0); b++, numCellsCheck--) {
    // laneChar = laneMap[mapToReadShift + maxWidth * (trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge] + trafficPersonVec[p].numOfLaneInEdge) + b];
    const uint posToSample = mapToReadShift + kMaxMapWidthM * (trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge] + (((int) (byteInLine / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + b % kMaxMapWidthM;
    laneChar = laneMap[posToSample];

    if (laneChar != 0xFF) {
      s = ((float)(b - byteInLine)); //m
      delta_v = trafficPersonVec[p].v - (laneChar /
                                         3.0f); //laneChar is in 3*ms (to save space in array)
      found = true;
      noFirstInLaneAfterSign = true;
      if (DEBUG_TRAFFIC == 1) {
        printf("    2. Person: %d AFTER Intersection: first obstacle -> car\n", p);
      }
      break;
    }
  }

  if (trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] == 0x0F && numCellsCheck > 0) { //stop 
    //check
    if (noFirstInLaneBeforeSign == false && byteInLine < numOfCells && //first before traffic
        trafficPersonVec[p].v == 0 && //stopped
        noFirstInLaneAfterSign == false) { // noone after the traffic light (otherwise wait before stop) !! TODO also check the beginning of next edge

      trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] = 0x00; //reset stop
      trafficPersonVec[p].posInLaneM = ceil(numOfCells) + 1; //move magicly after stop

      if (DEBUG_TRAFFIC == 1) {
        printf("    2. Person: %d move after stop\n", p);
      }

    } else { //stop before STOP
      if (noFirstInLaneBeforeSign == false) { //just update this if it was the first one before sign
        s = ((float)(numOfCells - byteInLine)); //m
        delta_v = trafficPersonVec[p].v - 0; //it should be treated as an obstacle
        nextVehicleIsATrafficLight = true;
        found = true;
        if (DEBUG_TRAFFIC == 1) {
          printf("    2. Person: %d just update this if it was the first one before sign\n", p);
        }
      }
    }
  }

  // NEXT LINE
  if (found == false && numCellsCheck > 0) { //check if in next line
    if ((nextEdge != -1) &&
        (trafficPersonVec[p].edgeNextInters != trafficPersonVec[p].end_intersection)) { // we haven't arrived to destination (check next line)
      if (DEBUG_TRAFFIC == 1) {
        printf("    2. Person: NEXT LINE\n", p);
      }
      ushort nextEdgeLaneToBe = trafficPersonVec[p].numOfLaneInEdge; //same lane

      //printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
      if (nextEdgeLaneToBe >= trafficPersonVec[p].nextEdgeNumLanes) {
        nextEdgeLaneToBe = trafficPersonVec[p].nextEdgeNumLanes -
                           1; //change line if there are less roads
      }

      //printf("2trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
      ushort numOfCells = ceil(trafficPersonVec[p].nextEdgeLength);

      for (ushort b = 0; (b < numOfCells) && (found == false) && (numCellsCheck > 0); b++, numCellsCheck--) {
        //laneChar = laneMap[mapToReadShift + maxWidth * (nextEdge + nextEdgeLaneToBe) + b];
        const uint posToSample = mapToReadShift + kMaxMapWidthM * (nextEdge + nextEdgeLaneToBe) + b; // b18 not changed since we check first width
        laneChar = laneMap[posToSample];

        if (laneChar != 0xFF) {
          s = ((float)(b)); //m
          delta_v = trafficPersonVec[p].v - (laneChar / 3.0f);  // laneChar is in 3*ms (to save space in array)
          found = true;
          break;
        }
      }
    }
  }


  float s_star;

  if (found == true) { //car in front and slower than us
    // 2.1.2 calculate dv_dt
    s_star = s_0 + std::max(0.0f,
                            (trafficPersonVec[p].v * trafficPersonVec[p].T + (trafficPersonVec[p].v *
                                delta_v) / (2 * std::sqrt(trafficPersonVec[p].a * trafficPersonVec[p].b))));
    thirdTerm = std::pow(((s_star) / (s)), 2);
    //printf(">FOUND s_star %f thirdTerm %f!!!!\n",s_star,thirdTerm);
  }

  float dv_dt = trafficPersonVec[p].a * (1.0f - std::pow((
      trafficPersonVec[p].v / trafficPersonVec[p].maxSpeedMperSec), 4) - thirdTerm);

  // 2.1.3 update values
  numMToMove = std::max(0.0f,
                        trafficPersonVec[p].v * deltaTime + 0.5f * (dv_dt) * deltaTime * deltaTime);

  //if(p==4524){
  //	printf("v %f v0 %f a %f dv_dt %f s %f s_star %f thirdTerm %f MOVE %f\n",trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec,trafficPersonVec[p].a,dv_dt,s,s_star,thirdTerm,numMToMove);
  //}
  //if(numMToMove==0){
  //printf("v %f v0 %f a %f dv_dt %f s %f s_star %f thirdTerm %f MOVE %f\n",trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec,trafficPersonVec[p].a,dv_dt,s,s_star,thirdTerm,numMToMove);
  //}
  if (DEBUG_TRAFFIC == 1) {
    printf("v %f v0 %f a %f dv_dt %f MOVE %f\n", trafficPersonVec[p].v,
           trafficPersonVec[p].maxSpeedMperSec, trafficPersonVec[p].a, dv_dt, numMToMove);
  }

  //printf("v %.10f v d %.10f\n",trafficPersonVec[p].v,trafficPersonVec[p].v+((dv_dt/(deltaTime)/deltaTime)));
  trafficPersonVec[p].v += dv_dt * deltaTime;

  if (trafficPersonVec[p].v < 0) {
    //printf("p %d v %f v0 %f a %f dv_dt %f s %f s_star %f MOVE %f\n",p,trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec,trafficPersonVec[p].a,dv_dt,s,s_star,numMToMove);
    trafficPersonVec[p].v = 0;
  }

  /////
  //CO2
  //if(trafficPersonVec[p].v>0)
  if (calculatePollution) {
    float speedMph = trafficPersonVec[p].v * 2.2369362920544; //mps to mph
    float gasStep = -0.064 + 0.0056 * speedMph + 0.00026 * (speedMph - 50.0f) *
                    (speedMph - 50.0f);

    if (gasStep > 0) {
      gasStep *= deltaTime;
      trafficPersonVec[p].gas += gasStep;
    }
  }
  //trafficPersonVec[p].gas+=numMToMove/1000.0f;
  //////////////////////////////////////////////

  if (trafficPersonVec[p].v == 0) { //if not moving not do anything else
    ushort posInLineCells = (ushort)(trafficPersonVec[p].posInLaneM);
    //laneMap[mapToWriteShift + maxWidth * (currentEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = 0;
    const uint posToSample = mapToWriteShift + kMaxMapWidthM * (currentEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;
    laneMap[posToSample] = 0;
    
    return;
  }

  //////////

  ///////////////////////////////
  // COLOR
  trafficPersonVec[p].color = p << 8;
  //if (clientMain->ui.b18RenderSimulationCheckBox->isChecked()) {
    //if(G::global().getInt("cuda_carInfoRendering_type")==0){
    //qsrand(p);
    
    /*}
    if(G::global().getInt("cuda_carInfoRendering_type")==1){
        uchar c=(uchar)(255*trafficPersonVec[p].v/15.0f);//84m/s is more than 300km/h
        trafficPersonVec[p].color=(c<<24)|(c<<16)|(c<<8);
    }
    if(G::global().getInt("cuda_carInfoRendering_type")==2){
        uchar c=255*trafficPersonVec[p].LC_stateofLaneChanging;
        trafficPersonVec[p].color=(c<<24)|(c<<16)|(c<<8);

    }*/
  //}

  ////////////////////////////////

  if (DEBUG_TRAFFIC == 1) {
    printf("2 v(t+de) %f\n\n", trafficPersonVec[p].v);
  }

  // STOP (check if it is a stop if it can go through)

  trafficPersonVec[p].posInLaneM = trafficPersonVec[p].posInLaneM + numMToMove;

  if (trafficPersonVec[p].posInLaneM >
      trafficPersonVec[p].length) { //reach intersection
    numMToMove = trafficPersonVec[p].posInLaneM - trafficPersonVec[p].length;
    getToNextEdge = true;
  } else { //does not research next intersection
    ////////////////////////////////////////////////////////
    // LANE CHANGING (happens when we are not reached the intersection)
    if (trafficPersonVec[p].v > 3.0f && //at least 10km/h to try to change lane
        trafficPersonVec[p].num_steps % 10 == 0 //just check every (10 steps) 5 seconds
       ) {
      //next thing is not a traffic light
      // skip if there is one lane (avoid to do this)
      // skip if it is the last edge
      if (nextVehicleIsATrafficLight == false &&
          trafficPersonVec[p].edgeNumLanes > 1 && nextEdge != -1) {

        ////////////////////////////////////////////////////
        // LC 1 update lane changing status
        if (trafficPersonVec[p].LC_stateofLaneChanging == 0) {
          // 2.2-exp((x-1)^2)
          float x = trafficPersonVec[p].posInLaneM / trafficPersonVec[p].length;

          if (x > 0.4f) { //just after 40% of the road
            float probabiltyMandatoryState = 2.2 - exp((x - 1) * (x - 1));

            if (((float)qrand() / RAND_MAX) < probabiltyMandatoryState) {
              trafficPersonVec[p].LC_stateofLaneChanging = 1;
            }
          }

        }

        ////////////////////////////////////////////////////
        // LC 2 NOT MANDATORY STATE
        if (trafficPersonVec[p].LC_stateofLaneChanging == 0) {
          //if(p==40)printf("LC v %f v0 %f a %f\n",trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec*0.5f,dv_dt);
          // discretionary change: v slower than the current road limit and deccelerating and moving
          if ((trafficPersonVec[p].v < (trafficPersonVec[p].maxSpeedMperSec * 0.7f)) &&
              (dv_dt < 0) && trafficPersonVec[p].v > 3.0f) {
            //printf(">>LANE CHANGE\n");
            
            //printf("LC 0 %u\n",trafficPersonVec[p].numOfLaneInEdge);
            bool leftLane = trafficPersonVec[p].numOfLaneInEdge >
                            0; //at least one lane on the left
            bool rightLane = trafficPersonVec[p].numOfLaneInEdge <
                             trafficPersonVec[p].edgeNumLanes - 1; //at least one lane

            if (leftLane == true && rightLane == true) {
              if (qrand() % 2 == 0) {
                leftLane = false;
              } else {
                rightLane = false;
              }
            }
            ushort laneToCheck;
            if (leftLane == true) {
              laneToCheck = trafficPersonVec[p].numOfLaneInEdge - 1;
            } else {
              laneToCheck = trafficPersonVec[p].numOfLaneInEdge + 1;
            }

            uchar v_a, v_b;
            float gap_a, gap_b;
            //printf("p %u LC 1 %u\n",p,laneToCheck);
            uchar trafficLightState = trafficLights[currentEdge +
                                                                trafficPersonVec[p].numOfLaneInEdge];
            calculateGapsLC(mapToReadShift, laneMap, trafficLightState,
                            currentEdge + laneToCheck, trafficPersonVec[p].edgeNumLanes, trafficPersonVec[p].posInLaneM,
                            trafficPersonVec[p].length, v_a, v_b, gap_a, gap_b);

            //printf("LC 2 %u %u %f %f\n",v_a,v_b,gap_a,gap_b);
            if (gap_a == 1000.0f && gap_b == 1000.0f) { //lag and lead car very far
              trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE

            } else { // NOT ALONE
              float b1A = 0.05f, b2A = 0.15f;
              float b1B = 0.15f, b2B = 0.40f;
              // s_0-> critical lead gap
              float g_na_D, g_bn_D;
              bool acceptLC = true;

              if (gap_a != 1000.0f) {
                g_na_D = std::max(s_0, s_0 + b1A * trafficPersonVec[p].v + b2A *
                                  (trafficPersonVec[p].v - v_a * 3.0f));

                if (gap_a < g_na_D) { //gap smaller than critical gap
                  acceptLC = false;
                }
              }

              if (acceptLC == true && gap_b != 1000.0f) {
                g_bn_D = std::max(s_0, s_0 + b1B * v_b * 3.0f + b2B * (v_b * 3.0f -
                                  trafficPersonVec[p].v));

                if (gap_b < g_bn_D) { //gap smaller than critical gap
                  acceptLC = false;
                }
              }

              if (acceptLC == true) {
                trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE
              }
            }

            //printf("<<LANE CHANGE\n");
          }


        }// Discretionary

        ////////////////////////////////////////////////////
        // LC 3 *MANDATORY* STATE
        if (trafficPersonVec[p].LC_stateofLaneChanging == 1) {
          // LC 3.1 Calculate the correct lanes
          if (trafficPersonVec[p].LC_endOKLanes == 0xFF) {
            calculateLaneCarShouldBe(currentEdge, nextEdge, intersections,
                                     trafficPersonVec[p].edgeNextInters, trafficPersonVec[p].edgeNumLanes,
                                     trafficPersonVec[p].LC_initOKLanes, trafficPersonVec[p].LC_endOKLanes);

            //printf("p%u num lanes %u min %u max %u\n",p,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
            if (trafficPersonVec[p].LC_initOKLanes == 0 &&
                trafficPersonVec[p].LC_endOKLanes == 0) {
              exit(0);
            }
          }


          //printf(">>LANE CHANGE\n");
          //printf("LC 0 %u\n",trafficPersonVec[p].numOfLaneInEdge);
          bool leftLane = false, rightLane = false;

          // LC 3.2 CORRECT LANES--> DICRETIONARY LC WITHIN
          if (trafficPersonVec[p].numOfLaneInEdge >= trafficPersonVec[p].LC_initOKLanes &&
              trafficPersonVec[p].numOfLaneInEdge < trafficPersonVec[p].LC_endOKLanes) {
            // for discretionary it should be under some circustances
            if ((trafficPersonVec[p].v < (trafficPersonVec[p].maxSpeedMperSec * 0.7f)) &&
                (dv_dt < 0) && trafficPersonVec[p].v > 3.0f) {
              leftLane =
                (trafficPersonVec[p].numOfLaneInEdge > 0) && //at least one lane on the left
                (trafficPersonVec[p].numOfLaneInEdge - 1 >= trafficPersonVec[p].LC_initOKLanes)
                &&
                (trafficPersonVec[p].numOfLaneInEdge - 1 < trafficPersonVec[p].LC_endOKLanes);
              rightLane =
                (trafficPersonVec[p].numOfLaneInEdge < trafficPersonVec[p].edgeNumLanes - 1) &&
                //at least one lane
                (trafficPersonVec[p].numOfLaneInEdge + 1 >= trafficPersonVec[p].LC_initOKLanes)
                &&
                (trafficPersonVec[p].numOfLaneInEdge + 1 < trafficPersonVec[p].LC_endOKLanes);
              //printf("D\n");
            }
          }
          // LC 3.3 INCORRECT LANES--> MANDATORY LC
          else {
            //printf("num lanes %u min %u max %u\n",trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
            //printf("p%u num lanes %u min %u max %u\n",p,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);

            if (trafficPersonVec[p].numOfLaneInEdge < trafficPersonVec[p].LC_initOKLanes) {
              rightLane = true;
            } else {
              leftLane = true;
            }

            if (rightLane == true &&
                trafficPersonVec[p].numOfLaneInEdge + 1 >= trafficPersonVec[p].edgeNumLanes) {
              printf("ERROR: RT laneToCheck>=trafficPersonVec[p].edgeNumLanes\n");
            }

            if (leftLane == true && trafficPersonVec[p].numOfLaneInEdge == 0) {
              printf("ERROR %u: LT laneToCheck>=trafficPersonVec[p].edgeNumLanes OK %u-%u NE %u\n",
                     p, trafficPersonVec[p].LC_initOKLanes, trafficPersonVec[p].LC_endOKLanes,
                     nextEdge);
              exit(0);
            }

            //printf("M L %d R %d nL %u\n",leftLane,rightLane,trafficPersonVec[p].numOfLaneInEdge);
          }

          if (leftLane == true || rightLane == true) {

            // choose lane (if necessary)
            if (leftLane == true && rightLane == true) {
              if (qrand() % 2 == 0) {
                leftLane = false;
              } else {
                rightLane = false;
              }
            }
            ushort laneToCheck;
            if (leftLane == true) {
              laneToCheck = trafficPersonVec[p].numOfLaneInEdge - 1;
            } else {
              laneToCheck = trafficPersonVec[p].numOfLaneInEdge + 1;
            }

            if (laneToCheck >= trafficPersonVec[p].edgeNumLanes) {
              printf("ERROR: laneToCheck>=trafficPersonVec[p].edgeNumLanes %u %u\n",
                     laneToCheck, trafficPersonVec[p].edgeNumLanes);
            }

            uchar v_a, v_b;
            float gap_a, gap_b;
            //printf("p %u LC 1 %u\n",p,laneToCheck);
            uchar trafficLightState = trafficLights[currentEdge +
                                                                trafficPersonVec[p].numOfLaneInEdge];
            calculateGapsLC(mapToReadShift, laneMap, trafficLightState,
                            currentEdge + laneToCheck, trafficPersonVec[p].edgeNumLanes, trafficPersonVec[p].posInLaneM,
                            trafficPersonVec[p].length, v_a, v_b, gap_a, gap_b);

            //printf("LC 2 %u %u %f %f\n",v_a,v_b,gap_a,gap_b);
            if (gap_a == 1000.0f && gap_b == 1000.0f) { //lag and lead car very far
              trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE

            } else { // NOT ALONE
              float b1A = 0.05f, b2A = 0.15f;
              float b1B = 0.15f, b2B = 0.40f;
              float gamma = 0.000025;
              // s_0-> critical lead gap
              float distEnd = trafficPersonVec[p].length - trafficPersonVec[p].posInLaneM;
              float expTerm = (1 - exp(-gamma * distEnd * distEnd));

              float g_na_M, g_bn_M;
              bool acceptLC = true;

              if (gap_a != 1000.0f) {
                g_na_M = std::max(s_0, s_0 + (b1A * trafficPersonVec[p].v + b2A *
                                              (trafficPersonVec[p].v - v_a * 3.0f)));

                if (gap_a < g_na_M) { //gap smaller than critical gap
                  acceptLC = false;
                }
              }

              if (acceptLC == true && gap_b != 1000.0f) {
                g_bn_M = std::max(s_0, s_0 + (b1B * v_b * 3.0f + b2B * (v_b * 3.0f -
                                              trafficPersonVec[p].v)));

                if (gap_b < g_bn_M) { //gap smaller than critical gap
                  acceptLC = false;
                }
              }

              if (acceptLC == true) {
                trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE
              }
            }


          }

        }// Mandatory

      }//at least two lanes and not stopped by traffic light

    }

    ///////////////////////////////////////////////////////
    if (DEBUG_TRAFFIC == 1) {
      printf("    2. Person: %d moving in edge %u pos %f of %f\n", p, currentEdge,
             trafficPersonVec[p].posInLaneM, trafficPersonVec[p].length);
    }

    uchar vInMpS = (uchar)(trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
    ushort posInLineCells = (ushort)(trafficPersonVec[p].posInLaneM);
    //laneMap[mapToWriteShift + maxWidth * (currentEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = vInMpS;
    const uint posToSample = mapToWriteShift + kMaxMapWidthM * (currentEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;
    laneMap[posToSample] = vInMpS;
    //printf("2<<LANE CHANGE\n");
    return;
  }

  //}
  //2.2 close to intersection

  //2.2 check if change intersection
  //!!!ALWAYS CHANGE
  //2.2.1 find next edge
  /*ushort curr_intersection=trafficPersonVec[p].edgeNextInters;
  ushort end_intersection=trafficPersonVec[p].end_intersection;
  //2.1 check if end*/
  if (nextEdge == -1) { //if(curr_intersection==end_intersection){
    if (DEBUG_TRAFFIC == 1) {
      printf("    2.1 Person: %d FINISHED\n", p);
    }

    trafficPersonVec[p].active = 2; //finished
    return;
  }

  //if(trafficPersonVec[p].nextPathEdge>=nextEdgeM.size())printf("AAAAAAAAAAAAAAAAA\n");
  /////////////
  // update edge
  /*// stop
  if(noFirstInLane==false&&trafficLights[currentEdge+trafficPersonVec[p].numOfLaneInEdge]==0x0F){
        // first in lane and stop--> update to avoid to pass another car
        trafficLights[currentEdge+trafficPersonVec[p].numOfLaneInEdge]=0x00;
  }*/
  //trafficPersonVec[p].curEdgeLane=trafficPersonVec[p].nextEdge;
  trafficPersonVec[p].currPathEdge++;
  trafficPersonVec[p].maxSpeedMperSec = trafficPersonVec[p].nextEdgemaxSpeedMperSec;
  trafficPersonVec[p].edgeNumLanes = trafficPersonVec[p].nextEdgeNumLanes;
  trafficPersonVec[p].edgeNextInters = trafficPersonVec[p].nextEdgeNextInters;
  trafficPersonVec[p].length = trafficPersonVec[p].nextEdgeLength;
  trafficPersonVec[p].posInLaneM = numMToMove;

  if (trafficPersonVec[p].numOfLaneInEdge >= trafficPersonVec[p].edgeNumLanes) {
    trafficPersonVec[p].numOfLaneInEdge = trafficPersonVec[p].edgeNumLanes - 1; //change line if there are less roads
  }

  ////////////
  // update next edge
  uint nextNEdge =
    trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge + 1];

  //trafficPersonVec[p].nextEdge=nextEdge;
  if (nextNEdge != -1) {
    //trafficPersonVec[p].nextPathEdge++;
    trafficPersonVec[p].LC_initOKLanes = 0xFF;
    trafficPersonVec[p].LC_endOKLanes = 0xFF;

    if (DEBUG_TRAFFIC == 1) {
      printf("    2.2.1 Person: %d curr %u end %u nextEdge %u\n", p,
             trafficPersonVec[p].edgeNextInters, trafficPersonVec[p].end_intersection,
             nextNEdge);
    }

    //2.2.3 update person edgeData
    //trafficPersonVec[p].nextEdge=nextEdge;
    trafficPersonVec[p].nextEdgemaxSpeedMperSec =
      edgesData[nextNEdge].maxSpeedMperSec;
    trafficPersonVec[p].nextEdgeNumLanes = edgesData[nextNEdge].numLines;
    trafficPersonVec[p].nextEdgeNextInters = edgesData[nextNEdge].nextInters;
    trafficPersonVec[p].nextEdgeLength = edgesData[nextNEdge].length;
  }

  trafficPersonVec[p].LC_stateofLaneChanging = 0;
  uchar vInMpS = (uchar)(trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
  ushort posInLineCells = (ushort)(trafficPersonVec[p].posInLaneM);

  // laneMap[mapToWriteShift + maxWidth * (nextEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = vInMpS;
  const uint posToSample = mapToWriteShift + kMaxMapWidthM * (nextEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;  // note the last % should not happen
  laneMap[posToSample] = vInMpS;


  if (DEBUG_TRAFFIC == 1) {
    printf("    2.2.4 Person: %d New Lane %u Pos %f\n", p, nextEdge, numMToMove);
  }

}//

void simulateOneSTOPIntersectionCPU(
  uint i,
  float deltaTime,
  float currentTime, std::vector<B18IntersectionData> &intersections,
  std::vector<uchar> &trafficLights,
  std::vector<LC::B18EdgeData> &edgesData,//for the length
  std::vector<uchar> &laneMap,//to check if there are cars
  uint mapToReadShift) {
  const float deltaEvent = 0; //10.0f;!!!! CHANGE

  //if(i==0)printf("i %d\n",i);
  if (currentTime > intersections[i].nextEvent && intersections[i].totalInOutEdges > 0) {
    uint edgeOT = intersections[i].edge[intersections[i].state];
    uchar numLinesO = edgeOT >> 24;
    uint edgeONum = edgeOT & kMaskLaneMap; // 0xFFFFF

    // red old traffic lights
    for (int nL = 0; nL < numLinesO; nL++) {
      trafficLights[edgeONum + nL] = 0x00; //red old traffic light
    }

    for (int iN = 0; iN <= intersections[i].totalInOutEdges + 1;
         iN++) { //to give a round
      intersections[i].state = (intersections[i].state + 1) %
                               intersections[i].totalInOutEdges;//next light

      if ((intersections[i].edge[intersections[i].state] & kMaskInEdge) == kMaskInEdge) {  // 0x800000
        uint edgeIT = intersections[i].edge[intersections[i].state];
        uint edgeINum = edgeIT & kMaskLaneMap; //get edgeI 0xFFFFF
        uchar numLinesI = edgeIT >> 24;
        /// check if someone in this edge
        int rangeToCheck = 5.0f; //5m
        ushort firstPosToCheck = edgesData[edgeINum].length - intersectionClearance; //last po
        bool atLeastOneStopped = false;

        for (int posCheck = firstPosToCheck; rangeToCheck >= 0 && posCheck >= 0; posCheck--, rangeToCheck--) { //as many cells as the rangeToCheck says
          for (int nL = 0; nL < numLinesI; nL++) {
            //int cellNum = mapToReadShift + maxWidth * (edgeINum + nL) + posCheck;
            const uint posToSample = mapToReadShift + kMaxMapWidthM * (edgeINum + (((int) (posCheck / kMaxMapWidthM)) * numLinesI) + nL) + posCheck % kMaxMapWidthM;


            if (laneMap[posToSample] == 0) { //car stopped
              trafficLights[edgeINum + nL] = 0x0F; // STOP SIGN 0x0F--> Let pass
              atLeastOneStopped = true;
            }
          }
        }

        if (atLeastOneStopped == true) {
          intersections[i].nextEvent = currentTime + deltaEvent; //just move forward time if changed (otherwise check in next iteration)
          break;
        }
      }
    }
  }
}//

void simulateOneIntersectionCPU(uint i, float currentTime,
                                std::vector<B18IntersectionData> &intersections,
                                std::vector<uchar> &trafficLights) {
  const float deltaEvent = 20.0f; /// !!!!

  //if (DEBUG_T_LIGHT == 1) printf("Inter %d CurrTime %.1f Next Event %.1f InOut %d\n", i, currentTime, intersections[i].nextEvent, intersections[i].totalInOutEdges);
  if (currentTime > intersections[i].nextEvent && intersections[i].totalInOutEdges > 0) {
    
    uint edgeOT = intersections[i].edge[intersections[i].state];
    uchar numLinesO = edgeOT >> 24;
    uint edgeONum = edgeOT & kMaskLaneMap; // 0xFFFFF;

    if (DEBUG_T_LIGHT == 1) {
      printf("Inter %d: Event happen: numLines %u laneMap %u\n", i, numLinesO, edgeONum);
    }

    // red old traffic lights
    if ((edgeOT&kMaskInEdge) == kMaskInEdge) {  // Just do it if we were in in
      for (int nL = 0; nL < numLinesO; nL++) {
        trafficLights[edgeONum + nL] = 0x00; //red old traffic light
        if (DEBUG_T_LIGHT == 1) {
          printf("Inter %d: Event happen: numLines %u laneMap %u --> Set red traffic light %u\n", i, numLinesO, edgeONum, edgeONum + nL);
        }
      }
    }

    for (int iN = 0; iN <= intersections[i].totalInOutEdges + 1; iN++) { //to give a round
      intersections[i].state = (intersections[i].state + 1) % intersections[i].totalInOutEdges;//next light

      if ((intersections[i].edge[intersections[i].state] & kMaskInEdge) == kMaskInEdge) {  // 0x800000
        // green new traffic lights
        uint edgeIT = intersections[i].edge[intersections[i].state];
        uint edgeINum = edgeIT & kMaskLaneMap; //  0xFFFFF; //get edgeI
        uchar numLinesI = edgeIT >> 24;

        for (int nL = 0; nL < numLinesI; nL++) {
          trafficLights[edgeINum + nL] = 0xFF;
          if (DEBUG_T_LIGHT == 1) {
            printf("Inter %d: Event happen: numLines %u laneMap %u --> Set green traffic light %u\n", i, numLinesO, edgeONum, edgeINum + nL);
          }
        }

        //trafficLights[edgeINum]=0xFF;
        break;
      }
    }//green new traffic light

    //printf("i %d CHANGE state %u of %d (Old edge %u New Edge %u)\n",i,intersections[i].state,intersections[i].totalInOutEdges,edgeO,edgeI);
    ////
    intersections[i].nextEvent = currentTime + deltaEvent;
  }
}//

void sampleTraffic(std::vector<CUDATrafficPerson> &trafficPersonVec,
                   std::vector<float> &accSpeedPerLinePerTimeInterval,
                   std::vector<float> &numVehPerLinePerTimeInterval, uint offset) {
  int numPeople = trafficPersonVec.size();

  //printf("offset %d\n",offset);
  for (int p = 0; p < numPeople; p++) {
    if (trafficPersonVec[p].active == 1) {
      int lineNum = trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge];
      accSpeedPerLinePerTimeInterval[lineNum + offset] += trafficPersonVec[p].v / 3.0f;
      numVehPerLinePerTimeInterval[lineNum + offset]++;
    }
  }//for people
}//

// numOfPasses-> define if just disjktra or iterative
// reGeneratePeopleLanes-> recompute lanes (it is used in MCMC that calls those func before)
void B18TrafficSimulator::simulateInCPU_MultiPass(int numOfPasses, float startTimeH, float endTimeH, bool useJohnsonRouting) {

  // create lane map
  if (DEBUG_SIMULATOR) {
    printf(">>simulateInCPU_MultiPass\n");
  }
  createLaneMap();

  if (DEBUG_SIMULATOR) {
    printf("MP Start multi pass simulation\n");
  }
  int weigthMode;
  float peoplePathSampling[] = {1.0f, 1.0f, 0.5f, 0.25f, 0.12f, 0.67f};

  for (int nP = 0; nP < numOfPasses; nP++) {
    resetPeopleJobANDintersections();//reset people
    weigthMode = 1; //first time run normal weights

    if (nP == 0) {
      weigthMode = 0;
    }

    if (DEBUG_SIMULATOR) {
      printf("Num edges %d\n", boost::num_edges(simRoadGraph->myRoadGraph_BI));
    }

    QTime pathTimer;
    pathTimer.start();
    
    if (useJohnsonRouting) {
      printf("***Start generateRoute Johnson\n");
      B18TrafficJohnson::generateRoutes(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, weigthMode, peoplePathSampling[nP]);
    } else {
      printf("***Start generateRoutesMulti Disktra\n");
      cudaTrafficPersonShortestPath.generateRoutesMulti(simRoadGraph->myRoadGraph_BI, trafficPersonVec, edgeDescToLaneMapNum, weigthMode, peoplePathSampling[nP]);
    }
    printf("***Path time %f\n",pathTimer.elapsed());
    // run simulation
    printf("***Start simulateInCPU \n");
    simulateInCPU(startTimeH, endTimeH);
    printf("***End simulateInCPU \n");
    //estimate traffic density
    calculateAndDisplayTrafficDensity();
  }
  if (DEBUG_SIMULATOR) {
    printf("<<simulateInCPU_MultiPass\n");
  }
}//

void B18TrafficSimulator::simulateInCPU_Onepass(float startTimeH,
  float endTimeH, bool useJohnsonRouting) {
  resetPeopleJobANDintersections();
  // create lane map
  createLaneMap();
  // find path for each vehicle
  generateCarPaths(useJohnsonRouting);
  // run simulation
  simulateInCPU(startTimeH, endTimeH);
  //estimate traffic density
  calculateAndDisplayTrafficDensity();
}//


void B18TrafficSimulator::simulateInCPU(float startTimeH, float endTimeH) {
  float startTime = startTimeH * 3600.0f; //7.0f
  float endTime = endTimeH * 3600.0f; //8.0f//10.0f

  if (DEBUG_SIMULATOR) {
    printf(">>simulateInCPU\n");
  }

  float currentTime = 23.99f * 3600.0f;

  for (int p = 0; p < trafficPersonVec.size() ; p++) {
    //for(int p=40;p<41;p++){
    if (currentTime > trafficPersonVec[p].time_departure) {
      currentTime = trafficPersonVec[p].time_departure;
    }
  }

  //round time to be delta divisible
  int numInt = currentTime / deltaTime; //floor

  if (DEBUG_SIMULATOR) {
    printf("currentTime %.2fs (%.2fh)--> %f\n", currentTime, currentTime/3600.0f, numInt * deltaTime);
  }

  currentTime = numInt * deltaTime;
  uint steps = 0;
  steps = (currentTime - startTime) / deltaTime;

  // start as early as starting time
  if (currentTime < startTime) { //as early as the starting time
    currentTime = startTime;
  }

  // 0. put data as references to simplified code

  int numPeople = trafficPersonVec.size();
  int numIntersec = intersections.size();

  if (DEBUG_SIMULATOR) {
    printf(">>Start Simulation %.2fs (%.2fh)\n", currentTime, currentTime/3600.0f);
  }

  QTime timer;
  timer.start();
  bool readFirstMap = true;
  uint mapToReadShift;
  uint mapToWriteShift;
  uint halfLaneMap = laneMap.size() / 2; //it is multiple of 2

  // sampling
  uint numStepsPerSample = 30.0f / deltaTime; //each min
  const uint numStepsTogether = 12; //change also in density (10 per hour)
  uint numSamples = ceil(((endTime - startTime) / (deltaTime * numStepsPerSample * numStepsTogether))) + 1; //!!!

  if (DEBUG_SIMULATOR) {
    printf("numSamples %d\n", numSamples);
  }

  uint numLines = trafficLights.size();
  accSpeedPerLinePerTimeInterval.resize(numSamples * numLines); //!!
  numVehPerLinePerTimeInterval.resize(numSamples * numLines);
  memset(accSpeedPerLinePerTimeInterval.data(), 0, accSpeedPerLinePerTimeInterval.size()*sizeof(float));
  memset(numVehPerLinePerTimeInterval.data(), 0, numVehPerLinePerTimeInterval.size()*sizeof(float));
  int samplingNumber = 0;

  if (DEBUG_SIMULATOR) {
    printf("Map -->hlafMap /2 %u delta %f\n", halfLaneMap, deltaTime);
    printf("\nStart Simulation\n");
  }

  int count = 0;

  while (currentTime < endTime) {
    count++;
    if (count % 180 == 0) {
      printf("Time %.2fh (%.2f --> %.2f): %.0f%%\n", (currentTime / 3600.0f), (startTime / 3600.0f), (endTime / 3600.0f), 100.0f-(100.0f * (endTime - currentTime) / (endTime-startTime)));
    }
    //bool anyActive=false;
    ////////////////////////////////////////////////////////////
    // 1. CHANGE MAP: set map to use and clean the other
    if (DEBUG_SIMULATOR) {
      printf("Clean Map\n");
    }
    if (readFirstMap == true) {
      mapToReadShift = 0;
      mapToWriteShift = halfLaneMap;
      memset(&laneMap.data()[halfLaneMap], -1,
             halfLaneMap * sizeof(unsigned char)); //clean second half
    } else {
      mapToReadShift = halfLaneMap;
      mapToWriteShift = 0;
      memset(&laneMap.data()[0], -1,
             halfLaneMap * sizeof(unsigned char)); //clean first half
    }

    readFirstMap = !readFirstMap; //next iteration invert use


    ////////////////////////////////////////////////////////////
    // 2. UPDATE INTERSECTIONS
    printf("Update Intersections\n");
    for (int i = 0; i < intersections.size(); i++) {
      simulateOneIntersectionCPU(i, currentTime, intersections, trafficLights);
    }

    //printf("Sim people\n");
    ////////////////////////////////////////////////////////////
    // 3. UPDATE PEOPLE
    printf("Update People\n");
    for (int p = 0; p < numPeople; p++) {
      simulateOnePersonCPU(p, deltaTime, currentTime, mapToReadShift, mapToWriteShift, trafficPersonVec, edgesData, laneMap, intersections, trafficLights);
    }//for people

    ////////////////////////////////////////////////////////////
    // 4. UPDATE SAMPLIG
    //if(steps%numStepsPerSample==0){
    printf("Update Sampling\n");
    if ((((float)((int)currentTime)) == (currentTime)) &&
        ((int)currentTime % ((int)30)) == 0) { //3min //(sample double each 3min)
      samplingNumber = (currentTime - startTime) / (30 * numStepsTogether);
      uint offset = numLines * samplingNumber;
      sampleTraffic(trafficPersonVec, accSpeedPerLinePerTimeInterval,
                    numVehPerLinePerTimeInterval, offset);
      //if((((float)((int)currentTime))==(currentTime))&&((int)currentTime%((int)30*60))==0)//each hour
      int cTH = currentTime / 3600;
      int cTM = (currentTime - cTH * 3600) / 60;

      //int cTS=(currentTime-cTH*3600-cTM*60);
      if (cTM % 20 == 0) {
        printf("T[%d] Time %f (%d:%02d) samplingNumber %d\n", threadNumber, currentTime,
               cTH, cTM, samplingNumber);
      }
    }

    // UPDATE POLLUTION
    printf("Update Pollution\n");
    if (calculatePollution == true &&
        (((float)((int)currentTime)) == (currentTime)) &&
        ((int)currentTime % ((int)60 * 6)) == 0) { //each 6 min
      cudaGridPollution.addValueToGrid(currentTime, trafficPersonVec, simRoadGraph, clientMain, laneMapNumToEdgeDesc);
    }

    currentTime += deltaTime;
    steps++;

    // RENDER IF NECESSARY
    printf("Render if necessary\n");
    if (clientMain!=nullptr && clientMain->ui.b18RenderSimulationCheckBox->isChecked()) {
      /*while (clientMain->ui.b18RenderStepSpinBox->value() == 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        clientMain->glWidget3D->updateGL();
        QApplication::processEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }*/

      if (steps % clientMain->ui.b18RenderStepSpinBox->value() == 0) { //each "steps" steps, render
        QString timeT;
        int timeH = currentTime / 3600.0f;
        int timeM = (currentTime - timeH * 3600.0f) / 60.0f;
        timeT.sprintf("%d:%02d", timeH, timeM);
        //clientMain->ui.timeLCD->display(timeT);
        printf("About to render\n");
        clientMain->glWidget3D->updateGL();
        QApplication::processEvents();
      }
    }
    printf("Finish one step\n");
  }

  if (DEBUG_SIMULATOR) {
    printf("<<End Simulation (count %d) TIME: %d ms.\n", count, timer.elapsed());
  }

  {
    // Total number of car steps
    uint totalNumSteps = 0;
    float totalGas = 0.0f;

    for (int p = 0; p < numPeople; p++) {
      totalNumSteps += trafficPersonVec[p].num_steps;
      totalGas += trafficPersonVec[p].gas;
    }

    avgTravelTime = (totalNumSteps) / (trafficPersonVec.size() * 60.0f); //in min
    printf("(Count %d) Total num steps %u Avg %f min Avg Gas %f. Calculated in %d ms\n",
           count, totalNumSteps, avgTravelTime, totalGas / trafficPersonVec.size(),
           timer.elapsed());
  }

  //clientMain->glWidget3D->updateGL();

  if (DEBUG_SIMULATOR) {
    printf("<<simulate\n");
  }
}//



void B18TrafficSimulator::render(VBORenderManager &rendManager) {
  ///////////////////////////////
  // RENDER POLUTION
  //printf("pol %d\n",clientMain->ui.pollutionCheckBox->isChecked());
  /*if(clientMain->ui.pollutionCheckBox->isChecked()){
        int indexToRead=LC::misctools::Global::global()->cuda_current_time_slider-70;
        cudaGridPollution.renderPollution(indexToRead);
  }*/
  rendManager.removeStaticGeometry("Sim_Points", false);
  rendManager.removeStaticGeometry("Sim_Box", false);
  /////////////////////////////
  // RENDER CARS
  const bool renderMultiEdge = false;
  const bool renderBoxes = false;
  const bool renderModels = false;

  // init render people
  if (b18TrafficSimulatorRender.size() != trafficPersonVec.size() && trafficPersonVec.size() > 0) {
    printf("SIm init render people");
    b18TrafficSimulatorRender.clear();
    b18TrafficSimulatorRender.resize(trafficPersonVec.size());
  }

  if (b18TrafficLightRender.size() != trafficLights.size() && trafficLights.size() > 0) {
    printf("SIm init traffic sim");
    b18TrafficLightRender.clear();
    b18TrafficLightRender.resize(trafficLights.size());
  }

  //////////////////////////
  // RENDER CAR AS MODELS
  /*if(renderModels==true&&trafficPersonVec.size()>0){
        int numPeople=trafficPersonVec.size();
        const float heightPoint=2.0f;
        //glDisable(GL_DEPTH_TEST);
        //glBegin(GL_QUADS);

        // SAVE TO FILE
        QTextStream out;
        QFile* file;
        int numPeopleActive=0;
        if(clientMain->ui.animationSaveCheckBox->isChecked()){
                QString fileName="data/animation.txt";
                file=new QFile(fileName);
                if(!file->open(QIODevice::WriteOnly|QIODevice::Append)) {
                        printf("ERROR: Not possible to read %s\n",fileName.toLatin1().constData());
                        exit(0);
                }
                out.setDevice(file);
                for(int p=0;p<numPeople;p++){
                        if(trafficPersonVec[p].active!=1){
                                continue;
                        }
                        numPeopleActive++;
                }
                out<<"NumCars:"<<numPeopleActive<<"\n";
        }


        for(int p=0;p<numPeople;p++){

                //0. check if finished
                if(trafficPersonVec[p].active!=1){
                        continue;
                }
                glColor3ub((trafficPersonVec[p].color>>24)&0xFF,(trafficPersonVec[p].color>>16)&0xFF,(trafficPersonVec[p].color>>8)&0xFF);
                //position
                if(laneMapNumToEdgeDesc.find(trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge])==laneMapNumToEdgeDesc.end()){
                        printf("ERROR\n");//edge not found in map
                        continue;
                }
                QVector3D p0,p1;
                float posInLaneM= trafficPersonVec[p].posInLaneM;
                float lenghtLane= trafficPersonVec[p].length;//meters?
                RoadGraph::roadGraphEdgeDesc_BI ei=laneMapNumToEdgeDesc[trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]];
                // MULTI EDGE
                if(renderMultiEdge==true&& simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size()>0){

                        for(int gN=0;gN<simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size()-1;gN++){//note -1
                                float currLe=(simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN]-simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN+1]).length();
                                if(posInLaneM>currLe){//next edge
                                        posInLaneM-=currLe;
                                        continue;
                                }
                                 p0 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN];
                                 p1 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN+1];
                                break;
                        }
                }else{
                        // ONE EDGE
                         p0 = simRoadGraph->myRoadGraph_BI[boost::source(ei, simRoadGraph->myRoadGraph_BI)].pt;
                         p1 = simRoadGraph->myRoadGraph_BI[boost::target(ei, simRoadGraph->myRoadGraph_BI)].pt;
                }
                QVector3D dir=(p1-p0).normalized();
                QVector3D per=(QVector3D::crossProduct(QVector3D(0,0,1.0f),dir).normalized());
                bool goodPoint=true;
                float range=-1.0f;
                if(
                        //posInLaneM<intersectionClearance+range||
                        posInLaneM<intersectionClearance+range*2.0f||
                        posInLaneM>(lenghtLane-intersectionClearance-range)){
                                goodPoint=false;
                }
                //printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
                float perShift=-0.5f*LC::misctools::Global::global()->roadLaneWidth*(1+2*trafficPersonVec[p].numOfLaneInEdge);
                QVector3D v=p0+dir*posInLaneM+perShift*per;// center of the back of the car
                ///////////////////////////////
                // INT
                //printf("Pos %f %f %f ",v.x(),v.y(),v.z());
                if(B18TrafficSimulatorRender[p].getInterpolated(goodPoint,v,dir,v,dir)==false){//copy and ref
                        // not render, just write a negative z
                        if(clientMain->ui.animationSaveCheckBox->isChecked()){//save to file
                                out<<"Car:"<<QString::number(p)<<":pos:"<<10000.0f<<":"<<10000.0f<<":-100.0f:dir:"<<dir.x()<<":"<<dir.y()<<":"<<dir.z()<<":\n";
                        }
                        continue;//no ready to render
                }
                //printf("Inter %f %f %f\n",v.x(),v.y(),v.z());
                ///////////////////////////////
                ModelSpec specR;
                specR.type=0;//car
                specR.transMatrix.setToIdentity();
                //specR.transMatrix.rotate(180.0f,0,0,1.0f);

                float angle;
                //angle=acos(QVector3D::dotProduct(QVector3D(1.0,0,0),dir))*57.2957795f+90.0f;
                {
                        QVector3D referenceForward =QVector3D(1.0,0,0);
                        QVector3D referenceRight= QVector3D::crossProduct(QVector3D(0,0,1.0f), referenceForward);
                        QVector3D newDirection = dir;
                        //float angle = Vector3.Angle(newDirection, referenceForward);
                        {
                                float s = QVector3D::crossProduct(referenceForward,newDirection).length();
                                float c = QVector3D::dotProduct(referenceForward,newDirection);
                                angle = atan2(s, c);
                        }
                        float sign = (QVector3D::dotProduct(newDirection, referenceRight) > 0.0f) ? 1.0f: -1.0f;
                        angle = sign * angle;
                        angle*=57.2957795f;//to defrees
                        angle+=90.0f;//model is original turned 90degrees
                }
                specR.transMatrix.translate(v.x(),v.y(),v.z());//-1.5f is the street level
                specR.transMatrix.rotate(angle,0,0,1.0f);

                clientMain->mGLWidget_3D->modelsVBO.renderModel(specR);
                if(clientMain->ui.animationSaveCheckBox->isChecked()){//save to file
                        out<<"Car:"<<QString::number(p)<<":pos:"<<v.x()<<":"<<v.y()<<":"<<v.z()<<":dir:"<<dir.x()<<":"<<dir.y()<<":"<<dir.z()<<":\n";
                }

        }//for people

        if(clientMain->ui.animationSaveCheckBox->isChecked()){
                out.flush();
                file->close();
                //file.close();
        }
        //glEnd();
        //glEnable(GL_DEPTH_TEST);
  }//model*/

  //////////////////////////
  // RENDER CAR AS BOXES
  if (renderBoxes == true && trafficPersonVec.size() > 0) {
    std::vector<Vertex> carQuads;
    QVector3D carColor;
    int numPeople = trafficPersonVec.size();
    //carQuads.resize(4*numPeople);//quad per car many not active

    const float heightPoint = 5.0f;

    for (int p = 0; p < numPeople; p++) {

      //0. check if finished
      if (trafficPersonVec[p].active != 1) {
        continue;
      }

      carColor = QVector3D(1, 0, 0);

      //position
      if (laneMapNumToEdgeDesc.find(
            trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]) ==
          laneMapNumToEdgeDesc.end()) {
        printf("ERROR\n");//edge not found in map
        continue;
      }

      QVector3D p0, p1;
      float posInLaneM = trafficPersonVec[p].posInLaneM;
      RoadGraph::roadGraphEdgeDesc_BI ei =
        laneMapNumToEdgeDesc[trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]];

      // MULTI EDGE
      if (renderMultiEdge == true &&
          simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size() > 0) {

        for (int gN = 0;
             gN < simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size() - 1;
             gN++) { //note -1
          float currLe = (simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN] -
                          simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN + 1]).length();

          if (posInLaneM > currLe) { //next edge
            posInLaneM -= currLe;
            continue;
          }

          p0 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN];
          p1 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN + 1];
          break;
        }
      } else {
        // ONE EDGE
        p0 = simRoadGraph->myRoadGraph_BI[boost::source(ei,
                                          simRoadGraph->myRoadGraph_BI)].pt;
        p1 = simRoadGraph->myRoadGraph_BI[boost::target(ei,
                                          simRoadGraph->myRoadGraph_BI)].pt;
      }

      QVector3D dir = (p1 - p0).normalized();
      QVector3D per = (QVector3D::crossProduct(QVector3D(0, 0, 1.0f),
                       dir).normalized());
      //printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
      float perShift = -0.5f * G::global().getFloat("roadLaneWidth") *
                       (1 + 2 * trafficPersonVec[p].numOfLaneInEdge);
      QVector3D v = p0 + dir * posInLaneM + perShift *
                    per; // center of the back of the car
      float widthHalf = 1.82f / 2.0f;
      float length = 4.12f;
      QVector3D v0 = v - widthHalf * per;
      QVector3D v1 = v - widthHalf * per + length * dir;
      QVector3D v2 = v + widthHalf * per + length * dir;
      QVector3D v3 = v + widthHalf * per;



      carQuads.push_back(Vertex(QVector3D(v0.x(), v0.y(), heightPoint), carColor,
                                QVector3D(), QVector3D()));
      carQuads.push_back(Vertex(QVector3D(v1.x(), v1.y(), heightPoint), carColor,
                                QVector3D(), QVector3D()));
      carQuads.push_back(Vertex(QVector3D(v2.x(), v2.y(), heightPoint), carColor,
                                QVector3D(), QVector3D()));
      carQuads.push_back(Vertex(QVector3D(v3.x(), v3.y(), heightPoint), carColor,
                                QVector3D(), QVector3D()));

    }//for people

    //printf("carQuads %d\n",carQuads.size());

    if (carQuads.size() > 0) {
      rendManager.addStaticGeometry("Sim_Box", carQuads, "", GL_QUADS,
                                    1 | mode_AdaptTerrain);
      glDisable(GL_DEPTH_TEST);
      rendManager.renderStaticGeometry("Sim_Box");
      glEnable(GL_DEPTH_TEST);
    }
  }//boxes


  //////////////////////////
  // RENDER CAR AS POINTS
  if (renderBoxes == false && renderModels == false && trafficPersonVec.size() > 0) {
    // printf("SIm render as points");
    //rendManager.removeStaticGeometry("Sim_Points",false);
    const float heightPoint = 5.0f;
    std::vector<Vertex> carPoints;
    QVector3D pointColor;
    int numPeople = trafficPersonVec.size();
    //carPoints.resize(numPeople);//many not active
    int activeCars = 0;

    for (int p = 0; p < numPeople; p++) {
      //printf(" 0. Person: %d\n",p);
      //0. check if finished
      if (trafficPersonVec[p].active == 1) {
        //printf("%x\n",trafficPersonVec[p].color);
        pointColor = QVector3D((trafficPersonVec[p].color >> 24) & 0xFF,
                               (trafficPersonVec[p].color >> 16) & 0xFF,
                               (trafficPersonVec[p].color >> 8) & 0xFF) / 255.0f;

        //position
        if (laneMapNumToEdgeDesc.find(trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]) == laneMapNumToEdgeDesc.end()) {
          printf("ERROR\n");//edge not found in map
          continue;
        }

        RoadGraph::roadGraphEdgeDesc_BI ei =
          laneMapNumToEdgeDesc[trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]];

        // MULTI EDGE
        if (renderMultiEdge == true &&
            simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size() > 0) {
          float posInLaneM = trafficPersonVec[p].posInLaneM;

          for (int gN = 0;
               gN < simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry.size() - 1;
               gN++) { //note -1
            float currLe = (simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN] -
                            simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN + 1]).length();

            if (posInLaneM > currLe) { //next edge
              posInLaneM -= currLe;
              continue;
            }

            QVector3D p0 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN];
            QVector3D p1 = simRoadGraph->myRoadGraph_BI[ei].roadSegmentGeometry[gN + 1];
            QVector3D dir = (p1 - p0).normalized();
            QVector3D per = (QVector3D::crossProduct(QVector3D(0, 0, 1.0f),
                             dir).normalized());
            //printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
            float perShift = -0.5f * G::global().getFloat("roadLaneWidth") *
                             (1 + 2 * trafficPersonVec[p].numOfLaneInEdge);
            QVector3D v = p0 + dir * posInLaneM + perShift * per;
            //glVertex3f(v.x(),v.y(),heightPoint);
            carPoints[p] = Vertex(QVector3D(v.x(), v.y(), heightPoint), pointColor,
                                  QVector3D(), QVector3D());
            activeCars++;
            break;
          }
        } else {
          // ONE EDGE
          QVector3D p0 = simRoadGraph->myRoadGraph_BI[boost::source(ei,
                         simRoadGraph->myRoadGraph_BI)].pt;
          QVector3D p1 = simRoadGraph->myRoadGraph_BI[boost::target(ei,
                         simRoadGraph->myRoadGraph_BI)].pt;
          QVector3D dir = (p1 - p0).normalized();
          QVector3D per = (QVector3D::crossProduct(QVector3D(0, 0, 1.0f),
                           dir).normalized());
          //printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
          float perShift = -0.5f * G::global().getFloat("roadLaneWidth") *
                           (1 + 2 * trafficPersonVec[p].numOfLaneInEdge);
          QVector3D v = p0 + dir * trafficPersonVec[p].posInLaneM + perShift * per;
          //glVertex3f(v.x(),v.y(),heightPoint);
          carPoints.push_back(Vertex(QVector3D(v.x(), v.y(), heightPoint), pointColor,
                                     QVector3D(), QVector3D()));
          activeCars++;
        }
      }
    }

    //glEnd();
    if (carPoints.size() > 0) {
      rendManager.addStaticGeometry("Sim_Points", carPoints, "", GL_POINTS, 1 | mode_AdaptTerrain);
      glDisable(GL_DEPTH_TEST);
      glPointSize(15.0f);
      rendManager.renderStaticGeometry("Sim_Points");
      glEnable(GL_DEPTH_TEST);
    }
  }

  /////////////////////////////
  // RENDER INTERSECTIONS
  RoadGraph::roadGraphEdgeIter_BI ei, eiEnd;
  QVector3D p0;
  QVector3D p1;

  //Render edges ===============
  if (false && edgeDescToLaneMapNum.size() > 0) {
    int numEdges = 0;

    for (boost::tie(ei, eiEnd) = boost::edges(simRoadGraph->myRoadGraph_BI);
         ei != eiEnd; ++ei) {


      p0 = simRoadGraph->myRoadGraph_BI[boost::source(*ei,
                                        simRoadGraph->myRoadGraph_BI)].pt;
      p1 = simRoadGraph->myRoadGraph_BI[boost::target(*ei,
                                        simRoadGraph->myRoadGraph_BI)].pt;
      /// N
      QVector3D dir = (p1 - p0).normalized();
      QVector3D per = (QVector3D::crossProduct(QVector3D(0, 0, 1.0f),
                       dir).normalized());

      ///
      int numL = simRoadGraph->myRoadGraph_BI[*ei].numberOfLanes;
      uint laneMapNum = edgeDescToLaneMapNum[*ei];

      for (int lN = 0; lN < numL; lN++) {
        uchar trafficLight = 0x01; //default

        //uint lN=0;
        if (trafficLights.size() > 0) {
          trafficLight = trafficLights[laneMapNum + lN];
        }

        b18TrafficLightRender[laneMapNum + lN].getInterpolated(trafficLight,
            trafficLight);

        switch (trafficLight) {
        case 0x00:
          glColor3ub(255, 0, 0);
          break;

        case 0xFF:
          glColor3ub(0, 255, 0);
          break;

        default:
          glColor3ub(0, 0, 0);
          //printf("Unknown traffic light state\n");
        }

        float perShift = -0.5f * G::global().getFloat("roadLaneWidth") * (1 + 2 * lN);
        QVector3D p1S = QVector3D(p1.x() + per.x() * perShift,
                                  p1.y() + per.y() * perShift, p1.z());
        p1S -= intersectionClearance * dir; //move to be brefore clearance
        /// CONE
        //if(LC::misctools::Global::global()->view_arterials_arrows_render==true){
        LC::misctools::drawCone(dir, p1S - dir * 0.5f, 2.0f, 1.0f, 16);
        //}
      }

      //
    }
  }
}//

void B18TrafficSimulator::calculateAndDisplayTrafficDensity(/*RoadGraph& inRoadGraph,*/std::vector<float>
    &accSpeedPerLinePerTimeInterval,
    std::vector<float> &numVehPerLinePerTimeInterval,
    std::map<RoadGraph::roadGraphEdgeDesc_BI, uint> &edgeDescToLaneMapNum,
    std::map<uint, RoadGraph::roadGraphEdgeDesc_BI> &laneMapNumToEdgeDesc,
    int tNumLanes) {

  const float numStepsTogether = 12;
  int numSampling = accSpeedPerLinePerTimeInterval.size() / tNumLanes;

  if (DEBUG_SIMULATOR) {
    printf(">>calculateAndDisplayTrafficDensity numSampling %d\n", numSampling);
  }

  RoadGraph::roadGraphEdgeIter_BI ei, eiEnd;

  //std::vector<int> numTotalCarsInLanes;
  //numTotalCarsInLanes.resize(numSampling);
  //memset(numTotalCarsInLanes.data(),0,numTotalCarsInLanes.size()*sizeof(int
  printPercentageMemoryUsed();
  if (DEBUG_SIMULATOR) {
    printf(">>calculateAndDisplayTrafficDensity Allocate memory numSampling %d\n", numSampling);
  }
  
  int count = 0;
  for (boost::tie(ei, eiEnd) = boost::edges(simRoadGraph->myRoadGraph_BI); ei != eiEnd; ++ei) {
    if (edgeDescToLaneMapNum.count(*ei) == 0) {
      continue;
    }
    printf("Edge allocation %d -> lane %d\n", count, edgeDescToLaneMapNum[*ei]);
    printPercentageMemoryUsed();
    simRoadGraph->myRoadGraph_BI[*ei].averageSpeed.resize(numSampling);
    simRoadGraph->myRoadGraph_BI[*ei].averageUtilization.resize(numSampling);
    if ((++count) % 200 == 0) {
      //printf("Edge allocation %d: ", count);
      
    }
  }
  printPercentageMemoryUsed();

  if (DEBUG_SIMULATOR) {
    printf(">>calculateAndDisplayTrafficDensity Process\n");
  }
  for (boost::tie(ei, eiEnd) = boost::edges(simRoadGraph->myRoadGraph_BI); ei != eiEnd; ++ei) {
    int numLanes = simRoadGraph->myRoadGraph_BI[*ei].numberOfLanes;

    if (numLanes == 0) {
      continue;//edges with zero lines just skip
    }
    int numLane = edgeDescToLaneMapNum[*ei];
    int offset;
    //0.8f to make easier to become red
    float maxVehicles = 0.5f * simRoadGraph->myRoadGraph_BI[*ei].edgeLength * simRoadGraph->myRoadGraph_BI[*ei].numberOfLanes / (s_0);

    if (maxVehicles < 1.0f) {
      maxVehicles = 1.0f;
    }

    for (int sa = 0; sa < numSampling - 1; sa++) {
      offset = sa * tNumLanes;

      ////////////////////////////////////////////////
      // avarage speed
      /*if(numVehPerLinePerTimeInterval[numLane+offset]>0){
        simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]=accSpeedPerLinePerTimeInterval[numLane+offset]/(numVehPerLinePerTimeInterval[numLane+offset]*numStepsTogether);
        if(simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]>1)
                simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]=1.0f;
      }else{
        simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]=simRoadGraph->myRoadGraph_BI[*ei].maxSpeedMperSec;
      }*/
      if (numVehPerLinePerTimeInterval[numLane + offset]*numStepsTogether > 0) {
        simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa] = (accSpeedPerLinePerTimeInterval[numLane + offset]) / ((float) numVehPerLinePerTimeInterval[numLane + offset]); //!!!!!
      } else {
        simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa] = 0;
      }

      //if(simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]!=0)
      //	printf("num %f acc %f -->avfg %f\n",numVehPerLinePerTimeInterval[numLane+offset]/numStepsTogether,accSpeedPerLinePerTimeInterval[numLane+offset]/numStepsTogether,simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]);
      //simRoadGraph->myRoadGraph_BI[*ei].averageSpeed[sa]=numVehPerLinePerTimeInterval[numLane+offset]/numStepsTogether;//!!!!!
      ///////////////////////////////
      // average utilization
      simRoadGraph->myRoadGraph_BI[*ei].averageUtilization[sa] = numVehPerLinePerTimeInterval[numLane + offset] / (maxVehicles * numStepsTogether);

      if (simRoadGraph->myRoadGraph_BI[*ei].averageUtilization[sa] > 1.0f) {
        //printf("numVehPerLinePerTimeInterval[numLane+offset] %d maxVe %f--> %f\n",numVehPerLinePerTimeInterval[numLane+offset],maxVehicles,simRoadGraph->myRoadGraph_BI[*ei].averageUtilization[sa]);
        simRoadGraph->myRoadGraph_BI[*ei].averageUtilization[sa] = 1;
      }
    }
  }

  /*for(int sa=0;sa<numSampling-1;sa++){
        printf("%d ",numTotalCarsInLanes[sa]);
  }*/
  if (DEBUG_SIMULATOR) {
    printf("\n");
  }

  //printf("\n<<CUDATrafficDesigner::calculateAndDisplayTrafficDensity numTotalCarsInLanes %d\n",numTotalCarsInLanes);
}//

//////////////////////////////
// control the number of samples
static const uint numElements = 200;

bool B18TrafficSimulatorRender::getInterpolated(bool goodPoint,
    QVector3D newPos, QVector3D newDir, QVector3D &interPos, QVector3D &interDir) {
  // initialization
  if (positions.size() < numElements) {
    positions.push_back(newPos);
    directions.push_back(newDir);
    goodPoints.push_back(true);
    interPos = positions[0];
    interDir = directions[0];
    indexToRead = 0;
    return false;
  }

  interPos = positions[indexToRead];
  interDir = directions[indexToRead];

  if (goodPoints[indexToRead] == false) {
    // reconstruct!!!
    // 1. find the following good point and length
    QVector3D p_0 = positions[indexToRead];
    QVector3D d_0 = directions[indexToRead];
    QVector3D p_1, d_1;
    float intersectionLength = 0;
    bool foundCorrect = false;
    int index;
    std::vector<float> lengthToStart;

    for (int i = 1; i < numElements - 1; i++) { //note i=1
      index = (indexToRead + i) % numElements;
      int formerIndex = index - 1;

      if (formerIndex < 0) {
        formerIndex = numElements - 1;
      }

      intersectionLength += (positions[index] - positions[formerIndex]).length();
      lengthToStart.push_back(intersectionLength);

      if (goodPoints[index] == true) {
        foundCorrect = true;
        break;
      }

    }

    /*if(foundCorrect==false){
        printf("No enough elements to reconstruct!!!\n");
        exit(0);
    }*/
    p_1 = positions[index];
    d_1 = directions[index];

    // 2. update values with curve
    for (int i = (indexToRead + 1) % numElements, j = 1; i != index;
         i = (i + 1) % numElements, j++) { //skip first because it stays
      float currenLength = lengthToStart[j - 1];
      float t = currenLength / intersectionLength; //0-1
      // en.wikipedia.org/wiki/Hermite_curve
      // {p}(t) = (2t^3-3t^2+1){p}_0 + (t^3-2t^2+t){m}_0 + (-2t^3+3t^2){p}_1 +(t^3-t^2){m}_1
      float t2 = t * t;
      float t3 = t2 * t;
      float d_0_factor = 2.0f;
      positions[i] = (2 * (t3) - 3 * (t2) + 1) * p_0 + ((t3) - 2 *
                     (t2) + t) * d_0 * d_0_factor + (-2 * (t3) + 3 * (t2)) * p_1 + ((t3) -
                         (t2)) * d_1;

      // dir
      //int beforeIndex=i-1;
      //if(beforeIndex<0)beforeIndex=numElements-1;
      //directions[i]=(positions[i]-positions[beforeIndex]).normalized();
      t = t - 0.1f;

      if (t < 0) {
        t = 0;
      }

      t2 = t * t;
      t3 = t2 * t;
      QVector3D beforePos = (2 * (t3) - 3 * (t2) + 1) * p_0 + ((t3) - 2 *
                            (t2) + t) * d_0 * d_0_factor + (-2 * (t3) + 3 * (t2)) * p_1 + ((t3) -
                                (t2)) * d_1;
      directions[i] = (positions[i] - beforePos).normalized();
      // !! direction
      goodPoints[i] = true;
    }
  }

  positions[indexToRead] = newPos;
  directions[indexToRead] = newDir;
  goodPoints[indexToRead] = goodPoint;
  indexToRead = (indexToRead + 1) % numElements;
  return true;
  //printf("Pos %f %f %f Inter %f %f %f\n",newPos.x(),newPos.y(),newPos.z(),interPos.x(),interPos.y(),interPos.z());
}//


void B18TrafficLightRender::getInterpolated(uchar newTrafficLight,
    uchar &interTrafficLight) {
  if (trafficLight.size() < numElements) {
    trafficLight.push_back(newTrafficLight);
    interTrafficLight = trafficLight[0];
    indexToRead = 0;
    return;
  }

  interTrafficLight = trafficLight[indexToRead];
  trafficLight[indexToRead] = newTrafficLight;
  indexToRead = (indexToRead + 1) % numElements;
}//

}
