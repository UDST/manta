//CUDA CODE
#include <stdio.h>
#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"

#include "b18TrafficPerson.h"
#include "b18EdgeData.h"
#include <vector>

#define ushort unsigned short
#define uint unsigned int
#define uchar unsigned char

///////////////////////////////
// CONSTANTS

__constant__ float intersectionClearance = 7.8f;
__constant__ float s_0 = 7.0f;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

////////////////////////////////
// VARIABLES
LC::B18TrafficPerson *trafficPersonVec_d;
uint *indexPathVec_d;
LC::B18EdgeData *edgesData_d;

__constant__ bool calculatePollution = true;
__constant__ float cellSize = 1.0f;

__constant__ float deltaTime = 0.5f;
const float deltaTimeH = 0.5f;

const uint numStepsPerSample = 30.0f / deltaTimeH; //each min
const uint numStepsTogether = 12; //change also in density (10 per hour)

uchar *laneMap_d;
bool readFirstMapC=true;
uint mapToReadShift;
uint mapToWriteShift;
uint halfLaneMap;
float startTime;


LC::B18IntersectionData *intersections_d;
uchar *trafficLights_d;

float* accSpeedPerLinePerTimeInterval_d;
float* numVehPerLinePerTimeInterval_d;

void b18InitCUDA(
  bool fistInitialization,
  std::vector<LC::B18TrafficPerson>& trafficPersonVec, 
  std::vector<uint> &indexPathVec, 
  std::vector<LC::B18EdgeData>& edgesData, 
  std::vector<uchar>& laneMap, 
  std::vector<uchar>& trafficLights, 
  std::vector<LC::B18IntersectionData>& intersections,
  float startTimeH, float endTimeH,
  std::vector<float>& accSpeedPerLinePerTimeInterval,
  std::vector<float>& numVehPerLinePerTimeInterval) {
  printf(">>b18InitCUDA\n");

  { // people
    size_t size = trafficPersonVec.size() * sizeof(LC::B18TrafficPerson);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &trafficPersonVec_d, size));   // Allocate array on device
    gpuErrchk(cudaMemcpy(trafficPersonVec_d, trafficPersonVec.data(), size, cudaMemcpyHostToDevice));
  }
  
  { // indexPathVec
    size_t sizeIn = indexPathVec.size() * sizeof(uint);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &indexPathVec_d, sizeIn));   // Allocate array on device
    gpuErrchk(cudaMemcpy(indexPathVec_d, indexPathVec.data(), sizeIn, cudaMemcpyHostToDevice));
  }
  {//edgeData
    size_t sizeD = edgesData.size() * sizeof(LC::B18EdgeData);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &edgesData_d, sizeD));   // Allocate array on device
    gpuErrchk(cudaMemcpy(edgesData_d, edgesData.data(), sizeD, cudaMemcpyHostToDevice));
  }
  {//laneMap
    size_t sizeL = laneMap.size() * sizeof(uchar);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &laneMap_d, sizeL));   // Allocate array on device
    gpuErrchk(cudaMemcpy(laneMap_d, laneMap.data(), sizeL, cudaMemcpyHostToDevice));
    halfLaneMap = laneMap.size() / 2;
  }
  {// intersections
    size_t sizeI = intersections.size() * sizeof(LC::B18IntersectionData);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &intersections_d, sizeI));   // Allocate array on device
    gpuErrchk(cudaMemcpy(intersections_d, intersections.data(), sizeI, cudaMemcpyHostToDevice));
    size_t sizeT = trafficLights.size() * sizeof(uchar);//total number of lanes
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &trafficLights_d, sizeT));   // Allocate array on device
    gpuErrchk(cudaMemcpy(trafficLights_d, trafficLights.data(), sizeT, cudaMemcpyHostToDevice));
  }
  {
    startTime = startTimeH * 3600.0f;
    uint numSamples = ceil(((endTimeH*3600.0f - startTimeH*3600.0f) / (deltaTimeH * numStepsPerSample * numStepsTogether))) + 1; //!!!
    accSpeedPerLinePerTimeInterval.clear();
    numVehPerLinePerTimeInterval.clear();
    accSpeedPerLinePerTimeInterval.resize(numSamples * trafficLights.size());
    numVehPerLinePerTimeInterval.resize(numSamples * trafficLights.size());
    size_t sizeAcc = accSpeedPerLinePerTimeInterval.size() * sizeof(float);
    if (fistInitialization)gpuErrchk(cudaMalloc((void **) &accSpeedPerLinePerTimeInterval_d, sizeAcc));   // Allocate array on device
    if (fistInitialization)gpuErrchk(cudaMalloc((void **) &numVehPerLinePerTimeInterval_d, sizeAcc));   // Allocate array on device
    gpuErrchk(cudaMemset(&accSpeedPerLinePerTimeInterval_d[0], 0, sizeAcc));
    gpuErrchk(cudaMemset(&numVehPerLinePerTimeInterval_d[0], 0, sizeAcc));
  }
  printf("<<b18InitCUDA\n");
}//

void b18FinishCUDA(void){
	//////////////////////////////
	// FINISH
	cudaFree(trafficPersonVec_d);
  cudaFree(indexPathVec_d);
	cudaFree(edgesData_d);
	cudaFree(laneMap_d);
	cudaFree(intersections_d);
	cudaFree(trafficLights_d);

  cudaFree(accSpeedPerLinePerTimeInterval_d);
  cudaFree(numVehPerLinePerTimeInterval_d);
}//

 void b18GetDataCUDA(std::vector<LC::B18TrafficPerson>& trafficPersonVec){
	 // copy back people
	 size_t size = trafficPersonVec.size() * sizeof(LC::B18TrafficPerson);
	 cudaMemcpy(trafficPersonVec.data(),trafficPersonVec_d,size,cudaMemcpyDeviceToHost);//cudaMemcpyHostToDevice
 }


 __device__ void calculateGapsLC(
   uint mapToReadShift,
   uchar* laneMap,
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
   for (ushort b = initShift - 1; (b < numOfCells) && (found == false); b++) { //NOTE -1 to make sure there is none in at the same level
     // laneChar = laneMap[mapToReadShift + maxWidth * (laneToCheck) + b];
     const uint posToSample = mapToReadShift + kMaxMapWidthM * (laneToCheck + (((int) (b / kMaxMapWidthM)) * numLinesEdge)) + b % kMaxMapWidthM;
     laneChar = laneMap[posToSample];

     if (laneChar != 0xFF) {
       gap_a = ((float) b - initShift); //m
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
       gap_b = ((float) initShift - b); //m
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

 __device__ void calculateLaneCarShouldBe(
   uint curEdgeLane,
   uint nextEdge,
   LC::B18IntersectionData* intersections,
   uint edgeNextInters,
   ushort edgeNumLanes,
   ushort &initOKLanes,
   ushort &endOKLanes) {
   initOKLanes = 0;
   endOKLanes = edgeNumLanes;
   bool currentEdgeFound = false;
   bool exitFound = false;
   ushort numExitToTake = 0;
   ushort numExists = 0;

   for (int eN = intersections[edgeNextInters].totalInOutEdges - 1; eN >= 0; eN--) {  // clockwise
     uint procEdge = intersections[edgeNextInters].edge[eN];

     if ((procEdge & kMaskLaneMap) == curEdgeLane) { //current edge 0xFFFFF
       currentEdgeFound = true;
       if (exitFound == false) {
         numExitToTake = 0;
       }
       continue;
     }

     if ((procEdge & kMaskInEdge) == 0x0) { //out edge 0x800000
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
     }
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

 // Kernel that executes on the CUDA device
__global__ void kernel_trafficSimulation(
	 int numPeople,
	 float currentTime,
   uint mapToReadShift,
   uint mapToWriteShift,
	 LC::B18TrafficPerson *trafficPersonVec,
   uint *indexPathVec,
   LC::B18EdgeData* edgesData,
	 uchar *laneMap,
   LC::B18IntersectionData *intersections,
	 uchar *trafficLights
	 )
 {
	 const bool DEBUG_TRAFFIC=0;
	 int p = blockIdx.x * blockDim.x + threadIdx.x;
	 //printf("p %d Numpe %d\n",p,numPeople);
   if (p < numPeople) {//CUDA check (inside margins)
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
         trafficPersonVec[p].indexPathCurr = trafficPersonVec[p].indexPathInit; // reset index.
         uint firstEdge = indexPathVec[trafficPersonVec[p].indexPathCurr];

         if (firstEdge == -1) {
           trafficPersonVec[p].active = 2;
           //printf("0xFFFF\n");
           return;
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
         ushort initShift = (ushort) (0.5f * numOfCells); //number of cells it should be placed (half of road)

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
           uchar vInMpS = (uchar) (trafficPersonVec[p].v *
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
         uint nextEdge = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];

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
     uint currentEdge = indexPathVec[trafficPersonVec[p].indexPathCurr];
     uint nextEdge = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];
     if (DEBUG_TRAFFIC == 1) {
       printf("    2. Person: %d Try to move current Edge %u Next %u\n", p, currentEdge, nextEdge);
     }
     // www.vwi.tu-dresden.de/~treiber/MicroApplet/IDM.html
     // IDM
     float thirdTerm = 0;
     ///////////////////////////////////////////////////
     // 2.1.1 Find front car
     int numCellsCheck = fmaxf(30.0f, trafficPersonVec[p].v * deltaTime * 2); //30 or double of the speed*time

     // a) SAME LINE (BEFORE SIGNALING)
     bool found = false;
     bool noFirstInLaneBeforeSign = false; //use for stop control (just let 1st to pass)
     bool noFirstInLaneAfterSign = false; //use for stop control (just let 1st to pass)
     float s;
     float delta_v;
     uchar laneChar;
     ushort byteInLine = (ushort) floor(trafficPersonVec[p].posInLaneM);
     ushort numOfCells = ceil((trafficPersonVec[p].length - intersectionClearance));

     for (ushort b = byteInLine + 2; (b < numOfCells) && (found == false) && (numCellsCheck > 0); b++, numCellsCheck--) {
       // laneChar = laneMap[mapToReadShift + maxWidth * ((indexPathVec[trafficPersonVec[p].indexPathCurr + trafficPersonVec[p].numOfLaneInEdge) + b];
       // ShiftRead + WIDTH * (width number * # edges + # laneInEdge) + b 
       const uint posToSample = mapToReadShift + kMaxMapWidthM * (indexPathVec[trafficPersonVec[p].indexPathCurr] + (((int) (byteInLine / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + b % kMaxMapWidthM;
       laneChar = laneMap[posToSample];

       if (laneChar != 0xFF) {
         s = ((float) (b - byteInLine)); //m
         delta_v = trafficPersonVec[p].v - (laneChar / 3.0f); //laneChar is in 3*ms (to save space in array)
         found = true;
         noFirstInLaneBeforeSign = true;
         break;
       }
     }

     // b) TRAFFIC LIGHT
     if (byteInLine < numOfCells && found == false && numCellsCheck > 0) { //before traffic signaling (and not cell limited)
       if (trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] == 0x00) { //red
         s = ((float) (numOfCells - byteInLine)); //m
         delta_v = trafficPersonVec[p].v - 0; //it should be treated as an obstacle
         nextVehicleIsATrafficLight = true;
         //printf("\nFOUND TL\n",s,delta_v);
         found = true;
       }
     }

     // c) SAME LINE (AFTER SIGNALING)
     for (ushort b = byteInLine + 2; (b < numOfCells) && (found == false) && (numCellsCheck > 0); b++, numCellsCheck--) {
       // laneChar = laneMap[mapToReadShift + maxWidth * ((indexPathVec[trafficPersonVec[p].indexPathCurr + trafficPersonVec[p].numOfLaneInEdge) + b];
       const uint posToSample = mapToReadShift + kMaxMapWidthM * (indexPathVec[trafficPersonVec[p].indexPathCurr] + (((int) (byteInLine / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + b % kMaxMapWidthM;
       laneChar = laneMap[posToSample];

       if (laneChar != 0xFF) {
         s = ((float) (b - byteInLine)); //m
         delta_v = trafficPersonVec[p].v - (laneChar /
           3.0f); //laneChar is in 3*ms (to save space in array)
         found = true;
         noFirstInLaneAfterSign = true;
         break;
       }
     }

     if (trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] == 0x0F && numCellsCheck > 0) { //stop 
       //check
       if (noFirstInLaneBeforeSign == false && byteInLine < numOfCells && //first before traffic
         trafficPersonVec[p].v == 0 && //stopped
         noFirstInLaneAfterSign == false) { // noone after the traffic light (otherwise wait before stop) !! TODO also check the beginning of next edge

         trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] = 0x00; //reset stop
         trafficPersonVec[p].posInLaneM = ceilf(numOfCells) + 1; //move magicly after stop

         if (DEBUG_TRAFFIC == 1) {
           printf("    2. Person: %d move after stop\n", p);
         }

       } else { //stop before STOP
         if (noFirstInLaneBeforeSign == false) { //just update this if it was the first one before sign
           s = ((float) (numOfCells - byteInLine)); //m
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
             s = ((float) (b)); //m
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
       s_star = s_0 + fmaxf(0.0f,
         (trafficPersonVec[p].v * trafficPersonVec[p].T + (trafficPersonVec[p].v *
         delta_v) / (2 * std::sqrt(trafficPersonVec[p].a * trafficPersonVec[p].b))));
       thirdTerm = std::pow(((s_star) / (s)), 2);
       //printf(">FOUND s_star %f thirdTerm %f!!!!\n",s_star,thirdTerm);
     }

     float dv_dt = trafficPersonVec[p].a * (1.0f - std::pow((
       trafficPersonVec[p].v / trafficPersonVec[p].maxSpeedMperSec), 4) - thirdTerm);

     // 2.1.3 update values
     numMToMove = fmaxf(0.0f,
       trafficPersonVec[p].v * deltaTime + 0.5f * (dv_dt) * deltaTime * deltaTime);

     if (DEBUG_TRAFFIC == 1) {
       printf("v %f v0 %f a %f dv_dt %f MOVE %f\n", trafficPersonVec[p].v,
         trafficPersonVec[p].maxSpeedMperSec, trafficPersonVec[p].a, dv_dt, numMToMove);
     }

     trafficPersonVec[p].v += dv_dt * deltaTime;

     if (trafficPersonVec[p].v < 0) {
       trafficPersonVec[p].v = 0;
     }

     //////////////////////////////////////////////

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

     if (trafficPersonVec[p].v == 0) { //if not moving not do anything else
       ushort posInLineCells = (ushort) (trafficPersonVec[p].posInLaneM);
       //laneMap[mapToWriteShift + maxWidth * (currentEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = 0;
       const uint posToSample = mapToWriteShift + kMaxMapWidthM * (currentEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;
       laneMap[posToSample] = 0;

       return;
     }


     ////////////////////////////////

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
               if ((((int) (x * 100) % 100) / 100.0f) < probabiltyMandatoryState) { // pseudo random number
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
                 if (int(trafficPersonVec[p].v) % 2 == 0) { // pseudo random
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
                   g_na_D = fmaxf(s_0, s_0 + b1A * trafficPersonVec[p].v + b2A *
                     (trafficPersonVec[p].v - v_a * 3.0f));

                   if (gap_a < g_na_D) { //gap smaller than critical gap
                     acceptLC = false;
                   }
                 }

                 if (acceptLC == true && gap_b != 1000.0f) {
                   g_bn_D = fmaxf(s_0, s_0 + b1B * v_b * 3.0f + b2B * (v_b * 3.0f -
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
                 printf("CATASTROFIC ERROR mandatory state");
               }
             }

             //printf(">>LANE CHANGE\n");
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
                 printf("CATASTROFIC ERROR mandatory lc");
               }

               //printf("M L %d R %d nL %u\n",leftLane,rightLane,trafficPersonVec[p].numOfLaneInEdge);
             }

             if (leftLane == true || rightLane == true) {

               // choose lane (if necessary)
               if (leftLane == true && rightLane == true) {
                 if ((int) (trafficPersonVec[p].posInLaneM) % 2 == 0) { //pseudo random
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
                   g_na_M = fmaxf(s_0, s_0 + (b1A * trafficPersonVec[p].v + b2A *
                     (trafficPersonVec[p].v - v_a * 3.0f)));

                   if (gap_a < g_na_M) { //gap smaller than critical gap
                     acceptLC = false;
                   }
                 }

                 if (acceptLC == true && gap_b != 1000.0f) {
                   g_bn_M = fmaxf(s_0, s_0 + (b1B * v_b * 3.0f + b2B * (v_b * 3.0f -
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

       uchar vInMpS = (uchar) (trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
       ushort posInLineCells = (ushort) (trafficPersonVec[p].posInLaneM);
       //laneMap[mapToWriteShift + maxWidth * (currentEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = vInMpS;
       const uint posToSample = mapToWriteShift + kMaxMapWidthM * (currentEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;
       laneMap[posToSample] = vInMpS;
       //printf("2<<LANE CHANGE\n");
       return;
     }

     //}
     //2.2 close to intersection
     if (nextEdge == -1) { //if(curr_intersection==end_intersection){
       trafficPersonVec[p].active = 2; //finished
       return;
     }

     trafficPersonVec[p].indexPathCurr++;
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
     uint nextNEdge = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];

     //trafficPersonVec[p].nextEdge=nextEdge;
     if (nextNEdge != -1) {
       //trafficPersonVec[p].nextPathEdge++;
       trafficPersonVec[p].LC_initOKLanes = 0xFF;
       trafficPersonVec[p].LC_endOKLanes = 0xFF;

       //2.2.3 update person edgeData
       //trafficPersonVec[p].nextEdge=nextEdge;
       trafficPersonVec[p].nextEdgemaxSpeedMperSec =
         edgesData[nextNEdge].maxSpeedMperSec;
       trafficPersonVec[p].nextEdgeNumLanes = edgesData[nextNEdge].numLines;
       trafficPersonVec[p].nextEdgeNextInters = edgesData[nextNEdge].nextInters;
       trafficPersonVec[p].nextEdgeLength = edgesData[nextNEdge].length;
     }

     trafficPersonVec[p].LC_stateofLaneChanging = 0;
     uchar vInMpS = (uchar) (trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
     ushort posInLineCells = (ushort) (trafficPersonVec[p].posInLaneM);

     // laneMap[mapToWriteShift + maxWidth * (nextEdge + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells] = vInMpS;
     const uint posToSample = mapToWriteShift + kMaxMapWidthM * (nextEdge + (((int) (posInLineCells / kMaxMapWidthM)) * trafficPersonVec[p].edgeNumLanes) + trafficPersonVec[p].numOfLaneInEdge) + posInLineCells % kMaxMapWidthM;  // note the last % should not happen
     laneMap[posToSample] = vInMpS;
   }
}//

__global__ void kernel_intersectionSTOPSimulation(
     uint numIntersections, 
     float currentTime, 
     LC::B18IntersectionData *intersections, 
     uchar *trafficLights,
     LC::B18EdgeData* edgesData,//for the length
     uchar* laneMap,//to check if there are cars
     uint mapToReadShift) {
     int i = blockIdx.x * blockDim.x + threadIdx.x;
     if (i<numIntersections) {//CUDA check (inside margins)

     const float deltaEvent = 0.0f; 

     //if(i==0)printf("i %d\n",i);
     if (currentTime > intersections[i].nextEvent && intersections[i].totalInOutEdges > 0) {
       uint edgeOT = intersections[i].edge[intersections[i].state];
       uchar numLinesO = edgeOT >> 24;
       uint edgeONum = edgeOT & kMaskLaneMap; // 0xFFFFF

       // red old traffic lights
       for (int nL = 0; nL < numLinesO; nL++) {
         trafficLights[edgeONum + nL] = 0x00; //red old traffic light
       }

       for (int iN = 0; iN <= intersections[i].totalInOutEdges + 1; iN++) { //to give a round
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
		 ///
	 }
	 
}//

__global__ void kernel_intersectionOneSimulation(
      uint numIntersections,
      float currentTime,
      LC::B18IntersectionData *intersections,
      uchar *trafficLights) {
	 int i = blockIdx.x * blockDim.x + threadIdx.x;
	 if(i<numIntersections){//CUDA check (inside margins)
     const float deltaEvent = 20.0f; /// !!!!
     if (currentTime > intersections[i].nextEvent && intersections[i].totalInOutEdges > 0) {

       uint edgeOT = intersections[i].edge[intersections[i].state];
       uchar numLinesO = edgeOT >> 24;
       uint edgeONum = edgeOT & kMaskLaneMap; // 0xFFFFF;

       // red old traffic lights
       if ((edgeOT&kMaskInEdge) == kMaskInEdge) {  // Just do it if we were in in
         for (int nL = 0; nL < numLinesO; nL++) {
           trafficLights[edgeONum + nL] = 0x00; //red old traffic light
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
           }

           //trafficLights[edgeINum]=0xFF;
           break;
         }
       }//green new traffic light

       intersections[i].nextEvent = currentTime + deltaEvent;
     }
		 //////////////////////////////////////////////////////
	 }
	 
 }//

// Kernel that executes on the CUDA device
__global__ void kernel_sampleTraffic(
  int numPeople,
  LC::B18TrafficPerson *trafficPersonVec,
  uint *indexPathVec,
  float *accSpeedPerLinePerTimeInterval,
  float *numVehPerLinePerTimeInterval, //this could have been int
  uint offset
  ) {
  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p < numPeople) {//CUDA check (inside margins)
    if (trafficPersonVec[p].active == 1) { // just active
      int edgeNum = indexPathVec[trafficPersonVec[p].indexPathCurr];
      accSpeedPerLinePerTimeInterval[edgeNum + offset] += trafficPersonVec[p].v / 3.0f;
      numVehPerLinePerTimeInterval[edgeNum + offset]++;
    }
  }
}
__global__ void kernel_resetPeople(
  int numPeople,
  LC::B18TrafficPerson *trafficPersonVec) {
  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p < numPeople) {//CUDA check (inside margins)
    trafficPersonVec[p].active = 0;
  }
}

void b18GetSampleTrafficCUDA(std::vector<float>& accSpeedPerLinePerTimeInterval, std::vector<float>& numVehPerLinePerTimeInterval) {
  // copy back people
  size_t size = accSpeedPerLinePerTimeInterval.size() * sizeof(float);
  cudaMemcpy(accSpeedPerLinePerTimeInterval.data(), accSpeedPerLinePerTimeInterval_d, size, cudaMemcpyDeviceToHost);

  size_t sizeI = numVehPerLinePerTimeInterval.size() * sizeof(uchar);
  cudaMemcpy(numVehPerLinePerTimeInterval.data(), numVehPerLinePerTimeInterval_d, sizeI, cudaMemcpyDeviceToHost);
}

void b18ResetPeopleLanesCUDA(uint numPeople) {
  kernel_resetPeople << < ceil(numPeople / 1024.0f), 1024 >> > (numPeople, trafficPersonVec_d);
  cudaMemset(&laneMap_d[0], -1, halfLaneMap*sizeof(unsigned char));
  cudaMemset(&laneMap_d[halfLaneMap], -1, halfLaneMap*sizeof(unsigned char));
}

void b18SimulateTrafficCUDA(float currentTime, uint numPeople, uint numIntersections) {

	////////////////////////////////////////////////////////////
	// 1. CHANGE MAP: set map to use and clean the other
	if(readFirstMapC==true){
		mapToReadShift=0;
		mapToWriteShift=halfLaneMap;
    gpuErrchk(cudaMemset(&laneMap_d[halfLaneMap], -1, halfLaneMap*sizeof(unsigned char)));//clean second half
	}else{
		mapToReadShift=halfLaneMap;
		mapToWriteShift=0;
    gpuErrchk(cudaMemset(&laneMap_d[0], -1, halfLaneMap*sizeof(unsigned char)));//clean first half
	}
	readFirstMapC=!readFirstMapC;//next iteration invert use

  // Simulate intersections.
	kernel_intersectionOneSimulation <<< numIntersections,1 >>> (numIntersections, currentTime,intersections_d,trafficLights_d);
  gpuErrchk(cudaPeekAtLastError());
  
  // Simulate people.
  kernel_trafficSimulation << < ceil(numPeople / 256.0f), 256 >> > (numPeople, currentTime, mapToReadShift, mapToWriteShift, trafficPersonVec_d, indexPathVec_d, edgesData_d, laneMap_d, intersections_d, trafficLights_d);
  gpuErrchk(cudaPeekAtLastError());

  // Sample if necessary.
  if ((((float) ((int) currentTime)) == (currentTime)) &&
    ((int) currentTime % ((int) 30)) == 0) { //3min //(sample double each 3min)
    int samplingNumber = (currentTime - startTime) / (30 * numStepsTogether);
    uint offset = numIntersections * samplingNumber;
    //printf("Sample %d\n", samplingNumber);
    kernel_sampleTraffic << < ceil(numPeople / 1024.0f), 1024 >> > (numPeople, trafficPersonVec_d, indexPathVec_d, accSpeedPerLinePerTimeInterval_d, numVehPerLinePerTimeInterval_d, offset);
    gpuErrchk(cudaPeekAtLastError());
  }
}//
