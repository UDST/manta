//CUDA CODE
#include <stdio.h>
#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"

#include "b18TrafficPerson.h"
#include "b18EdgeData.h"
#include <vector>
#include <iostream>

#ifndef ushort
#define ushort uint16_t
#endif
#ifndef uint
#define uint uint32_t
#endif
#ifndef uchar
#define uchar uint8_t
#endif

///////////////////////////////
// CONSTANTS

__constant__ float intersectionClearance = 7.8f;
// `s_0` refers to the minimum spacing distance used in the Intelligent Driver Model (IDM)
__constant__ float s_0 = 7.0f;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}
inline void printMemoryUsage() {
  // show memory usage of GPU
  size_t free_byte;
  size_t total_byte;
  cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);
  if (cudaSuccess != cuda_status) {
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
    exit(1);
  }
  double free_db = (double) free_byte;
  double total_db = (double) total_byte;
  double used_db = total_db - free_db;
  printf("GPU memory usage: used = %.0f, free = %.0f MB, total = %.0f MB\n", used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
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

LC::Connection *deviceConnections;
size_t amountOfConnections;
LC::Intersection *deviceIntersections;
size_t amountOfIntersections;

void b18InitCUDA(
    bool fistInitialization,
    std::vector<LC::B18TrafficPerson>& trafficPersonVec,
    std::vector<uint> &indexPathVec,
    std::vector<LC::B18EdgeData>& edgesData,
    std::vector<uchar>& laneMap,
    std::vector<uchar>& trafficLights,
    std::vector<LC::B18IntersectionData>& b18Intersections,
    float startTimeH, float endTimeH,
    std::vector<float>& accSpeedPerLinePerTimeInterval,
    std::vector<float>& numVehPerLinePerTimeInterval,
    const std::vector<LC::Connection> & hostConnections,
    const std::vector<LC::Intersection> & hostIntersections) {

  { // Connections
    amountOfConnections = hostConnections.size();
    size_t size = hostConnections.size() * sizeof(LC::Connection);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &deviceConnections, size));   // Allocate array on device
    gpuErrchk(cudaMemcpy(deviceConnections, hostConnections.data(), size, cudaMemcpyHostToDevice));
  }

  { // Intersections
    amountOfIntersections = hostIntersections.size();
    size_t size = hostIntersections.size() * sizeof(LC::Intersection);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &deviceIntersections, size));   // Allocate array on device
    gpuErrchk(cudaMemcpy(deviceIntersections, hostIntersections.data(), size, cudaMemcpyHostToDevice));
  }

  printMemoryUsage();
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
  {// b18Intersections
    size_t sizeI = b18Intersections.size() * sizeof(LC::B18IntersectionData);
    if (fistInitialization) gpuErrchk(cudaMalloc((void **) &intersections_d, sizeI));   // Allocate array on device
    gpuErrchk(cudaMemcpy(intersections_d, b18Intersections.data(), sizeI, cudaMemcpyHostToDevice));
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
  printMemoryUsage();
}//

void b18FinishCUDA(void){
  cudaFree(deviceConnections);
  cudaFree(deviceIntersections);
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
   bool found = false;

   // CHECK FORWARD
   //printf("initShift %u numOfCells %u\n",initShift,numOfCells);
   for (ushort b = initShift - 1; (b < numOfCells) && (found == false); b++) { //NOTE -1 to make sure there is none in at the same level
     // laneChar = laneMap[mapToReadShift + maxWidth * (laneToCheck) + b];
     const uint posToSample = mapToReadShift + kMaxMapWidthM * (laneToCheck + (((int) (b / kMaxMapWidthM)) * numLinesEdge)) + b % kMaxMapWidthM;
     const uchar laneChar = laneMap[posToSample];

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
     const uint posToSample = mapToReadShift + kMaxMapWidthM * (laneToCheck + (((int) (b / kMaxMapWidthM)) * numLinesEdge)) + b % kMaxMapWidthM;
     const uchar laneChar = laneMap[posToSample];
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
   LC::B18IntersectionData* b18Intersections,
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

   for (int eN = b18Intersections[edgeNextInters].totalInOutEdges - 1; eN >= 0; eN--) {  // clockwise
     // retrieve
     uint procEdge = b18Intersections[edgeNextInters].edge[eN];

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
   const int numPeople,
   float currentTime,
   uint mapToReadShift,
   uint mapToWriteShift,
   LC::B18TrafficPerson *trafficPersonVec,
   uint *indexPathVec,
   LC::B18EdgeData* edgesData,
   uchar *laneMap,

   LC::B18IntersectionData *b18Intersections,
   // TODO: Remove usage of this data
   uchar *trafficLights,

   // TODO: Use only this data for computations
   LC::Connection *connections,
   size_t amountOfConnections,
   LC::Intersection *intersections,
   size_t amountOfIntersections)
 {
   const int p = blockIdx.x * blockDim.x + threadIdx.x;
   // Only proceed if the computed index `p` is valid
   if (p < numPeople) {
     /**
      * First ensure this person's car's info is initialized and whether is it active or not
      */
     if (trafficPersonVec[p].active == 2) {
       // Return if this person has reached its destiny
       return;
     }

     if (trafficPersonVec[p].active == 0){
       if (trafficPersonVec[p].time_departure > currentTime) {
       // Return if it's not yet the time for this person
         return;
       }

       const uint firstEdge = indexPathVec[trafficPersonVec[p].indexPathInit];
       if (firstEdge == -1) {
         // Return if this person's path has length zero
         trafficPersonVec[p].active = 2;
         return;
       }

       // Else initialize this person's data
       trafficPersonVec[p].indexPathCurr = trafficPersonVec[p].indexPathInit;
       trafficPersonVec[p].edgeNumLanes = edgesData[firstEdge].numLines;
       trafficPersonVec[p].edgeNextInters = edgesData[firstEdge].nextInters;
       trafficPersonVec[p].length = edgesData[firstEdge].length;
       trafficPersonVec[p].maxSpeedMperSec = edgesData[firstEdge].maxSpeedMperSec;

       // Find the starting position of the current person
       // At least `requiredAmountOfEmptyCells` are needed before the position where the car will be placed
       const ushort requiredAmountOfEmptyCells = s_0;
       const ushort startingRoadAmountOfCells = ceil(trafficPersonVec[p].length);
       // We will start to search from the middle of the starting road
       const ushort initShift = static_cast<ushort>(0.5f * startingRoadAmountOfCells);
       bool placed = false;
       ushort amountOfEmptySells = 0;
       for (ushort position = initShift; (position < startingRoadAmountOfCells) && (placed == false); position++) {
         const ushort numberOfRightLane = trafficPersonVec[p].edgeNumLanes - 1;
         const uchar laneChar = laneMap[mapToReadShift + kMaxMapWidthM * (firstEdge + numberOfRightLane) + position];
         if (laneChar != 0xFF) {
           // If the cell is not empty reset the empty-cells counter
           amountOfEmptySells = 0;
           continue;
         }

         // Keep advancing until enough empty cells have been found
         amountOfEmptySells++;
         if (amountOfEmptySells < requiredAmountOfEmptyCells) { continue; }

         // If we get to this point we can place the car
         trafficPersonVec[p].numOfLaneInEdge = numberOfRightLane;
         trafficPersonVec[p].posInLaneM = position; //m
         const uchar vInMpS = (uchar) (trafficPersonVec[p].v * 3); //speed in m/s *3 (to keep more precision
         laneMap[mapToWriteShift + kMaxMapWidthM * (firstEdge + numberOfRightLane) + position] = vInMpS;
         placed = true;
         break;
       }

       if (!placed) {
         // Return if the current road is too busy
         return;
       }

       trafficPersonVec[p].v = 0;
       trafficPersonVec[p].LC_stateofLaneChanging = 0;
       trafficPersonVec[p].active = 1;
       trafficPersonVec[p].isInIntersection = 0;
       trafficPersonVec[p].num_steps = 1;
       trafficPersonVec[p].co = 0.0f;
       trafficPersonVec[p].gas = 0.0f;

       const uint nextEdge = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];
       if (nextEdge != -1) {
         // TODO: Storing a reference to the next edge data would be safer (instead of plainly copying it)
         trafficPersonVec[p].nextEdgemaxSpeedMperSec = edgesData[nextEdge].maxSpeedMperSec;
         trafficPersonVec[p].nextEdgeNumLanes = edgesData[nextEdge].numLines;
         trafficPersonVec[p].nextEdgeNextInters = edgesData[nextEdge].nextInters;
         trafficPersonVec[p].nextEdgeLength = edgesData[nextEdge].length;
         trafficPersonVec[p].LC_initOKLanes = 0xFF;
         trafficPersonVec[p].LC_endOKLanes = 0xFF;
       }
       return;
     }

     // At this point we can assume trafficPersonVec[p].active == 1
     if (float(currentTime) == int(currentTime)) { // assuming deltatime = 0.5f --> each second
       trafficPersonVec[p].num_steps++;
     }

     /**
      * Gather enough information to know how the current car should be updated, using the Intelligent Driver Model (IDM)
      */
     const uint currentEdge = indexPathVec[trafficPersonVec[p].indexPathCurr];
     const uint nextEdge = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];
     float numMToMove;
     bool getToNextEdge = false;
     bool nextVehicleIsATrafficLight = false;
     float thirdTerm = 0;
     int remainingCellsToCheck = max(30.0f, trafficPersonVec[p].v * deltaTime * 2); //30 or double of the speed*time

     bool obstacleFound = false;
     bool noFirstInLaneBeforeIntersection = false; //use for stop control (just let 1st to pass)
     bool noFirstInLaneAfterSign = false; //use for stop control (just let 1st to pass)
     float s;
     float delta_v;
     const ushort byteInLine = (ushort) floor(trafficPersonVec[p].posInLaneM);
     const ushort numOfCells = ceil((trafficPersonVec[p].length - intersectionClearance));

     // Check if there is another car in the same lane
     for (ushort b = byteInLine + 2; (b < numOfCells) && (obstacleFound == false) && (remainingCellsToCheck > 0); b++, remainingCellsToCheck--) {
       const uint posToSample =
          mapToReadShift
          + kMaxMapWidthM * (
            indexPathVec[trafficPersonVec[p].indexPathCurr]
            + static_cast<int>(byteInLine / kMaxMapWidthM) * trafficPersonVec[p].edgeNumLanes
            + trafficPersonVec[p].numOfLaneInEdge)
          + b % kMaxMapWidthM;
       const uchar laneChar = laneMap[posToSample];

       if (laneChar != 0xFF) {
         s = ((float) (b - byteInLine)); //m
         delta_v = trafficPersonVec[p].v - (laneChar / 3.0f); //laneChar is in 3*ms (to save space in array)
         obstacleFound = true;
         noFirstInLaneBeforeIntersection = true;
         break;
       }
     }

     // At this point we found an obstacle or we reached the end of the current edge
     // If we are at the end of the current edge, check if this car's lane's connections are enabled
     if (byteInLine < numOfCells && !obstacleFound && remainingCellsToCheck > 0) {
       const int dstVertexNumber = edgesData[currentEdge].originalTargetVertexIndex;
       const auto currentLaneNumber = currentEdge + trafficPersonVec[p].numOfLaneInEdge;
       const auto nextEdgeNumber = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];
       bool atLeastOneEnabledConnection = false;
       for (int connectionIdx = intersections[dstVertexNumber].connectionGraphStart; connectionIdx < intersections[dstVertexNumber].connectionGraphEnd; ++connectionIdx) {
         // Check if a least one connection is enabled between the current edge and the following one
         const LC::Connection & connection = connections[connectionIdx];
         if (connection.inLaneNumber == currentLaneNumber
             && connection.outEdgeNumber == nextEdgeNumber
             && connection.enabled) {
           // TODO: Here I could store the available connection so that I don't need to make this cycle again later on
           atLeastOneEnabledConnection = true;
           break;
         }
       }

       // If no connection to the needed edge is enabled, then that intersection will be treated as an obstacle
       if (!atLeastOneEnabledConnection && nextEdgeNumber != -1) {
         s = ((float) (numOfCells - byteInLine));  // In meters
         delta_v = trafficPersonVec[p].v - 0;
         nextVehicleIsATrafficLight = true;
         obstacleFound = true;
       }
     }

     // Check if there is another car in the same lane after the intersection
     // TODO: At this point we need the information about which lane would be taken in case the intersection needs to be crossed
     for (ushort b = byteInLine + 2; (b < numOfCells) && (obstacleFound == false) && (remainingCellsToCheck > 0); b++, remainingCellsToCheck--) {
       const uint posToSample =
         mapToReadShift
         + kMaxMapWidthM * (
           indexPathVec[trafficPersonVec[p].indexPathCurr]
           + static_cast<int>(byteInLine / kMaxMapWidthM) * trafficPersonVec[p].edgeNumLanes
           + trafficPersonVec[p].numOfLaneInEdge)
         + b % kMaxMapWidthM;
       const uchar laneChar = laneMap[posToSample];

       if (laneChar != 0xFF) {
         s = ((float) (b - byteInLine)); //m
         delta_v = trafficPersonVec[p].v - (laneChar / 3.0f); //laneChar is in 3*ms (to save space in array)
         obstacleFound = true;
         noFirstInLaneAfterSign = true;
         break;
       }
     }

     // TODO: Confirm if all this if statement can be removed
     if (trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] == 0x0F && remainingCellsToCheck > 0) { //stop
       //check
       if (!noFirstInLaneBeforeIntersection
           && byteInLine < numOfCells //first before traffic
           && trafficPersonVec[p].v == 0 //stopped
           && !noFirstInLaneAfterSign) { // noone after the traffic light (otherwise wait before stop) !! Todo also check the beginning of next edge
         trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge] = 0x00; //reset stop
         trafficPersonVec[p].posInLaneM = ceilf(numOfCells) + 1; //move magicly after stop
       } else { //stop before STOP
         if (noFirstInLaneBeforeIntersection == false) { //just update this if it was the first one before sign
           s = ((float) (numOfCells - byteInLine)); //m
           delta_v = trafficPersonVec[p].v - 0; //it should be treated as an obstacle
           nextVehicleIsATrafficLight = true;
           obstacleFound = true;
         }
       }
     }

     // NEXT LINE
     if (obstacleFound == false && remainingCellsToCheck > 0) { //check if in next line
       if ((nextEdge != -1) && (trafficPersonVec[p].edgeNextInters != trafficPersonVec[p].end_intersection)) { // we haven't arrived to destination (check next line)
         ushort nextEdgeLaneToBe = trafficPersonVec[p].numOfLaneInEdge; //same lane

         if (nextEdgeLaneToBe >= trafficPersonVec[p].nextEdgeNumLanes) {
           nextEdgeLaneToBe = trafficPersonVec[p].nextEdgeNumLanes - 1; //change line if there are less roads
         }

         ushort numOfCells = ceil(trafficPersonVec[p].nextEdgeLength);

         for (ushort b = 0; (b < numOfCells) && (obstacleFound == false) && (remainingCellsToCheck > 0); b++, remainingCellsToCheck--) {
           const uint posToSample = mapToReadShift + kMaxMapWidthM * (nextEdge + nextEdgeLaneToBe) + b; // b18 not changed since we check first width
           const uchar laneChar = laneMap[posToSample];

           if (laneChar != 0xFF) {
             s = ((float) (b)); //m
             delta_v = trafficPersonVec[p].v - (laneChar / 3.0f);  // laneChar is in 3*ms (to save space in array)
             obstacleFound = true;
             break;
           }
         }
       }
     }

     /**
      * Update car's information
      */
     float s_star;
     if (obstacleFound) {
       s_star = s_0 + max(
         0.0f,
         trafficPersonVec[p].v * trafficPersonVec[p].T + (trafficPersonVec[p].v * delta_v) / (2 * sqrtf(trafficPersonVec[p].a * trafficPersonVec[p].b))
       );
       thirdTerm = powf(s_star / s, 2);
     }
     float dv_dt = trafficPersonVec[p].a * (1.0f - std::pow((trafficPersonVec[p].v / trafficPersonVec[p].maxSpeedMperSec), 4) - thirdTerm);

     numMToMove = max(0.0f, trafficPersonVec[p].v * deltaTime + 0.5f * (dv_dt) * deltaTime * deltaTime);
     trafficPersonVec[p].v += dv_dt * deltaTime;
     if (trafficPersonVec[p].v < 0) {
       trafficPersonVec[p].v = 0;
       dv_dt = 0.0f;
     }

     if (calculatePollution && ((float(currentTime) == int(currentTime)))) { // enabled and each second (assuming deltaTime 0.5f)
       // Note: compute CO and Gas values each second

       // CO Calculation
       const float speedMph = trafficPersonVec[p].v * 2.2369362920544; //mps to mph
       const float COStepPerSecond = -0.064 + 0.0056 * speedMph + 0.00026 * (speedMph - 50.0f) * (speedMph - 50.0f);
       if (COStepPerSecond > 0) { trafficPersonVec[p].co += COStepPerSecond; }

       // Gas Consumption
       const float a = dv_dt;
       const float v = trafficPersonVec[p].v; // in mps
       const float Pea = a > 0.0f ? (0.472f*1.680f*a*a*v) : 0.0f;
       const float gasStepPerSecond = 0.666f + 0.072f*(0.269f*v + 0.000672f*(v*v*v) + 0.0171f*(v*v) + 1.680f*a*v + Pea);
       trafficPersonVec[p].gas += gasStepPerSecond;
     }

     if (trafficPersonVec[p].v == 0) {
       const ushort posInLineCells = static_cast<ushort>(trafficPersonVec[p].posInLaneM);
       const uint posToSample =
         mapToWriteShift
         + kMaxMapWidthM * (
             currentEdge
             + static_cast<int>(posInLineCells / kMaxMapWidthM) * trafficPersonVec[p].edgeNumLanes
             + trafficPersonVec[p].numOfLaneInEdge)
         + posInLineCells % kMaxMapWidthM;
       laneMap[posToSample] = 0;

       return;
     }

     trafficPersonVec[p].color = p << 8;
     trafficPersonVec[p].posInLaneM = trafficPersonVec[p].posInLaneM + numMToMove;

     if (trafficPersonVec[p].posInLaneM > trafficPersonVec[p].length) { //reach intersection
       numMToMove = trafficPersonVec[p].posInLaneM - trafficPersonVec[p].length;
       getToNextEdge = true;
     } else { //does not research next intersection
       // If the intersection has not been reached try to changed lane if:
       //   - The car is going at least 10 km per hour
       //   - 5 seconds have happened since the last lane change
       if (trafficPersonVec[p].v > 3.0f && trafficPersonVec[p].num_steps % 5 == 0) {
         // next thing is not a traffic light
         // skip if there is one lane (avoid to do this)
         // skip if it is the last edge
         if (nextVehicleIsATrafficLight == false && trafficPersonVec[p].edgeNumLanes > 1 && nextEdge != -1) {
           ////////////////////////////////////////////////////
           // LC 1 update lane changing status
           if (trafficPersonVec[p].LC_stateofLaneChanging == 0) {
             // 2.2-exp((x-1)^2)
             const float x = trafficPersonVec[p].posInLaneM / trafficPersonVec[p].length;
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
               bool leftLane = trafficPersonVec[p].numOfLaneInEdge > 0; //at least one lane on the left
               bool rightLane = trafficPersonVec[p].numOfLaneInEdge < trafficPersonVec[p].edgeNumLanes - 1; //at least one lane

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
               // TODO: Replace the following line by the a value indicating whether the corresponding intersection is enabled
               uchar trafficLightState = trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge];
               calculateGapsLC(mapToReadShift, laneMap, trafficLightState,
                 currentEdge + laneToCheck, trafficPersonVec[p].edgeNumLanes, trafficPersonVec[p].posInLaneM,
                 trafficPersonVec[p].length, v_a, v_b, gap_a, gap_b);

               if (gap_a == 1000.0f && gap_b == 1000.0f) { //lag and lead car very far
                 trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE
               } else { // NOT ALONE
                 float b1A = 0.05f, b2A = 0.15f;
                 float b1B = 0.15f, b2B = 0.40f;
                 // s_0-> critical lead gap
                 float g_na_D, g_bn_D;
                 bool acceptLC = true;

                 if (gap_a != 1000.0f) {
                   g_na_D = max(s_0, s_0 + b1A * trafficPersonVec[p].v + b2A *
                     (trafficPersonVec[p].v - v_a * 3.0f));

                   if (gap_a < g_na_D) { //gap smaller than critical gap
                     acceptLC = false;
                   }
                 }

                 if (acceptLC == true && gap_b != 1000.0f) {
                   g_bn_D = max(s_0, s_0 + b1B * v_b * 3.0f + b2B * (v_b * 3.0f - trafficPersonVec[p].v));

                   if (gap_b < g_bn_D) { //gap smaller than critical gap
                     acceptLC = false;
                   }
                 }

                 if (acceptLC == true) {
                   trafficPersonVec[p].numOfLaneInEdge = laneToCheck; // CHANGE LINE
                 }
               }
             }
           }// Discretionary

           ////////////////////////////////////////////////////
           // LC 3 *MANDATORY* STATE
           if (trafficPersonVec[p].LC_stateofLaneChanging == 1) {
             // LC 3.1 Calculate the correct lanes
             if (trafficPersonVec[p].LC_endOKLanes == 0xFF) {
               calculateLaneCarShouldBe(currentEdge, nextEdge, b18Intersections,
                 trafficPersonVec[p].edgeNextInters, trafficPersonVec[p].edgeNumLanes,
                 trafficPersonVec[p].LC_initOKLanes, trafficPersonVec[p].LC_endOKLanes);

               //printf("p%u num lanes %u min %u max %u\n",p,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
               if (trafficPersonVec[p].LC_initOKLanes == 0 &&
                 trafficPersonVec[p].LC_endOKLanes == 0) {
                 //exit(0);
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
                 //exit(0);
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
               // TODO: Replace the following line by the a value indicating whether the corresponding intersection is enabled
               uchar trafficLightState = trafficLights[currentEdge + trafficPersonVec[p].numOfLaneInEdge];
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
                   g_na_M = max(s_0, s_0 + (b1A * trafficPersonVec[p].v + b2A *
                     (trafficPersonVec[p].v - v_a * 3.0f)));

                   if (gap_a < g_na_M) { //gap smaller than critical gap
                     acceptLC = false;
                   }
                 }

                 if (acceptLC == true && gap_b != 1000.0f) {
                   g_bn_M = max(s_0, s_0 + (b1B * v_b * 3.0f + b2B * (v_b * 3.0f -
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

       // Update person' speed
       const uchar vInMpS = (uchar) (trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
       const ushort posInLineCells = (ushort) (trafficPersonVec[p].posInLaneM);
       const uint posToSample =
         mapToWriteShift
         + kMaxMapWidthM * (
             currentEdge
             + static_cast<int>(posInLineCells / kMaxMapWidthM) * trafficPersonVec[p].edgeNumLanes
             + trafficPersonVec[p].numOfLaneInEdge)
         + posInLineCells % kMaxMapWidthM;
       laneMap[posToSample] = vInMpS;
       return;
     }

     if (nextEdge == -1) {
       trafficPersonVec[p].active = 2;
       return;
     }

     // Update current edge information
     trafficPersonVec[p].indexPathCurr++;
     trafficPersonVec[p].maxSpeedMperSec = trafficPersonVec[p].nextEdgemaxSpeedMperSec;
     trafficPersonVec[p].edgeNumLanes = trafficPersonVec[p].nextEdgeNumLanes;
     trafficPersonVec[p].edgeNextInters = trafficPersonVec[p].nextEdgeNextInters;
     trafficPersonVec[p].length = trafficPersonVec[p].nextEdgeLength;
     trafficPersonVec[p].posInLaneM = numMToMove;

     if (trafficPersonVec[p].numOfLaneInEdge >= trafficPersonVec[p].edgeNumLanes) {
       trafficPersonVec[p].numOfLaneInEdge = trafficPersonVec[p].edgeNumLanes - 1; //change line if there are less roads
     }

     // Update person's next edge
     const uint nextEdgeIdx = indexPathVec[trafficPersonVec[p].indexPathCurr + 1];

     if (nextEdgeIdx != -1) {
       trafficPersonVec[p].LC_initOKLanes = 0xFF;
       trafficPersonVec[p].LC_endOKLanes = 0xFF;

       trafficPersonVec[p].nextEdgemaxSpeedMperSec = edgesData[nextEdgeIdx].maxSpeedMperSec;
       trafficPersonVec[p].nextEdgeNumLanes = edgesData[nextEdgeIdx].numLines;
       trafficPersonVec[p].nextEdgeNextInters = edgesData[nextEdgeIdx].nextInters;
       trafficPersonVec[p].nextEdgeLength = edgesData[nextEdgeIdx].length;
     }

     trafficPersonVec[p].LC_stateofLaneChanging = 0;
     const uchar vInMpS = static_cast<uchar>(trafficPersonVec[p].v * 3); //speed in m/s to fit in uchar
     const ushort posInLineCells = static_cast<ushort>(trafficPersonVec[p].posInLaneM);

     const uint posToSample =
       mapToWriteShift
       + kMaxMapWidthM * (
           nextEdge
           + static_cast<int>(posInLineCells / kMaxMapWidthM) * trafficPersonVec[p].edgeNumLanes
           + trafficPersonVec[p].numOfLaneInEdge)
       + posInLineCells % kMaxMapWidthM;  // note the last % should not happen
     laneMap[posToSample] = vInMpS;
   }
}

__global__ void kernel_intersectionOneSimulation(
    uint numIntersections,
    float currentTime,
    LC::B18IntersectionData *b18Intersections,
    uchar *trafficLights) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i<numIntersections){
    const float deltaEvent = 20.0f;  // 20 seconds between each change in the traffic lights
    if (currentTime > b18Intersections[i].nextEvent && b18Intersections[i].totalInOutEdges > 0) {
      uint edgeOT = b18Intersections[i].edge[b18Intersections[i].state];
      uchar numLinesO = edgeOT >> 24;
      uint edgeONum = edgeOT & kMaskLaneMap; // 0xFFFFF;

      // red old traffic lights
      if ((edgeOT&kMaskInEdge) == kMaskInEdge) {  // Just do it if we were in in
        for (int nL = 0; nL < numLinesO; nL++) {
          trafficLights[edgeONum + nL] = 0x00; //red old traffic light
        }
      }

      for (int iN = 0; iN <= b18Intersections[i].totalInOutEdges + 1; iN++) { //to give a round
        b18Intersections[i].state = (b18Intersections[i].state + 1) % b18Intersections[i].totalInOutEdges;//next light
        if ((b18Intersections[i].edge[b18Intersections[i].state] & kMaskInEdge) == kMaskInEdge) {  // 0x800000
          // green new traffic lights
          uint edgeIT = b18Intersections[i].edge[b18Intersections[i].state];
          uint edgeINum = edgeIT & kMaskLaneMap; //  0xFFFFF; //get edgeI
          uchar numLinesI = edgeIT >> 24;

          for (int nL = 0; nL < numLinesI; nL++) {
            trafficLights[edgeINum + nL] = 0xFF;
          }

          //trafficLights[edgeINum]=0xFF;
          break;
        }
      }//green new traffic light
      b18Intersections[i].nextEvent = currentTime + deltaEvent;
    }
  }
}

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

  // Update intersections.
  kernel_intersectionOneSimulation<<<ceil(numIntersections / 512.0f), 512>>>(
    numIntersections,
    currentTime,
    intersections_d,
    trafficLights_d);
  gpuErrchk(cudaPeekAtLastError());

  // Simulate people.
  kernel_trafficSimulation<<<ceil(numPeople / 384.0f), 384>>>(
    numPeople,
    currentTime,
    mapToReadShift,
    mapToWriteShift,
    trafficPersonVec_d,
    indexPathVec_d,
    edgesData_d,
    laneMap_d,
    intersections_d,
    trafficLights_d,
    deviceConnections,
    amountOfConnections,
    deviceIntersections,
    amountOfIntersections);
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
