//CUDA CODE
#include <stdio.h>
#include <vector>
#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"

#include "cudaTrafficPerson.h"
#include "cudaEdgeData.h"

#define ushort unsigned short
#define uint unsigned int
#define uchar unsigned char

///////////////////////////////
// CONSTANTS
//__constant__ float deltaTimeG=1.0f*(1.0f/3600.0f);
__constant__ float s_0=7.0f;
//const float deltaTimeC=1.0f*(1.0f/3600.0f);
////////////////////////////////
// VARIABLES
LC::CUDATrafficPerson *trafficPersonVec_d;
//ushort *nextEdgeM_d;
LC::edgeData *edgesData_d;
float cellSize;
float deltaTime;
uchar *laneMap_d;
bool readFirstMapC=true;
uint mapToReadShift;
uint mapToWriteShift;
uint halfLaneMap;
ushort maxWidth;

LC::intersectionData *intersections_d;
uchar *trafficLights_d;

void initCUDA(float _cellSize,float _deltaTime,ushort _maxWidth,std::vector<LC::CUDATrafficPerson>& trafficPersonVec,std::vector<LC::edgeData>& edgesData,std::vector<uchar>& laneMap,std::vector<LC::intersectionData>& intersections){
	cellSize=_cellSize;
	deltaTime=_deltaTime;
	maxWidth=_maxWidth;
	// people
	size_t size = trafficPersonVec.size() * sizeof(LC::CUDATrafficPerson);
	cudaError err ;
	err=cudaMalloc((void **) &trafficPersonVec_d, size);   // Allocate array on device
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	err =cudaMemcpy(trafficPersonVec_d,trafficPersonVec.data(),size,cudaMemcpyHostToDevice);
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	
	/*//nextEdge
	size_t sizeE = nextEdgeM.size() * sizeof(ushort);
	cudaMalloc((void **) &nextEdgeM_d, sizeE);   // Allocate array on device
	cudaMemcpy(nextEdgeM_d,nextEdgeM.data(),sizeE,cudaMemcpyHostToDevice);*/
	//edgeData
	size_t sizeD = edgesData.size() * sizeof(LC::edgeData);
	err=cudaMalloc((void **) &edgesData_d, sizeD);   // Allocate array on device
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	err=cudaMemcpy(edgesData_d,edgesData.data(),sizeD,cudaMemcpyHostToDevice);
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	//laneMap
	size_t sizeL = laneMap.size() * sizeof(uchar);
	err=cudaMalloc((void **) &laneMap_d, sizeL);   // Allocate array on device
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	err=cudaMemcpy(laneMap_d,laneMap.data(),sizeL,cudaMemcpyHostToDevice);
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	halfLaneMap=laneMap.size()/2;
	// intersections
	size_t sizeI = intersections.size() * sizeof(LC::intersectionData);
	err=cudaMalloc((void **) &intersections_d, sizeI);   // Allocate array on device
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	err=cudaMemcpy(intersections_d,intersections.data(),sizeI,cudaMemcpyHostToDevice);
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
	size_t sizeT = (halfLaneMap/_maxWidth) * sizeof(uchar);//total number of lanes
	err=cudaMalloc((void **) &trafficLights_d, sizeT);   // Allocate array on device
	if ( cudaSuccess != err )fprintf( stderr, "Cuda error: %s.\n",cudaGetErrorString( err) );
}//

void finishCUDA(void){
	//////////////////////////////
	// FINISH
	cudaFree(trafficPersonVec_d);
	//cudaFree(nextEdgeM_d);
	cudaFree(edgesData_d);
	cudaFree(laneMap_d);
	cudaFree(intersections_d);
	cudaFree(trafficLights_d);

}//

 void getDataCUDA(std::vector<LC::CUDATrafficPerson>& trafficPersonVec,std::vector<uchar>& trafficLights){
	 // copy back people
	 size_t size = trafficPersonVec.size() * sizeof(LC::CUDATrafficPerson);
	 cudaMemcpy(trafficPersonVec.data(),trafficPersonVec_d,size,cudaMemcpyDeviceToHost);//cudaMemcpyHostToDevice

	 size_t sizeI = trafficLights.size() * sizeof(uchar);
	 cudaMemcpy(trafficLights.data(),trafficLights_d,sizeI,cudaMemcpyDeviceToHost);
 }


 __device__ void calculateGapsLC(
		float cellSize,
		ushort maxWidth,
		uint mapToReadShift,
		uchar *laneMap,
		uchar trafficLightState,
		ushort laneToCheck,
		float posInMToCheck,
		float length,
		uchar &v_a,
		uchar &v_b,
		float &gap_a,
		float &gap_b){
		ushort numOfCells=ceil(length/cellSize);
		ushort initShift=ceil(posInMToCheck/cellSize);
		uchar laneChar;
		bool found=false;
		// CHECK FORWARD
		//printf("initShift %u numOfCells %u\n",initShift,numOfCells);
		for(ushort b=initShift-1;(b<numOfCells)&&(found==false);b++){//NOTE -1 to make sure there is none in at the same level
			laneChar=laneMap[mapToReadShift+maxWidth*(laneToCheck)+b];
			if(laneChar!=0xFF){
				gap_a=((float)b-initShift)*cellSize;//m
				v_a=laneChar;//laneChar is in 3*ms (to save space in array)
				found=true;
				break;
			}
		}
		if(found==false){
			if(trafficLightState==0x00){//red
				//gap_a=((float)numOfCells-initShift)*cellSize;
				//found=true;
				gap_a=gap_b=1000.0f;//force to change to the line without vehicle
				v_a=v_b=0xFF;
				return;
			}
		}
		if(found==false){
			gap_a=1000.0f;
		}
		// CHECK BACKWARDS
		found=false;
		//printf("2initShift %u numOfCells %u\n",initShift,numOfCells);
		for(int b=initShift+1;(b>=0)&&(found==false);b--){//NOTE +1 to make sure there is none in at the same level
			laneChar=laneMap[mapToReadShift+maxWidth*(laneToCheck)+b];
			if(laneChar!=0xFF){
				gap_b=((float)initShift-b)*cellSize;//m
				v_b=laneChar;//laneChar is in 3*ms (to save space in array)
				found=true;
				break;
			}
		}
		//printf("3initShift %u numOfCells %u\n",initShift,numOfCells);
		if(found==false){
			gap_b=1000.0f;
		}

	}//

 __device__ void calculateLaneCarShouldBe(
		ushort curEdgeLane,
		ushort nextEdge,
		LC::intersectionData *intersections,
		ushort edgeNextInters,
		ushort edgeNumLanes,
		ushort &initOKLanes,
		ushort &endOKLanes){
		initOKLanes=0;
		endOKLanes=edgeNumLanes;
		/*if(DEBUG_TRAFFIC==1)printf("curEdgeLane %05x nextEdge %05x\n",curEdgeLane,nextEdge);
		if(DEBUG_TRAFFIC==1){
			for(int eN=0;eN<intersections[edgeNextInters].totalInOutEdges;eN++){
				printf("* procEdge %05x\n",intersections[edgeNextInters].edge[eN]);
			}
		}*/
		bool currentEdgeFound=false;
		bool exitFound=false;
		ushort numExitToTake=0;
		ushort numExists=0;
		for(int eN=intersections[edgeNextInters].totalInOutEdges-1;eN>=0;eN--){//clockwise
			uint procEdge=intersections[edgeNextInters].edge[eN];
			if((procEdge&0xFFFF)==curEdgeLane){//current edge
				//if(DEBUG_TRAFFIC==1)printf("CE procEdge %05x\n",procEdge);
				currentEdgeFound=true;
				if(exitFound==false)
					numExitToTake=0;
				continue;
			}
			
			
			if((procEdge&0x010000)==0x0){//out edge
				//if(DEBUG_TRAFFIC==1)printf("   procEdge %05x\n",procEdge);
				numExists++;
				if(currentEdgeFound==true){
					numExitToTake++;
				}
				if(currentEdgeFound==false&&exitFound==false){
					numExitToTake++;
				}
			}
			if((procEdge&0xFFFF)==nextEdge){
				exitFound=true;
				currentEdgeFound=false;
				//if(DEBUG_TRAFFIC==1)printf("NE procEdge %05x\n",procEdge);
			}
		}
		//if(DEBUG_TRAFFIC==1)printf("Num extis %u Num exit to take %u%\n",numExists,numExitToTake);
		if(edgeNumLanes==0){
			return;//printf("ERRRROR\n");
		}
		switch(edgeNumLanes){
			/// ONE LANE
		case 1:
			initOKLanes=0;
			endOKLanes=1;
			break;
			/// TWO LANE
		case 2:
			switch(numExists){
			case 1:
			case 2://all okay
				initOKLanes=0;
				endOKLanes=2;
				break;
			case 3:
				if(numExitToTake>2){//left
					initOKLanes=0;
					endOKLanes=1;
					break;
				}
				initOKLanes=1;
				endOKLanes=2;
				break;
			default:
				if(numExitToTake>=numExists-1){
					initOKLanes=0;
					endOKLanes=1;
					break;
				}
				initOKLanes=1;
				endOKLanes=2;
				break;
			}
			break;
			/// THREE LANE
		case 3:
			switch(numExists){
			case 1:
			case 2://all okay
				initOKLanes=0;
				endOKLanes=3;
				break;
			case 3:
				if(numExitToTake>2){//left
					initOKLanes=0;
					endOKLanes=1;
					break;
				}
				initOKLanes=1;
				endOKLanes=3;
				break;
			default:
				if(numExitToTake>=numExists-1){
					initOKLanes=0;
					endOKLanes=1;
					break;
				}
				initOKLanes=1;
				endOKLanes=2;
				break;
			}
			break;
		case 4:
			switch(numExists){
			case 1:
			case 2://all okay
				initOKLanes=0;
				endOKLanes=4;
				break;
			case 3:
				if(numExitToTake==1){//right
					initOKLanes=3;
					endOKLanes=4;
				}
				if(numExitToTake>3){//left
					initOKLanes=0;
					endOKLanes=1;
					break;
				}
				initOKLanes=1;
				endOKLanes=4;
				break;
			default:
				if(numExitToTake==1){//right
					initOKLanes=edgeNumLanes-1;
					endOKLanes=edgeNumLanes;
				}
				if(numExitToTake>=numExists-2){
					initOKLanes=0;
					endOKLanes=2;
					break;
				}
				initOKLanes=1;//also lane 2
				endOKLanes=edgeNumLanes;
			}
			break;
		default:
			switch(numExists){
			case 1:
			case 2://all okay
				initOKLanes=0;
				endOKLanes=edgeNumLanes;
				break;
			case 3:
				if(numExitToTake==1){//right
					initOKLanes=edgeNumLanes-1;
					endOKLanes=edgeNumLanes;
				}
				if(numExitToTake>edgeNumLanes-2){//left
					initOKLanes=0;
					endOKLanes=2;
					break;
				}
				initOKLanes=1;
				endOKLanes=edgeNumLanes;
				break;
			default:
				if(numExitToTake<2){//right
					initOKLanes=edgeNumLanes-2;
					endOKLanes=edgeNumLanes;
				}
				if(numExitToTake>=numExists-2){
					initOKLanes=0;
					endOKLanes=2;
					break;
				}
				initOKLanes=1;//also lane 2
				endOKLanes=edgeNumLanes-1;
			}
			break;

		}
	}//

__device__ int cuda_qrand(){
	return 10;
}

 // Kernel that executes on the CUDA device
 __global__ void kernel_trafficSimulation(
	 int numPeople,
	 float currentTime,
	 float cellSize,
	 float deltaTime,
	 LC::CUDATrafficPerson *trafficPersonVec,
	 //ushort *nextEdgeM,
	 LC::edgeData* edgesData,
	 uchar *laneMap,
	 LC::intersectionData *intersections,
	 uchar *trafficLights,
	 uint mapToReadShift,
	 uint mapToWriteShift,
	 ushort maxWidth)
 {
	 const bool DEBUG_TRAFFIC=0;
	 int p = blockIdx.x * blockDim.x + threadIdx.x;
	 //printf("p %d Numpe %d\n",p,numPeople);
	 if(p<numPeople){//CUDA check (inside margins)
		 ///
		 ///////////////////////////////
		//2.0. check if finished
		if(trafficPersonVec[p].active==2){
				return;
			}
			///////////////////////////////
			//2.1. check if person should still wait or should start
			if(trafficPersonVec[p].active==0){

				//printf("  1. Person: %d active==0\n",p);
				if(trafficPersonVec[p].time_departure>currentTime){//wait
					//1.1 just continue waiting
					//printf("   1.1 Person: %d wait\n",p);
					return;
				}else{//start
					//1.2 find first edge
					trafficPersonVec[p].currPathEdge=0;
					ushort firstEdge=trafficPersonVec[p].personPath[0];
					if(firstEdge==0xFFFF){
						trafficPersonVec[p].active=2;
						//printf("0xFFFF\n");
						return;
					}

					
					//1.3 update person edgeData
					//if(DEBUG_TRAFFIC==1)printf("   1.3 Person: %d put in first edge %u\n",p,firstEdge);
					//printf("edgesData %d\n",edgesData);

					// COPY DATA FROM EDGE TO PERSON
					trafficPersonVec[p].edgeNumLanes=edgesData[firstEdge].numLines;
					trafficPersonVec[p].edgeNextInters=edgesData[firstEdge].nextInters;

					trafficPersonVec[p].length=edgesData[firstEdge].length;
					trafficPersonVec[p].maxSpeedMperSec=edgesData[firstEdge].maxSpeedMperSec;
					//printf("edgesData %.10f\n",edgesData[firstEdge].maxSpeedCellsPerDeltaTime);
					//1.4 try to place it in middle of edge
					ushort numOfCells=ceil(trafficPersonVec[p].length/cellSize);
					ushort initShift=(ushort)(0.5f*numOfCells);//number of cells it should be placed (half of road)

					uchar laneChar;
					bool placed=false;

					ushort numCellsEmptyToBePlaced=s_0/cellSize;
					ushort countEmptyCells=0;
					for(ushort b=initShift;(b<numOfCells)&&(placed==false);b++){
						//for(int lN=trafficPersonVec[p].edgeNumLanes-1;lN>=0;lN--){
						//ushort lN=0;//just first LANE !!!!!!!
						ushort lN=trafficPersonVec[p].edgeNumLanes-1;//just right LANE !!!!!!!
						laneChar=laneMap[mapToReadShift+maxWidth*(firstEdge+lN)+b];//get byte of edge (proper line)
						if(laneChar!=0xFF){
							countEmptyCells=0;
							continue;
						}
						countEmptyCells++;// ensure there is enough room to place the car
						if(countEmptyCells<numCellsEmptyToBePlaced){
							continue;
						}
						trafficPersonVec[p].numOfLaneInEdge=lN;
						trafficPersonVec[p].posInLaneM=b*cellSize;//m
						uchar vInMpS=(uchar)(trafficPersonVec[p].v*3);//speed in m/s *3 (to keep more precision
						laneMap[mapToWriteShift+maxWidth*(firstEdge+lN)+b]=vInMpS;
						placed=true;
						//printf("Placed\n");
						break;
						//}
					}
					if(placed==false){//not posible to start now
						return;
					}
					trafficPersonVec[p].v=0;//trafficPersonVec[p].maxSpeedCellsPerDeltaTime;//(20000.0f*deltaTime)/cellSize;//20km/h-->cell/delta time
					trafficPersonVec[p].LC_stateofLaneChanging=0;
					//1.5 active car
					
					trafficPersonVec[p].active=1;
					trafficPersonVec[p].num_steps=1;
					trafficPersonVec[p].gas=0;
					//trafficPersonVec[p].nextPathEdge++;//incremet so it continues in next edge
					// set up next edge info
					ushort nextEdge=trafficPersonVec[p].personPath[1];
					//trafficPersonVec[p].nextEdge=nextEdge;
					if(nextEdge!=0xFFFF){
						trafficPersonVec[p].nextEdgemaxSpeedMperSec=edgesData[nextEdge].maxSpeedMperSec;
						trafficPersonVec[p].nextEdgeNumLanes=edgesData[nextEdge].numLines;
						trafficPersonVec[p].nextEdgeNextInters=edgesData[nextEdge].nextInters;
						trafficPersonVec[p].nextEdgeLength=edgesData[nextEdge].length;
						//trafficPersonVec[p].nextPathEdge++;
						trafficPersonVec[p].LC_initOKLanes=0xFF;
						trafficPersonVec[p].LC_endOKLanes=0xFF;
					}
					return;
				}
			}
			
			///////////////////////////////
			//2. it is moving
			trafficPersonVec[p].num_steps++;
			//2.1 try to move
			float numMToMove;
			bool getToNextEdge=false;
			bool nextVehicleIsATrafficLight=false;
			ushort currentEdge=trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge];
			ushort nextEdge=trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge+1];
			//if(trafficPersonVec[p].posInLaneM<trafficPersonVec[p].length){
			// www.vwi.tu-dresden.de/~treiber/MicroApplet/IDM.html
			// IDM
			float thirdTerm=0;
			///////////////////////////////////////////////////
			// 2.1.1 Find front car
			//int numCellsCheck=fmax<float>(15.0f,trafficPersonVec[p].v*deltaTime*2);//15 or double of the speed*time
			// SAME LINE
			bool found=false;
			float s;
			float delta_v;
			uchar laneChar;
			ushort byteInLine=(ushort)floor(trafficPersonVec[p].posInLaneM/cellSize);
			ushort numOfCells=ceil(trafficPersonVec[p].length/cellSize);
			for(ushort b=byteInLine+2;(b<numOfCells)&&(found==false)/*&&(numCellsCheck>0)*/;b++/*,numCellsCheck--*/){
				laneChar=laneMap[mapToReadShift+maxWidth*(trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge]+trafficPersonVec[p].numOfLaneInEdge)+b];
				if(laneChar!=0xFF){
					s=((float)(b-byteInLine)*cellSize);//m
					delta_v=trafficPersonVec[p].v-(laneChar/3.0f);//laneChar is in 3*ms (to save space in array)
					found=true;
					//printf("\nFOUND Car Same Lane s %f delta_v %f!!!!\n",s,delta_v);
					break;
				}
			}
			// TRAFFIC LIGHT
			if(found==false){//check if traffic light is red
				if(trafficLights[currentEdge]==0x00){//red
					s=((float)(numOfCells-byteInLine)*cellSize);//m
					delta_v=trafficPersonVec[p].v-0;//it should be treated as an obstacle
					nextVehicleIsATrafficLight=true;
					//printf("\nFOUND TL\n",s,delta_v);
					found=true;
				}
			}
			// NEXT LINE
			if(found==false){//check if in next line
				if((nextEdge!=0xFFFF)&&(trafficPersonVec[p].edgeNextInters!=trafficPersonVec[p].end_intersection)){// we haven't arrived to destination (check next line)

					ushort nextEdgeLaneToBe=trafficPersonVec[p].numOfLaneInEdge;//same lane
					//printf("trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
					if(nextEdgeLaneToBe>=trafficPersonVec[p].nextEdgeNumLanes){
						nextEdgeLaneToBe=trafficPersonVec[p].nextEdgeNumLanes-1;//change line if there are less roads
					}
					//printf("2trafficPersonVec[p].numOfLaneInEdge %u\n",trafficPersonVec[p].numOfLaneInEdge);
					ushort numOfCells=ceil(trafficPersonVec[p].nextEdgeLength/cellSize);
					for(ushort b=0;(b<numOfCells)&&(found==false)/*&&(numCellsCheck>0)*/;b++/*,numCellsCheck--*/){
						laneChar=laneMap[mapToReadShift+maxWidth*(nextEdge+nextEdgeLaneToBe)+b];
						if(laneChar!=0xFF){
							s=((float)(b)*cellSize);//m
							delta_v=trafficPersonVec[p].v-(laneChar/3.0f);//laneChar is in 3*ms (to save space in array)
							found=true;
							//printf("\FOUND Car Next Lane s %f delta_v %f!!!!\n",s,delta_v);
							//printf("\FOUND Car Next Lane next Edge %u NumLanes %u Be in Lane %u!!!!\n",trafficPersonVec[p].nextEdge,trafficPersonVec[p].nextEdgeNumLanes,nextEdgeLaneToBe);
							break;
						}
					}
				}

			}
			float s_star;
			if(found==true){//car in front and slower than us
				// 2.1.2 calculate dv_dt
				s_star=s_0+fmax(0.0f,(trafficPersonVec[p].v*trafficPersonVec[p].T+(trafficPersonVec[p].v*delta_v)/(2*std::sqrt(trafficPersonVec[p].a*trafficPersonVec[p].b))));
				thirdTerm=std::pow(((s_star)/(s)),2);
				//printf(">FOUND s_star %f thirdTerm %f!!!!\n",s_star,thirdTerm);
			}

			float dv_dt= trafficPersonVec[p].a*(1.0f-std::pow((trafficPersonVec[p].v/trafficPersonVec[p].maxSpeedMperSec),4)-thirdTerm);

			// 2.1.3 update values
			numMToMove=fmax(0.0f,trafficPersonVec[p].v*deltaTime+0.5f*(dv_dt)*deltaTime*deltaTime);

			
			//printf("v %.10f v d %.10f\n",trafficPersonVec[p].v,trafficPersonVec[p].v+((dv_dt/(deltaTime)/deltaTime)));
			trafficPersonVec[p].v+=dv_dt*deltaTime;
			if(trafficPersonVec[p].v<0){
				//printf("p %d v %f v0 %f a %f dv_dt %f s %f s_star %f MOVE %f\n",p,trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec,trafficPersonVec[p].a,dv_dt,s,s_star,numMToMove);
				trafficPersonVec[p].v=0;
			}
			/////
			//CO2
			//if(trafficPersonVec[p].v>0)
			{
				float speedMph=trafficPersonVec[p].v*2.2369362920544;//mps to mph
				float gasStep=-0.064+0.0056*speedMph+0.00026*(speedMph-50.0f)*(speedMph-50.0f);
				if(gasStep>0){
					gasStep*=deltaTime;
					trafficPersonVec[p].gas+=gasStep;
				}
			}
			//trafficPersonVec[p].gas+=numMToMove/1000.0f;
			//////////////////////////////////////////////

			if(trafficPersonVec[p].v==0){//if not moving not do anything else
				ushort posInLineCells=(ushort)(trafficPersonVec[p].posInLaneM/cellSize);
				laneMap[mapToWriteShift+maxWidth*(currentEdge+trafficPersonVec[p].numOfLaneInEdge)+posInLineCells]=0;
				return;
			}
			//////////

			///////////////////////////////
			// COLOR
			////////////////////////////////

			
			//numCellsToMove=trafficPersonVec[p].v;
			trafficPersonVec[p].posInLaneM=trafficPersonVec[p].posInLaneM+numMToMove;
			if(trafficPersonVec[p].posInLaneM>trafficPersonVec[p].length){//research intersection
				numMToMove=trafficPersonVec[p].posInLaneM-trafficPersonVec[p].length;
				getToNextEdge=true;
			}else{//does not research next intersection
				////////////////////////////////////////////////////////
				// LANE CHANGING (happens when we are not reached the intersection)
				if(trafficPersonVec[p].v>3.0f&&//at least 10km/h to try to change lane
					trafficPersonVec[p].num_steps%10==0//just check every (10 steps) 5 seconds	
					){
					//next thing is not a traffic light 
					// skip if there is one lane (avoid to do this)
					// skip if it is the last edge
					if(nextVehicleIsATrafficLight==false&&trafficPersonVec[p].edgeNumLanes>1&&nextEdge!=0xFFFF){

						////////////////////////////////////////////////////
						// LC 1 update lane changing status
						if(trafficPersonVec[p].LC_stateofLaneChanging==0){
							// 2.2-exp((x-1)^2)
							float x=trafficPersonVec[p].posInLaneM/trafficPersonVec[p].length;
							if(x>0.4f){//just after 40% of the road
								float probabiltyMandatoryState=2.2-exp((x-1)*(x-1));
								{
									trafficPersonVec[p].LC_stateofLaneChanging=1;
								}
							}

						}

						////////////////////////////////////////////////////
						// LC 2 NOT MANDATORY STATE
						if(trafficPersonVec[p].LC_stateofLaneChanging==0){
							//if(p==40)printf("LC v %f v0 %f a %f\n",trafficPersonVec[p].v,trafficPersonVec[p].maxSpeedMperSec*0.5f,dv_dt);
							// discretionary change: v slower than the current road limit and deccelerating and moving
							if((trafficPersonVec[p].v<(trafficPersonVec[p].maxSpeedMperSec*0.7f))&&(dv_dt<0)&&trafficPersonVec[p].v>3.0f){
								//printf(">>LANE CHANGE\n");
								ushort laneToCheck;//!!!!
								//printf("LC 0 %u\n",trafficPersonVec[p].numOfLaneInEdge);
								bool leftLane=trafficPersonVec[p].numOfLaneInEdge>0;//at least one lane on the left
								bool rightLane=trafficPersonVec[p].numOfLaneInEdge<trafficPersonVec[p].edgeNumLanes-1;//at least one lane
								if(leftLane==true&&rightLane==true){
									if((int(currentTime))%2==0){
										leftLane=false;
									}else{
										rightLane=false;
									}
								}
								if(leftLane==true){
									laneToCheck=trafficPersonVec[p].numOfLaneInEdge-1;
								}else{
									laneToCheck=trafficPersonVec[p].numOfLaneInEdge+1;
								}

								uchar v_a,v_b;float gap_a,gap_b;
								//printf("p %u LC 1 %u\n",p,laneToCheck);
								uchar trafficLightState=trafficLights[currentEdge];
								calculateGapsLC(cellSize,maxWidth,mapToReadShift,laneMap,trafficLightState,currentEdge+laneToCheck,trafficPersonVec[p].posInLaneM,trafficPersonVec[p].length,v_a,v_b,gap_a,gap_b);
								//printf("LC 2 %u %u %f %f\n",v_a,v_b,gap_a,gap_b);
								if(gap_a==1000.0f&&gap_b==1000.0f){//lag and lead car very far
									trafficPersonVec[p].numOfLaneInEdge=laneToCheck;// CHANGE LINE

								}else{// NOT ALONE
									float b1A=0.05f,b2A=0.15f;
									float b1B=0.15f,b2B=0.40f;
									// s_0-> critical lead gap
									float g_na_D,g_bn_D;
									bool acceptLC=true;
									if(gap_a!=1000.0f){
										g_na_D=fmax(s_0,s_0+b1A*trafficPersonVec[p].v+b2A*(trafficPersonVec[p].v-v_a*3.0f));
										if(gap_a<g_na_D)//gap smaller than critical gap
											acceptLC=false;
									}
									if(acceptLC==true&&gap_b!=1000.0f){
										g_bn_D=fmax(s_0,s_0+b1B*v_b*3.0f+b2B*(v_b*3.0f-trafficPersonVec[p].v));
										if(gap_b<g_bn_D)//gap smaller than critical gap
											acceptLC=false;
									}
									if(acceptLC==true){
										trafficPersonVec[p].numOfLaneInEdge=laneToCheck;// CHANGE LINE
									}
								}
								//printf("<<LANE CHANGE\n");
							}


						}// Discretionary

						////////////////////////////////////////////////////
						// LC 3 *MANDATORY* STATE
						if(trafficPersonVec[p].LC_stateofLaneChanging==1){
							// LC 3.1 Calculate the correct lanes
							if(trafficPersonVec[p].LC_endOKLanes==0xFF){
								calculateLaneCarShouldBe(currentEdge,nextEdge,intersections,trafficPersonVec[p].edgeNextInters,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
								//printf("p%u num lanes %u min %u max %u\n",p,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
								if(trafficPersonVec[p].LC_initOKLanes==0&&trafficPersonVec[p].LC_endOKLanes==0)
									return;
							}


							//printf(">>LANE CHANGE\n");
							ushort laneToCheck;//!!!!
							//printf("LC 0 %u\n",trafficPersonVec[p].numOfLaneInEdge);
							bool leftLane=false,rightLane=false;
							// LC 3.2 CORRECT LANES--> DICRETIONARY LC WITHIN
							if(trafficPersonVec[p].numOfLaneInEdge>=trafficPersonVec[p].LC_initOKLanes&&trafficPersonVec[p].numOfLaneInEdge<trafficPersonVec[p].LC_endOKLanes){
								// for discretionary it should be under some circustances
								if((trafficPersonVec[p].v<(trafficPersonVec[p].maxSpeedMperSec*0.7f))&&(dv_dt<0)&&trafficPersonVec[p].v>3.0f){
									leftLane=
										(trafficPersonVec[p].numOfLaneInEdge>0) &&//at least one lane on the left
										(trafficPersonVec[p].numOfLaneInEdge-1>=trafficPersonVec[p].LC_initOKLanes)&&
										(trafficPersonVec[p].numOfLaneInEdge-1<trafficPersonVec[p].LC_endOKLanes);
									rightLane=
										(trafficPersonVec[p].numOfLaneInEdge<trafficPersonVec[p].edgeNumLanes-1)&&//at least one lane
										(trafficPersonVec[p].numOfLaneInEdge+1>=trafficPersonVec[p].LC_initOKLanes)&&
										(trafficPersonVec[p].numOfLaneInEdge+1<trafficPersonVec[p].LC_endOKLanes);
									//printf("D\n");
								}
							}
							// LC 3.3 INCORRECT LANES--> MANDATORY LC
							else{
								//printf("num lanes %u min %u max %u\n",trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);
								//printf("p%u num lanes %u min %u max %u\n",p,trafficPersonVec[p].edgeNumLanes,trafficPersonVec[p].LC_initOKLanes,trafficPersonVec[p].LC_endOKLanes);

								if(trafficPersonVec[p].numOfLaneInEdge<trafficPersonVec[p].LC_initOKLanes){
									rightLane=true;
								}else{
									leftLane=true;
								}
								if(rightLane==true&&trafficPersonVec[p].numOfLaneInEdge+1>=trafficPersonVec[p].edgeNumLanes){
									
								}
								if(leftLane==true&&trafficPersonVec[p].numOfLaneInEdge==0){
									
									return;
								}
								//printf("M L %d R %d nL %u\n",leftLane,rightLane,trafficPersonVec[p].numOfLaneInEdge);
							}
							if(leftLane==true||rightLane==true){

								// choose lane (if necessary)
								if(leftLane==true&&rightLane==true){
									if((int(currentTime))%2==0){
										leftLane=false;
									}else{
										rightLane=false;
									}
								}
								if(leftLane==true){
									laneToCheck=trafficPersonVec[p].numOfLaneInEdge-1;
								}else{
									laneToCheck=trafficPersonVec[p].numOfLaneInEdge+1;
								}
								if(laneToCheck>=trafficPersonVec[p].edgeNumLanes){
									
								}
								uchar v_a,v_b;float gap_a,gap_b;
								//printf("p %u LC 1 %u\n",p,laneToCheck);
								uchar trafficLightState=trafficLights[currentEdge];
								calculateGapsLC(cellSize,maxWidth,mapToReadShift,laneMap,trafficLightState,currentEdge+laneToCheck,trafficPersonVec[p].posInLaneM,trafficPersonVec[p].length,v_a,v_b,gap_a,gap_b);
								//printf("LC 2 %u %u %f %f\n",v_a,v_b,gap_a,gap_b);
								if(gap_a==1000.0f&&gap_b==1000.0f){//lag and lead car very far
									trafficPersonVec[p].numOfLaneInEdge=laneToCheck;// CHANGE LINE

								}else{// NOT ALONE
									float b1A=0.05f,b2A=0.15f;
									float b1B=0.15f,b2B=0.40f;
									float gamma=0.000025;
									// s_0-> critical lead gap
									float distEnd=trafficPersonVec[p].length-trafficPersonVec[p].posInLaneM;
									float expTerm=(1-exp(-gamma*distEnd*distEnd));

									float g_na_M,g_bn_M;
									bool acceptLC=true;
									if(gap_a!=1000.0f){
										g_na_M=fmax(s_0,s_0+(b1A*trafficPersonVec[p].v+b2A*(trafficPersonVec[p].v-v_a*3.0f)));
										if(gap_a<g_na_M)//gap smaller than critical gap
											acceptLC=false;
									}
									if(acceptLC==true&&gap_b!=1000.0f){
										g_bn_M=fmax(s_0,s_0+(b1B*v_b*3.0f+b2B*(v_b*3.0f-trafficPersonVec[p].v)));
										if(gap_b<g_bn_M)//gap smaller than critical gap
											acceptLC=false;
									}
									if(acceptLC==true){
										trafficPersonVec[p].numOfLaneInEdge=laneToCheck;// CHANGE LINE
									}
								}


							}

						}// Mandatory

					}//at least two lanes and not stopped by traffic light

				}
				///////////////////////////////////////////////////////
				
				uchar vInMpS=(uchar)(trafficPersonVec[p].v*3);//speed in m/s to fit in uchar
				ushort posInLineCells=(ushort)(trafficPersonVec[p].posInLaneM/cellSize);
				laneMap[mapToWriteShift+maxWidth*(currentEdge+trafficPersonVec[p].numOfLaneInEdge)+posInLineCells]=vInMpS;
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
			if(nextEdge==0xFFFF){//if(curr_intersection==end_intersection){
				
				trafficPersonVec[p].active=2;//finished
				return;
			}
			//if(trafficPersonVec[p].nextPathEdge>=nextEdgeM.size())printf("AAAAAAAAAAAAAAAAA\n");
			/////////////
			// update edge
			//trafficPersonVec[p].curEdgeLane=trafficPersonVec[p].nextEdge;
			trafficPersonVec[p].currPathEdge++;
			trafficPersonVec[p].maxSpeedMperSec=trafficPersonVec[p].nextEdgemaxSpeedMperSec;
			trafficPersonVec[p].edgeNumLanes=trafficPersonVec[p].nextEdgeNumLanes;
			trafficPersonVec[p].edgeNextInters=trafficPersonVec[p].nextEdgeNextInters;
			trafficPersonVec[p].length=trafficPersonVec[p].nextEdgeLength;
			trafficPersonVec[p].posInLaneM=numMToMove;
			if(trafficPersonVec[p].numOfLaneInEdge>=trafficPersonVec[p].edgeNumLanes){
				trafficPersonVec[p].numOfLaneInEdge=trafficPersonVec[p].edgeNumLanes-1;//change line if there are less roads
			}

			////////////
			// update next edge
			ushort nextNEdge=trafficPersonVec[p].personPath[trafficPersonVec[p].currPathEdge+1];
			//trafficPersonVec[p].nextEdge=nextEdge;
			if(nextNEdge!=0xFFFF){
				//trafficPersonVec[p].nextPathEdge++;
				trafficPersonVec[p].LC_initOKLanes=0xFF;
				trafficPersonVec[p].LC_endOKLanes=0xFF;
				
				//2.2.3 update person edgeData
				//trafficPersonVec[p].nextEdge=nextEdge;
				trafficPersonVec[p].nextEdgemaxSpeedMperSec=edgesData[nextNEdge].maxSpeedMperSec;
				trafficPersonVec[p].nextEdgeNumLanes=edgesData[nextNEdge].numLines;
				trafficPersonVec[p].nextEdgeNextInters=edgesData[nextNEdge].nextInters;
				trafficPersonVec[p].nextEdgeLength=edgesData[nextNEdge].length;
			}

			trafficPersonVec[p].LC_stateofLaneChanging=0;
			uchar vInMpS=(uchar)(trafficPersonVec[p].v*3);//speed in m/s to fit in uchar
			ushort posInLineCells=(ushort)(trafficPersonVec[p].posInLaneM/cellSize);
			laneMap[mapToWriteShift+maxWidth*(nextEdge+trafficPersonVec[p].numOfLaneInEdge)+posInLineCells]=vInMpS;

		 ///
	 }
	 
}//

__global__ void kernel_intersectionSimulation(uint numIntersections,float currentTime,LC::intersectionData *intersections,uchar *trafficLights) {
	 int i = blockIdx.x * blockDim.x + threadIdx.x;
	 if(i<numIntersections){//CUDA check (inside margins)

		 //////////////////////////////////////////////////////
		 const float deltaEvent=20.0f;
		 //if(i==0)printf("i %d\n",i);
		 if(currentTime>intersections[i].nextEvent&&intersections[i].totalInOutEdges>0){


			 uint edgeOT=intersections[i].edge[intersections[i].state];
			 uchar numLinesO=edgeOT>>24;
			 ushort edgeONum=edgeOT&0xFFFF;
			 // red old traffic lights
			 for(int nL=0;nL<numLinesO;nL++){
				 trafficLights[edgeONum+nL]=0x00;//red old traffic light
			 }
			 for(int iN=0;iN<=intersections[i].totalInOutEdges;iN++){//to give a round
				 intersections[i].state=(intersections[i].state+1)%intersections[i].totalInOutEdges;//next light
				 if((intersections[i].edge[intersections[i].state]&0x010000)==0x010000){
					 // green new traffic lights
					 uint edgeIT=intersections[i].edge[intersections[i].state];
					 ushort edgeINum=edgeIT&0xFFFF;//get edgeI
					 uchar numLinesI=edgeIT>>24;
					 for(int nL=0;nL<numLinesI;nL++){
						 trafficLights[edgeINum+nL]=0xFF;
					 }
					 break;
				 }
			 }//green new traffic light
			 //printf("i %d CHANGE state %u of %d (Old edge %u New Edge %u)\n",i,intersections[i].state,intersections[i].totalInOutEdges,edgeO,edgeI);
			 ////
			 intersections[i].nextEvent=currentTime+deltaEvent;
		 }
	 

		 //////////////////////////////////////////////////////
	 }
	 
 }//


void simulateTrafficCUDA(float currentTime,uint numPeople,uint numIntersections){

	////////////////////////////////////////////////////////////
	// 1. CHANGE MAP: set map to use and clean the other
	if(readFirstMapC==true){
		mapToReadShift=0;
		mapToWriteShift=halfLaneMap;
		cudaMemset(&laneMap_d[halfLaneMap],-1,halfLaneMap*sizeof (unsigned char));//clean second half
	}else{
		mapToReadShift=halfLaneMap;
		mapToWriteShift=0;
		cudaMemset(&laneMap_d[0],-1,halfLaneMap*sizeof (unsigned char));//clean first half
	}
	readFirstMapC=!readFirstMapC;//next iteration invert use

	//float currentTime=3.0f;
	//printf("!! >>Call kernel\n");
	//while(currentTime<24.0f){

	kernel_intersectionSimulation <<< numIntersections,1 >>> (numIntersections, currentTime,intersections_d,trafficLights_d);
	//int blocks=1;
	//if(numPeople>10000){
		//blocks=1024;
		//numPeople=ceil(numPeople/1024.0f);
	//}
	kernel_trafficSimulation <<< ceil(numPeople/1024.0f),1024 >>> (numPeople,currentTime,cellSize,deltaTime,trafficPersonVec_d,edgesData_d,laneMap_d,intersections_d,trafficLights_d, mapToReadShift,mapToWriteShift,maxWidth);// n_blocks, block_size >>> (a_d, N);
		//currentTime+=deltaTimeC;
	//}
	//printf("!! <<End kernel\n");

	

	//exampleCUDA();
	//const float deltaTime=1.0f*(1.0f/3600.0f);//1*(1hours/3600 seconds)
	//float currentTime=3.0f;//lets start at 300 AM
	// 0. put data as references to simplified code

	/*
	int numVertex=cudaTrafficRoutes.numVertex;
	int numPeople=trafficPersonVec.size();

	printf(">>Start Simulation\n");
	QTime timer;
	timer.start();
	while(currentTime<24.0f){
		bool anyActive=false;
		//printf("time %f\n",currentTime);
		for(int p=0;p<numPeople;p++){
			//printf(" 0. Person: %d\n",p);
			//0. check if finished
			

		}
		currentTime+=deltaTime;
		if(anyActive){
			clientMain->mGLWidget_3D->updateGL();
			Sleep(1);
			QApplication::processEvents();
		}
	}*/
	

	//printf("<<simulate\n");

}//
