QT += core gui opengl

LIBS += -L/opt/local/lib -lopencv_imgcodecs -lopencv_core -lopencv_imgproc -lGLEW

INCLUDEPATH += /opt/local/include/GL/ /opt/local/include/

RESOURCES += \
    LC_UrbanMain_.qrc \
    LC_UrbanMain.qrc

FORMS += \
    LC_UrbanMain.ui

HEADERS += \
    global.h \
    LC_camera_3d.h \
    LC_GLWidget3D_Shadows.h \
    LC_GLWidget3D.h \
    LC_Layer.h \
    LC_UrbanMain.h \
    roadGraphDynameq.h \
    VBOBlocks.h \
    VBOBuilding.h \
    VBOGUI.h \
    VBOModel_StreetElements.h \
    VBOModel.h \
    VBOPeopleJobInfoLayer.h \
    VBORenderManager.h \
    VBORoadGraph.h \
    VBORoadLabels.h \
    VBOShader.h \
    VBOSkyBox.h \
    VBOTerrain.h \
    VBOText.h \
    VBOUtil.h \
    VBOVegetation.h \
    VBOWater.h \
    bTraffic/bCPUTrafficThread.h \
    bTraffic/bCUDA_trafficSimulator.h \
    bTraffic/bEdgeIntersectionData.h \
    bTraffic/bGenerateTest.h \
    bTraffic/bPMTrafficPerson.h \
    bTraffic/bTrafficConstants.h \
    bTraffic/bTrafficDijkstra.h \
    bTraffic/bTrafficIntersection.h \
    bTraffic/bTrafficJohnson.h \
    bTraffic/bTrafficLaneMap.h \
    bTraffic/bTrafficPeople.h \
    bTraffic/bTrafficSimulator.h \
    Geometry/block.h \
    Geometry/building.h \
    Geometry/client_geometry.h \
    Geometry/parcel.h \
    Geometry/parcelBuildingAttributes.h \
    Geometry/placeTypeInstances.h \
    Geometry/zone.h \
    misctools/bounding_box.h \
    misctools/common.h \
    misctools/misctools.h \
    misctools/polygon_3D.h \
    nvModel/nvMath.h \
    nvModel/nvMatrix.h \
    nvModel/nvModel.h \
    nvModel/nvQuaternion.h \
    nvModel/nvVector.h \
    PM/pmBlocks.h \
    PM/pmBuildings.h \
    PM/pmMain.h \
    PM/pmParcels.h \
    PM/pmRoads.h \
    RoadGraph/roadGraph.h \
    RoadGraph/roadGraphEdge.h \
    RoadGraph/roadGraphVertex.h \
    traffic/cudaEdgeData.h \
    traffic/cudaGridPollution.h \
    traffic/cudaPmTrafficPersonJob.h \
    traffic/cudaTrafficDesigner.h \
    traffic/cudaTrafficLaneMap.h \
    traffic/cudaTrafficMCMC.h \
    traffic/cudaTrafficPerson.h \
    traffic/cudaTrafficPersonPath.h \
    traffic/cudaTrafficPersonShortestPath.h \
    traffic/cudaTrafficRoutes.h \
    traffic/cudaTrafficSimulator.h \
    triangle/triangle.h

SOURCES += \
    global.cpp \
    LC_camera_3d.cpp \
    LC_GLWidget3D_Shadows.cpp \
    LC_GLWidget3D.cpp \
    LC_Layer.cpp \
    LC_main.cpp \
    LC_UrbanMain.cpp \
    roadGraphDynameq.cpp \
    VBOBlocks.cpp \
    VBOBuilding.cpp \
    VBOGUI.cpp \
    VBOModel_StreetElements.cpp \
    VBOModel.cpp \
    VBOPeopleJobInfoLayer.cpp \
    VBORenderManager.cpp \
    VBORoadGraph.cpp \
    VBORoadLabels.cpp \
    VBOShader.cpp \
    VBOSkyBox.cpp \
    VBOTerrain.cpp \
    VBOText.cpp \
    VBOUtil.cpp \
    VBOVegetation.cpp \
    VBOWater.cpp \
    bTraffic/bCPUTrafficThread.cpp \
    bTraffic/bGenerateTest.cpp \
    bTraffic/bPMTrafficPerson.cpp \
    bTraffic/bTrafficDijkstra.cpp \
    bTraffic/bTrafficIntersection.cpp \
    bTraffic/bTrafficJohnson.cpp \
    bTraffic/bTrafficLaneMap.cpp \
    bTraffic/bTrafficSimulator.cpp \
    Geometry/block.cpp \
    Geometry/building.cpp \
    Geometry/client_geometry.cpp \
    Geometry/parcel.cpp \
    Geometry/parcelBuildingAttributes.cpp \
    Geometry/placeTypeInstances.cpp \
    Geometry/zone.cpp \
    misctools/bounding_box.cpp \
    misctools/misctools.cpp \
    misctools/polygon_3D.cpp \
    nvModel/nvModel.cpp \
    nvModel/nvModelObj.cpp \
    nvModel/nvModelQuery.cpp \
    PM/pmBlocks.cpp \
    PM/pmBuildings.cpp \
    PM/pmMain.cpp \
    PM/pmParcels.cpp \
    PM/pmRoads.cpp \
    RoadGraph/roadGraph.cpp \
    RoadGraph/roadGraphEdge.cpp \
    RoadGraph/roadGraphVertex.cpp \
    traffic/cudaGridPollution.cpp \
    traffic/cudaPmTrafficPersonJob.cpp \
    traffic/cudaTrafficDesigner.cpp \
    traffic/cudaTrafficLaneMap.cpp \
    traffic/cudaTrafficMCMC.cpp \
    traffic/cudaTrafficPersonPath.cpp \
    traffic/cudaTrafficPersonShortestPath.cpp \
    traffic/cudaTrafficRoutes.cpp \
    traffic/cudaTrafficSimulator.cpp \
    triangle/triangle.c

DISTFILES += \
    bTraffic/bCUDA_trafficSimulator.cu
