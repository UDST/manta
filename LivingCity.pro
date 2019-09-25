QT += core

TARGET = LivingCity.bin
# Project build directories
DESTDIR     = $$PWD
OBJECTS_DIR = $$DESTDIR/obj
CONFIG += c++17

unix {
    LIBS += -L/opt/local/lib -lopencv_imgcodecs -lopencv_core -lopencv_imgproc
	# -L/Developer/NVIDIA/CUDA-7.5/lib -lcudart -lcublas -lgomp
    INCLUDEPATH += \
      /usr/include/opencv2/ \
      /opt/local/include/ \
      /usr/local/boost_1_59_0/ \
      $$PWD/glew/include/

    CONFIG += debug
}
win32{
    # Note: OpenCV uses 2.4.12 since I compiled with VS 2013 (vc12)
    LIBS+= \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_core2412.lib \ # e.g., OPENCV_BUILD=D:\opencv\build
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_imgproc2412.lib \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_highgui2412.lib \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_legacy2412.lib \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_ml2412.lib \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_photo2412.lib \
        $$(OPENCV_BUILD)/x64/vc12/lib/opencv_video2412.lib \


    INCLUDEPATH += \
        $$(OPENCV_BUILD)/include/opencv \ # e.g., OPENCV_BUILD=D:\opencv\build
        $$(OPENCV_BUILD)/include/opencv2 \
        $$(OPENCV_BUILD)/include \
        $$(BOOST_ROOT) \ # e.g., BOOST_ROOT = D:\boost\boost_1_59_0
        $$PWD/glew/include/

    CONFIG += console # show printf in terminal
}

HEADERS += \
    ./LivingCity/dataExporter.h \
    ./LivingCity/Geometry/block.h \
    ./LivingCity/Geometry/building.h \
    ./LivingCity/Geometry/client_geometry.h \
    ./LivingCity/Geometry/parcel.h \
    ./LivingCity/Geometry/parcelBuildingAttributes.h \
    ./LivingCity/Geometry/placeTypeInstances.h \
    ./LivingCity/Geometry/zone.h \
    ./LivingCity/RoadGraph/roadGraph.h \
    ./LivingCity/RoadGraph/roadGraphEdge.h \
    ./LivingCity/RoadGraph/roadGraphVertex.h \
    ./LivingCity/global.h \
    ./LivingCity/misctools/bounding_box.h \
    ./LivingCity/misctools/common.h \
    ./LivingCity/misctools/misctools.h \
    ./LivingCity/misctools/polygon_3D.h \
    ./LivingCity/roadGraphB2018Loader.h \
    ./LivingCity/traffic/b18CUDA_trafficSimulator.h \
    ./LivingCity/traffic/b18CommandLineVersion.h \
    ./LivingCity/traffic/b18EdgeData.h \
    ./LivingCity/traffic/b18GridPollution.h \
    ./LivingCity/traffic/b18TrafficDijkstra.h \
    ./LivingCity/traffic/b18TrafficJohnson.h \
    ./LivingCity/traffic/b18TrafficLaneMap.h \
    ./LivingCity/traffic/b18TrafficOD.h \
    ./LivingCity/traffic/b18TrafficPerson.h \
    ./LivingCity/traffic/b18TrafficSP.h \
    ./LivingCity/traffic/b18TrafficSimulator.h \
    ./LivingCity/traffic/boostGeometry.h \
    ./LivingCity/traffic/laneCoordinatesComputer.h \
    ./LivingCity/traffic/simulatorConfiguration.h \
    ./LivingCity/traffic/sp/config.h \
    ./LivingCity/traffic/sp/external/csv.h \
    ./LivingCity/traffic/sp/external/tsl/robin_growth_policy.h \
    ./LivingCity/traffic/sp/external/tsl/robin_hash.h \
    ./LivingCity/traffic/sp/external/tsl/robin_map.h \
    ./LivingCity/traffic/sp/external/tsl/robin_set.h \
    ./LivingCity/traffic/sp/graph.h \
    ./LivingCity/traffic/sp/mpi_wrapper.h \
    ./LivingCity/traffic/sp/unordered_map_tuple_hash.h \
    ./LivingCity/OSMConstants.h \
    ./LivingCity/trafficControl.h \

SOURCES += \
    ./LivingCity/dataExporter.cpp \
    ./LivingCity/Geometry/block.cpp \
    ./LivingCity/Geometry/building.cpp \
    ./LivingCity/Geometry/client_geometry.cpp \
    ./LivingCity/Geometry/parcel.cpp \
    ./LivingCity/Geometry/parcelBuildingAttributes.cpp \
    ./LivingCity/Geometry/placeTypeInstances.cpp \
    ./LivingCity/Geometry/zone.cpp \
    ./LivingCity/LC_main.cpp \
    ./LivingCity/RoadGraph/roadGraph.cpp \
    ./LivingCity/RoadGraph/roadGraphEdge.cpp \
    ./LivingCity/RoadGraph/roadGraphVertex.cpp \
    ./LivingCity/global.cpp \
    ./LivingCity/misctools/bounding_box.cpp \
    ./LivingCity/misctools/misctools.cpp \
    ./LivingCity/misctools/polygon_3D.cpp \
    ./LivingCity/roadGraphB2018Loader.cpp \
    ./LivingCity/traffic/b18CommandLineVersion.cpp \
    ./LivingCity/traffic/b18GridPollution.cpp \
    ./LivingCity/traffic/b18TrafficDijkstra.cpp \
    ./LivingCity/traffic/b18TrafficJohnson.cpp \
    ./LivingCity/traffic/b18TrafficLaneMap.cpp \
    ./LivingCity/traffic/b18TrafficOD.cpp \
    ./LivingCity/traffic/b18TrafficSP.cpp \
    ./LivingCity/traffic/b18TrafficSimulator.cpp \
    ./LivingCity/traffic/boostGeometry.cpp \
    ./LivingCity/traffic/laneCoordinatesComputer.cpp \
    ./LivingCity/traffic/simulatorConfiguration.cpp \
    ./LivingCity/traffic/sp/graph.cc \
    ./LivingCity/OSMConstants.cpp \
    ./LivingCity/trafficControl.cpp \

OTHER_FILES += \
        ./LivingCity/traffic/b18CUDA_trafficSimulator.cu \


###################################################################
## CUDA
###################################################################
win32{
    # Cuda sources
    CUDA_SOURCES += ./LivingCity/traffic/b18CUDA_trafficSimulator.cu

    # Path to cuda toolkit install
    CUDA_DIR      = "D:/CUDA"

    # Path to header and libs files
    INCLUDEPATH  += $$CUDA_DIR/include
    QMAKE_LIBDIR += $$CUDA_DIR/lib/x64

    SYSTEM_TYPE = 64            # '32' or '64', depending on your system

    # libs used in your code
    LIBS += -lcuda -lcudart
    CUDA_LIBS += -lcuda -lcudart # LIBS

    # GPU architecture
    CUDA_ARCH     = sm_50

    # Here are some NVCC flags I've always used by default.
    NVCCFLAGS     = --use_fast_math


    # Prepare the extra compiler configuration (taken from the nvidia forum - i'm not an expert in this part)
    CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')


    # MSVCRT link option (static or dynamic, it must be the same with your Qt SDK link option)
    MSVCRT_LINK_FLAG_DEBUG = "/MDd"
    MSVCRT_LINK_FLAG_RELEASE = "/MD"

    QMAKE_EXTRA_COMPILERS += cuda

    # Configuration of the Cuda compiler
    CONFIG(debug, debug|release) {
        # Debug mode
        cuda_d.input = CUDA_SOURCES
        cuda_d.output = $$OBJECTS_DIR/${QMAKE_FILE_BASE}.obj
        cuda_d.commands = $$CUDA_DIR/bin/nvcc.exe -D_DEBUG $$NVCC_OPTIONS $$CUDA_INC $$CUDA_LIBS --machine $$SYSTEM_TYPE \
                         -arch=$$CUDA_ARCH -c -Xcompiler $$MSVCRT_LINK_FLAG_DEBUG -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda_d.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda_d
    }
    else {
        # Release mode
        cuda.input = CUDA_SOURCES
        cuda.output = $$OBJECTS_DIR/${QMAKE_FILE_BASE}.obj
        cuda.commands = $$CUDA_DIR/bin/nvcc.exe $$NVCC_OPTIONS $$CUDA_INC $$CUDA_LIBS --machine $$SYSTEM_TYPE \
                       -arch=$$CUDA_ARCH -c -Xcompiler $$MSVCRT_LINK_FLAG_RELEASE -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda
    }
}

unix {
  # Cuda sources
  CUDA_SOURCES += ./LivingCity/traffic/b18CUDA_trafficSimulator.cu
  # Path to cuda toolkit install
  CUDA_DIR = /usr/local/cuda-9.0
  INCLUDEPATH += $$CUDA_DIR/include
  QMAKE_LIBDIR += $$CUDA_DIR/lib64
  # GPU architecture
  CUDA_ARCH = sm_50
  NVCCFLAGS = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v -Xcompiler -fopenmp
  # Path to libraries
  LIBS += -lcudart -lcuda -lgomp
  QMAKE_CXXFLAGS += -fopenmp
  #LIBS += -fopenmp
  # join the includes in a line
  CUDA_INC = $$join(INCLUDEPATH,' -I','-I',' ')
  cuda.commands = $$CUDA_DIR/bin/nvcc -m64 -O3 -arch=$$CUDA_ARCH -c $$NVCCFLAGS $$CUDA_INC $$LIBS ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}
  cuda.dependcy_type = TYPE_C
  cuda.depend_command = $$CUDA_DIR/bin/nvcc -O3 -M $$CUDA_INC $$NVCCFLAGS      ${QMAKE_FILE_NAME}

  cuda.input = CUDA_SOURCES
  cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o
  # Tell Qt that we want add more stuff to the Makefile
  QMAKE_EXTRA_COMPILERS += cuda
}

INCLUDEPATH += $$CUDA_DIR/include ./LivingCity . ./src

QMAKE_CXXFLAGS += -Wno-sign-compare -Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-parameter
