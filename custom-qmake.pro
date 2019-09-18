TARGET = LivingCity.bin

# Define output directories
DESTDIR = ./
OBJECTS_DIR = release/obj
CUDA_OBJECTS_DIR = release/cuda
CONFIG += c++17

HEADERS += \
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
    ./LivingCity/traffic/b18TrafficSimulator.h \

SOURCES += \
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
    ./LivingCity/traffic/b18TrafficSimulator.cpp \

CUDA_SOURCES += ./LivingCity/traffic/b18CUDA_trafficSimulator.cu

# This makes the .cu files appear in your project
OTHER_FILES += ./LivingCity/traffic/b18CUDA_trafficSimulator.cu

CUDA_DIR = /usr/local/cuda
SYSTEM_TYPE = 64
CUDA_ARCH = sm_50
NVCC_OPTIONS = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v

INCLUDEPATH += $$CUDA_DIR/include ./LivingCity . ./src 

LIBS += -lcuda -lcudart

# The following makes sure all path names (which often include spaces) are put between quotation marks
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')

cuda.input = CUDA_SOURCES
cuda.output = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.o
cuda.commands = $$CUDA_DIR/bin/nvcc $$NVCC_OPTIONS $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
cuda.dependency_type = TYPE_C
QMAKE_EXTRA_COMPILERS += cuda
