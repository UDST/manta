#FROM nvidia/cuda:11.2.1-devel

# Using 18.04 for compatibility with the qmake Makefile
FROM  nvidia/cuda:11.2.0-devel-ubuntu18.04

# Workarond so the setup doesnt get stuck
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# libraries
RUN apt-get update
RUN apt update
RUN apt-get install qtchooser -y
RUN apt-get install qt5-default -y
RUN apt-get install libglew-dev -y
RUN apt-get install build-essential -y
RUN apt-get install libfontconfig1 -y
RUN apt-get install mesa-common-dev -y
RUN apt-get install wget -y
RUN apt-get install pciutils -y
RUN apt-get install git -y
RUN apt-get install vim -y

# boost
RUN apt-get install wget
RUN wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz
RUN tar xf boost_1_59_0.tar.gz -C /usr/local


# Pandana installation
RUN cd /usr/include/ && \
    git clone https://github.com/UDST/pandana && \
    cd pandana/src && \
    echo "CC = gcc  # C compiler" > Makefile && \
    echo "CXX = g++" >> Makefile && \
    echo "CPPFLAGS = -DLINUX -DMAC -std=c++0x -c -fPIC -g -O3 -Wall -pedantic -fopenmp -w # C flags" >> Makefile && \
    echo "LDFLAGS = -shared   # linking flags" >> Makefile && \
    echo "RM = rm -f   # rm command" >> Makefile && \
    echo "TARGET_LIB = libchrouting.so  # target lib" >> Makefile && \ 
    echo "SRCS =  accessibility.cpp graphalg.cpp contraction_hierarchies/src/libch.cpp" >> Makefile && \
    echo "OBJS = \$(SRCS:.cpp=.o)" >> Makefile && \
    echo ".PHONY: all" >> Makefile && \
    echo "all: \${TARGET_LIB}" >> Makefile && \
    echo "\$(TARGET_LIB): \$(OBJS)" >> Makefile && \
    echo "\t\$(CXX) \${LDFLAGS} -o \$@ $^" >> Makefile && \
    echo ".PHONY: clean" >> Makefile && \
    echo "clean:" >> Makefile && \
    echo "\t-\${RM} \${TARGET_LIB} \${OBJS}" >> Makefile && \
    make
     

# OpenCV dependencies
RUN apt -y install build-essential checkinstall cmake pkg-config yasm
RUN apt -y install git gfortran
RUN apt -y install libjpeg8-dev libpng-dev

RUN apt -y install software-properties-common
RUN add-apt-repository "deb http://security.ubuntu.com/ubuntu xenial-security main"
RUN apt -y update

RUN apt -y install libjasper1
RUN apt -y install libtiff-dev

RUN apt -y install libavcodec-dev libavformat-dev libswscale-dev libdc1394-22-dev
RUN apt -y install libxine2-dev libv4l-dev
RUN cd /usr/include/linux
RUN ln -s -f ../libv4l1-videodev.h videodev.h

RUN apt -y install libgstreamer1.0-dev libgstreamer-plugins-base1.0-dev
RUN apt -y install libgtk2.0-dev libtbb-dev qt5-default
RUN apt -y install libatlas-base-dev
RUN apt -y install libfaac-dev libmp3lame-dev libtheora-dev
RUN apt -y install libvorbis-dev libxvidcore-dev
RUN apt -y install libopencore-amrnb-dev libopencore-amrwb-dev
RUN apt -y install libavresample-dev

RUN apt -y install unzip


# OpenCV installation
RUN git clone https://github.com/opencv/opencv.git && \
    cd opencv && \
    mkdir build && \
    cd build && \
    cmake  \
    -D BUILD_TIFF=ON  \
    -D WITH_CUDA=OFF \
    -D ENABLE_AVX=OFF \
    -D WITH_OPENGL=OFF \
    -D WITH_OPENCL=OFF \
    -D WITH_IPP=OFF \
    -D WITH_TBB=ON \
    -D BUILD_TBB=ON \
    -D WITH_EIGEN=OFF \
    -D WITH_V4L=OFF \
    -D WITH_VTK=OFF \
    -D BUILD_TESTS=OFF \
    -D BUILD_PERF_TESTS=OFF \
    -D OPENCV_GENERATE_PKGCONFIG=ON \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D CMAKE_INSTALL_PREFIX=/usr .. && \
#     cmake-D CMAKE_BUILD_TYPE=RELEASE \
#    -D CMAKE_INSTALL_PREFIX=/usr/local \
#     -D INSTALL_C_EXAMPLES=ON \
#     -D WITH_TBB=ON \
#     -D WITH_V4L=ON \
#     -D OPENCV_PYTHON3_INSTALL_PATH=$cwd/OpenCV-$cvVersion-py3/lib/python3.5/site-packages \
#     -D WITH_QT=ON \
#     -D WITH_OPENGL=ON \
#     -D BUILD_EXAMPLES=ON .. && \
     make -j 12 && \
     make install





# CUDA paths
ENV PATH="/usr/local/cuda-11.2/bin:${PATH}"
#RUN export PATH=/usr/local/cuda-11.2/bin:$PATH
#RUN export LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LIBRARY_PATH
ENV LIBRARY_PATH="/usr/local/cuda-11.2/lib64:${LIBRARY_PATH}"
#RUN export LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH="/usr/local/cuda-11.2/lib64:${LD_LIBRARY_PATH}"

# Pandana path - modify it in case Pandana is not in your home directory
#RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/include/pandana/src
ENV LD_LIBRARY_PATH="/usr/include/pandana/src:${LD_LIBRARY_PATH}"

# Python libraries
RUN apt install python3-pip

RUN pip3 install -r requirements.txt

# Check if CUDA is properly installed
CMD nvidia-smi
