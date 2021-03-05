# ========== OpenCV subimage ==========
FROM ubuntu:18.04 AS opencvbuilder

# OpenCV dependencies
RUN apt update
RUN apt -y install build-essential checkinstall cmake pkg-config yasm \
    git gfortran libjpeg8-dev libpng-dev

RUN apt -y install software-properties-common
RUN add-apt-repository "deb http://security.ubuntu.com/ubuntu xenial-security main"
RUN apt -y update

RUN apt -y install libjasper1 libtiff-dev \
    libavcodec-dev libavformat-dev libswscale-dev \
    libdc1394-22-dev libxine2-dev libv4l-dev

RUN cd /usr/include/linux
RUN ln -s -f ../libv4l1-videodev.h videodev.h

RUN apt -y install libgstreamer1.0-dev libgstreamer-plugins-base1.0-dev \
    libgtk2.0-dev libtbb-dev qt5-default \
    libatlas-base-dev \
    libfaac-dev libmp3lame-dev libtheora-dev \
    libvorbis-dev libxvidcore-dev \
    libopencore-amrnb-dev libopencore-amrwb-dev \
    libavresample-dev \
    unzip


# OpenCV compilation and installation
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
     make -j 12 && \
     make install


# ========== Pandana subimage ==========
FROM ubuntu:18.04 AS pandanabuilder

WORKDIR /usr/include/

RUN apt update && \
    apt install -y qtchooser \
    qt5-default \
    libglew-dev \
    build-essential \
    libfontconfig1 \
    mesa-common-dev \
    wget \
    pciutils \
    git

RUN git clone https://github.com/UDST/pandana

COPY /PandanaMakefile /usr/include/pandana/src/Makefile

WORKDIR /usr/include/pandana/src

RUN make

# ========== MANTA image ==========
FROM nvidia/cuda:11.2.0-devel-ubuntu18.04 as mantabuilder

COPY --from=opencvbuilder /usr/include/opencv4/ /usr/include/opencv4/

#COPY --from=opencvbuilder /opt/local/lib/ /opt/local/lib/

COPY --from=opencvbuilder /usr/lib/x86_64-linux-gnu/libopencv_core.* /usr/lib/x86_64-linux-gnu/
COPY --from=opencvbuilder /usr/lib/x86_64-linux-gnu/libopencv_imgproc.* /usr/lib/x86_64-linux-gnu/
COPY --from=opencvbuilder /usr/lib/x86_64-linux-gnu/libopencv_imgcodecs.* /usr/lib/x86_64-linux-gnu/

COPY --from=pandanabuilder /usr/include/pandana/ /usr/include/pandana/

# libraries
RUN apt update && \
    apt install qtchooser \
    qt5-default \
    libglew-dev \
    build-essential \
    libfontconfig1 \
    mesa-common-dev \
    wget \
    pciutils -y

# boost
RUN wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz && \
    tar xf boost_1_59_0.tar.gz -C /usr/local

# CUDA paths
ENV PATH="/usr/local/cuda-11.2/bin:${PATH}"
ENV LIBRARY_PATH="/usr/local/cuda-11.2/lib64:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="/usr/local/cuda-11.2/lib64:${LD_LIBRARY_PATH}"

# Pandana path - modify it in case Pandana is not in your home directory
ENV LD_LIBRARY_PATH="/usr/include/pandana/src:${LD_LIBRARY_PATH}"

# Python libraries
RUN apt install python3-pip -y

ADD . ./

RUN pip3 install -r requirements.txt


# Check if CUDA is properly installed
CMD nvidia-smi
