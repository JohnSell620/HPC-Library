FROM jhnns/ubuntu-cuda-mpich:onbuild

USER root

#### Install Google Test ####
RUN apt-get update && apt-get install -y libgtest-dev cmake && \
    rm -rf /var/lib/apt/lists/* && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make && \
    cp *.a /usr/lib

# Copy the content of `HPCLibrary` directory in the host machine to
# the current working directory in this Docker image
COPY HPCLibrary .

# Build HPC code
RUN rm -rf dev graphs && make

USER ${USER}
