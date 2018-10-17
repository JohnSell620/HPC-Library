FROM jsellers/gpu-mpi:onbuild

USER root

#### Install Google Test ####
RUN apt-get install -y libgtest-dev && \
    apt-get install -y cmake && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make && \
    cp *.a /usr/lib

# Copy the content of `HPCLibrary` directory in the host machine to
# the current working directory in this Docker image
COPY HPCcode/ .

# Build HPC code
RUN make

USER ${USER}
