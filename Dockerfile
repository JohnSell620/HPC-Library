FROM jhnns/ubuntu-cuda-mpich:onbuild

USER root

#### Copy HPC-Library and make MPI_EXECS ####
COPY HPCLibrary .
RUN rm -rf dev graphs && make
RUN echo "export PATH=/HPCLibrary/exe:$PATH" >> /etc/profile

USER ${USER}
