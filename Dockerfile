FROM jsellers/hpc

# Copy the `HPC-Library` directory from host to container
COPY HPC-Library/ .

# Build HPC-Library objects and executables
RUN make all
