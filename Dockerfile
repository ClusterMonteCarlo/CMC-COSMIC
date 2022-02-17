FROM ubuntu:hirsute-20220113 AS spython-base
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get -y install libopenmpi-dev openmpi-bin libhdf5-serial-dev cmake git python3-mpi4py python3-pip python3-numpy libgsl-dev libgslcblas0
RUN git clone --recurse-submodules https://github.com/ClusterMonteCarlo/CMC-COSMIC.git /CMC-COSMIC
RUN cd /CMC-COSMIC && mkdir build && cd build && FC=mpifort CC=mpicc cmake .. -DBUILD_COSMIC=ON -DCMAKE_INSTALL_PREFIX=/opt/local/ && make install
CMD /opt/local/bin/cmc
