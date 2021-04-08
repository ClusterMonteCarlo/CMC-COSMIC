FROM ubuntu:groovy-20210225 AS spython-base
RUN apt-get -y update && apt-get -y install libopenmpi-dev openmpi-bin libhdf5-serial-dev cmake git python3-mpi4py python3-pip python3-numpy libgsl-dev libgslcblas0
RUN git clone https://github.com/ClusterMonteCarlo/CMC-COSMIC.git /CMC-COSMIC
RUN cd /CMC-COSMIC && mkdir build && cd build && FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/ && make install
CMD /usr/local/bin/cmc
