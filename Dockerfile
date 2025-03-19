FROM ubuntu:focal

# set a directory for the app
WORKDIR /usr/mesa

# copy mesa files to the container
COPY . .

# prereqs for MESASDK
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      binutils \
      curl \
      make \
      perl \
      libx11-dev \
      zlib1g \
      zlib1g-dev \
      tcsh \
      git \
      git-lfs \
      python3 \
      ruby \
      && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
# MESASDK from Zenodo
RUN curl -L https://zenodo.org/records/13768913/files/mesasdk-x86_64-linux-24.7.1.tar.gz | \
    tar xzf - -C /opt/
ENV MESASDK_ROOT=/opt/mesasdk

# Need to set shell to bash for source to work as part of MESA install
SHELL ["/bin/bash", "-c"]

# set OMP_NUM_THREADS to be the number of cores on your machine
ENV OMP_NUM_THREADS=4

# set MESA_DIR to be the directory to which you downloaded MESA
ENV MESA_DIR=/usr/mesa

# add shmesa (the MESA command line tool) to your PATH
ENV PATH=$PATH:$MESA_DIR/scripts/shmesa

# GYRE 
ENV GYRE_DIR=$MESA_DIR/gyre/gyre

# Print info
RUN echo $MESASDK_ROOT
RUN echo $MESA_DIR
RUN echo $OMP_NUM_THREADS
RUN echo $PATH
RUN pwd
RUN ls

# Source MESASDK and install MESA
RUN cd $MESA_DIR && source $MESASDK_ROOT/bin/mesasdk_init.sh && ./install
CMD ["/bin/bash", "-c", "source /opt/mesasdk/bin/mesasdk_init.sh && exec /bin/bash"]
