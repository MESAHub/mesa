FROM --platform=linux/amd64 ubuntu:focal

# set a directory for the app
WORKDIR /usr/mesa

# need to set shell to bash for `source` to work as part of MESA install
SHELL ["/bin/bash", "-c"]

# Support to build on Mac silicon
RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install g++-x86-64-linux-gnu libc6-dev-amd64-cross
RUN DEBIAN_FRONTEND=noninteractive apt-get -y clean

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

# download MESASDK from Zenodo
RUN curl -L https://zenodo.org/records/13768913/files/mesasdk-x86_64-linux-24.7.1.tar.gz | \
    tar xzf - -C /opt/

# set environment variable for MESASDK
ENV MESASDK_ROOT=/opt/mesasdk

# copy mesa files to the container
COPY . .

# set MESA_DIR to be the directory to which you downloaded MESA
ENV MESA_DIR=/usr/mesa

# add shmesa (the MESA command line tool) to your PATH
ENV PATH=$PATH:$MESA_DIR/scripts/shmesa

# GYRE 
ENV GYRE_DIR=$MESA_DIR/gyre/gyre

# set environment variable to indicate we are in a docker container
ENV IN_DOCKER=true

# print info
RUN echo $MESASDK_ROOT
RUN echo $MESA_DIR
RUN echo $PATH
RUN pwd
RUN ls

# skip architecture check
RUN touch "${MESASDK_ROOT}/etc/check_arch.done"

# source MESASDK and install MESA
RUN cd $MESA_DIR && source $MESASDK_ROOT/bin/mesasdk_init.sh && ./install
CMD ["/bin/bash", "-c", "source /opt/mesasdk/bin/mesasdk_init.sh && exec /bin/bash"]
