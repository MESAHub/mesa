# http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Linux_.01Intel.01
FROM ubuntu:focal
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
RUN curl -L https://zenodo.org/records/13768913/files/mesasdk-x86_64-linux-24.7.1.tar.gz | \
    tar xzf - -C /opt/
ENV MESASDK_ROOT=/opt/mesasdk
RUN bash $MESASDK_ROOT/bin/mesasdk_init.sh
RUN gem install mesa_test -v 1.1.11
