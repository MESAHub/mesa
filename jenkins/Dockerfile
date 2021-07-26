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
      ruby \
      && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
RUN curl -so /tmp/mesasdk.tar.gz http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-21.4.1.tar.gz && \
    tar xzf /tmp/mesasdk.tar.gz -C /opt/ && rm /tmp/mesasdk.tar.gz
ENV MESASDK_ROOT=/opt/mesasdk
RUN bash $MESASDK_ROOT/bin/mesasdk_init.sh
RUN gem install mesa_test
