FROM ubuntu:14.04
MAINTAINER Berke Toptas <berke.toptas@sbgdinc.com>

# Install build essential and g++
RUN apt-get update \
 && apt-get install -y --force-yes --no-install-recommends\
      groff \
      g++ \
      wget \
      build-essential \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
 && rm -rf /var/lib/apt/lists/*;

# Add htslib-1.6
WORKDIR /home
RUN wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 --no-check-certificate
RUN tar xvjf htslib-1.6.tar.bz2
WORKDIR /home/htslib-1.6

# Compile htslib-1.6
RUN ./configure
RUN make
RUN make install

# Add vbt folder
ADD vbt /home/varbenchtools
WORKDIR /home/varbenchtools

RUN cp -R /home/htslib-1.6/htslib /home/varbenchtools/
RUN cp /home/htslib-1.6/libhts.a /home/varbenchtools/lib/
RUN cp /home/htslib-1.6/libhts.so /home/varbenchtools/lib/
RUN cp /home/htslib-1.6/libhts.so.2 /home/varbenchtools/lib/

# Compile vbt
RUN make

RUN ldconfig

# Test if vbt is working
#RUN ./vbt varcomp --help