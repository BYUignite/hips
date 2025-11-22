# Using the terminal. 
#   Install Docker desktop
#   Make sure docker executables are on the path (e.g., on mac: /Applications/Docker.app/Contents/Resources/bin)
# 
# Build the image: 
#     docker build -t hips_img .
#
# Create and run in a container: (BE SURE TO CHANGE /PATH/TO/HIPS/DOCKER e.g., to the directory of this file)
#     docker run --name hips_container -v /path/to/hips/docker:/home --hostname docker -it hips_img
# 
#     This is setup for an interactive shell. 
#     The folder ~/docker will be created in your home (host) if not already there and connected to /home in the docker container
#     This makes it easy to access and share files back and forth
#     add the --rm flag to delete the container on exit, like docker run --rm --name ...
#
# Connect to an existing container: 
#     docker start -ai hips_container
# 
#     -a attaches, and -i is interactive
#
# When entering the docker container you will start in /app where hips is installed.
# You can do "cp -r hips /home" to copy the install to /home/hips.
# This location is mounted to the host at /path/to/hips/docker (in the docker run command above).
# You will see that /path/to/hips/docker/hips is now present on the host 
# and matches /home/hips in the container.
# Inside the container you can run the hips code, 
# and then view/plot/edit the files on your local host.
#
#################################################################################

FROM alpine:3.22.0

# Install the necessary packages
RUN apk add --no-cache wget tar cmake make g++ lapack-dev vim python3-dev catch2-3

############## Build Sundials

WORKDIR /app/

RUN wget https://github.com/LLNL/sundials/releases/download/v7.3.0/sundials-7.3.0.tar.gz \
&&  tar -xvzf sundials-7.3.0.tar.gz \
&&  rm sundials-7.3.0.tar.gz \
&&  cd sundials-7.3.0 \
&&  mkdir build \
&&  cd build \
&&  cmake .. \
&&  make -j8 \
&&  make install \
&&  cd ../../ \
&&  rm -r sundials-7.3.0

############## Build Yaml

RUN apk add --no-cache git

RUN git clone https://github.com/jbeder/yaml-cpp.git -b 0.8.0 --depth=1 \
&&  cd yaml-cpp \
&&  mkdir build \
&&  cd build \
&&  cmake .. \
&&  make -j8 \
&&  make install \
&&  cd ../../ \
&&  rm -r yaml-cpp

############## Build Cantera

RUN apk add --no-cache scons python3 boost-dev hdf5-dev
RUN apk add --no-cache py3-setuptools

RUN git clone https://github.com/Cantera/cantera.git -b 3.1 --depth=1 \
&&  cd cantera \
&&  scons build -j8 python_package=n \
&&  scons install \
&&  cd ../ \
&&  rm -r cantera

############## Build hips
# note,
# git clone https://github.com/BYUignite/hips.git --depth=1 --branch v1.0 \
# will download version tagged v1.0, and similarly for other tags

# RUN git clone https://github.com/BYUignite/hips.git --depth=1 \
# && cd hips \
# &&  mkdir build \
# &&  cd build \
# &&  cmake -D REACTIONS_ENABLED=ON -D HIPS_BUILD_DOCS=ON -D HIPS_BUILD_EXAMPLES=ON -D HIPB_BUILD_TESTS=ON .. \
# &&  make -j \
# &&  make install \
# &&  cd ../

# Set the working directory
WORKDIR /app/hips

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/sh"]


