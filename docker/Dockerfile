
#################################################################################
# Dockerfile for building a Docker image with Sundials, Yaml, Cantera, and hips
#
# See README.md for usage instructions.
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

############## Build hips with entrypoint script
WORKDIR /app/hips

RUN cat > /entrypoint.sh << 'EOF' && chmod +x /entrypoint.sh
#!/bin/sh

if [ ! -d "installed" ]; then

    mkdir -p build
    cd build
    cmake -D REACTIONS_ENABLED=ON -D HIPS_BUILD_DOCS=ON -D HIPS_BUILD_EXAMPLES=ON -D HIPB_BUILD_TESTS=ON ..
    make -j
    make install
    cd ../

fi

exec "$@"
EOF

# Set the working directory
WORKDIR /app/hips

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/sh"]


