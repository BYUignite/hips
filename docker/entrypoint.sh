#!/bin/sh

if [ ! -d "installed" ]; then

    mkdir build
    cd build
    cmake -D REACTIONS_ENABLED=ON -D HIPS_BUILD_DOCS=ON -D HIPS_BUILD_EXAMPLES=ON -D HIPB_BUILD_TESTS=ON ..
    make -j
    make install
    cd ../

fi

exec "$@"