#!/bin/bash

# Build the Docker image. If the image exists it will use cached
# instructions.
docker build -t hips_img .

# If a container named 'hips' is already running, stop it. 
if [ -n "$(docker ps -q -f name=hips -f status=running)" ]; then
    docker stop hips
fi

# Run the Docker container with a bind mount to the host's hips directory.
docker run --rm --name hips -it -v "$(pwd)/..":/app/hips hips_img   