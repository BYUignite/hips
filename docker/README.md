# Using HiPS with Docker

## Quick Start
1. Install [Docker](https://www.docker.com/get-started/) and make sure it is up and running.
2. Clone [HiPS](https://github.com/BYUignite/hips) to desired host folder.
3. Navigate to `path_to_hips/hips/docker/`.
4. Run `./launch_hips.sh` which will:
     
     - Build the docker image if it hasn't been built yet or the Dockerfile has been changed
     - Open the HiPS container in the terminal, ready to run executables
  
5. Inside the container, run `cd run && ./ex_1.x` to verify the Docker setup is working. Data will be kept in `hips/data/`.
6. [Data processing](../README.md#post-processing) and analysis can be done with the scripts in `hips/post/` but should be done on the host. The container doesn't include the dependencies.
7. Run `exit` to close the container.

______________

## How It Works

**Docker Setup:**
- The `launch_hips.sh` script initializes the docker setup and can be used to enter the container at any time
- The first run produces an `entrypoint.sh` file that is run anytime the container is open and ensures HiPS is installed in the container; however, it doesn't rebuild HiPS if source code is edited, so you must rebuild manually after making changes.
- The container provides the complete build environment (dependencies, compilers) for HiPS

**File Sharing & Workflow:**
- All HiPS files including source code are shared between host and container via bind mount
- Data processing and analysis should be done on the host (container lacks post-processing dependencies for python)

**Building & Rebuilding:**
- HiPS builds automatically on first run and skips building on subsequent runs with `launch_hips.sh`
- If the Docker is modified, running `launch_hips.sh` will always update the image with the changes.
- **When source code is modified:** Run `cd build && make -j8 && make install` inside the container to rebuild