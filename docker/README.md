# Using HiPS with Docker

## Quick Start
1. Install Docker
2. Navigate to `hips/docker/`
3. Run `./launch_hips.sh`

## How It Works
- Your HiPS source code stays on your host machine
- The container provides the build environment (dependencies, compilers)
- Files are shared between host and container via bind mount
- HiPS builds automatically on first run, skips building on subsequent runs

## Development Workflow
1. Edit code on your host with any editor
2. In container: `cd build && make -j8 && make install` (to recompile)
3. In container: `cd run && ./ex_1.x` (to run examples)
4. Results appear on your host for analysis

## Notes
- First build takes ~10 minutes, subsequent starts are instant
- Container is deleted on exit but your work persists on the host
- Modify `Dockerfile` and rerun script to update the environment