# HiPS (Hierarchical Parcel Swapping)

This code implements the Hierarchal Parcel Swapping (HiPS) model for turbulent reacting or nonreacting flows. 

## Documentation is available at [ignite.byu.edu/hips_documentation](https://ignite.byu.edu/hips_documentation)

* Publications

  * Behrang et al. (2025): *A C++ library for turbulent mixing simulation using hierarchical parcel swapping (HiPS)*, SoftwareX, accepted for publication.
  * Behrang et al. (2025): *Turbulent Mixing of Scalars with Nonunity Schmidt Numbers Using Hierarchical Parcel-Swapping*, *Journal of Fluid Mechanics*, accepted for publication. [Paper](https://ignite.byu.edu/publications/).
  * [Kerstein (2021)](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.6.044611): *Hierarchical Parcel-Swapping Representation of Turbulent Mixing. Part 3. Origins of Correlation Patterns Observed in Turbulent Boundary Layers*, *Physical Review Fluids*, 6:044611.
  * [Kerstein (2014)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/hierarchical-parcelswapping-%20%20representation-of-turbulent-mixing-part-2-application-to-channel-flow/19D6D1CAC4D2FAFFC67A67925D7E527B): *Hierarchical Parcel-Swapping Representation of Turbulent Mixing. Part 2. Application to Channel Flow*, *Journal of Fluid Mechanics*, 750:421–463.
  * [Kerstein (2013)](https://link.springer.com/content/pdf/10.1007/s10955-013-0811-z): *Hierarchical Parcel-Swapping Representation of Turbulent Mixing. Part 1. Formulation and Scaling Properties*, *Journal of Statistical Physics*, 153:142–161.


## Directory Structure

* `build/`: This directory is not included by default and must be created by the user.
    * To compile the code, navigate to this directory and run:
      ```bash
      cmake ..
      make
      make install
      make documentation  # Optional: only if documentation is desired
      ```

* `src/`: Contains all source code for the HiPS implementation.

* `example/`: Contains three default example cases that users can run to get started.

* `run/`: Contains the compiled executables generated during the build process.

* `post/`: Contains Python scripts for post-processing.
    * When users run an example case, a file named `parameters.dat` is automatically saved in the `post/` folder. This file records the parameters used for the run, which may be required by the post-processing scripts.
    * Post-processing results are saved into a subfolder named `processed_data/`, created automatically inside the corresponding `post/` directory.

* `data/`: Created only if data writing is enabled during execution.
    * This happens when the function:
      ```cpp
      void calculateSolution(const double tRun, bool shouldWriteData = true);
      ```
      is called with `shouldWriteData` set to `true`.
    * A subfolder is generated with a name corresponding to the case name specified in the run script (e.g., `rlz_00001` by default).
    * The number of realization folders corresponds to the number of realizations specified by the user.
* `docker/`: Contains the Dockerfile along with a bash script to easily use HiPS with Docker.  

* `docs/`: Contains documentation built with Doxygen.

* `test/`: Contains unit/integration tests using the Catch2 library. There is a git hook in `.githooks/pre-push` that runs the tests before the code can be pushed.
    * configure the hook location by running `git config core.hooksPath .githooks`


## Dependencies

### HiPS code
* [CMake ≥ 3.15](https://cmake.org): build configuration and code compilation.
* C++11-compatible compiler: required to compile the code.
* (OPTIONAL) [Cantera](http://cantera.org): open-source suite for chemical kinetics, thermodynamics, and transport.
  * Only needed if reactions are enabled via `-DREACTIONS_ENABLED=ON` in the CMake configuration.
* (OPTIONAL) [SUNDIALS](https://computing.llnl.gov/projects/sundials): required only for the CVODE integrator used in reaction-enabled builds.
* (OPTIONAL) [Doxygen](https://www.doxygen.nl/): used to build code documentation from annotated source files.
* (OPTIONAL) [Catch2](https://github.com/catchorg/Catch2): library for unit tests. 
    * Only available if the `HIPS_BUILD_TESTS` flag is on, in which case the Catch2 library should be available of the system.

### Post-processing
Post-processing of simulation data is performed using Python 3 scripts. We recommend Python 3.6 or higher. The following packages are required and can be installed via `pip`:

* numpy
* scipy
* matplotlib
* glob *(built-in)*
* os *(built-in)*
* sys *(built-in)*

