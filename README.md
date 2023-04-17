i# HiPS (Hierarchical Parcel Swapping)

This code implements the Hierarchal Parcel Swapping (HiPS) model for turbulent reacting or nonreacting flows. 

## Documentation
 * [A. Kerstein, Hierarchical Parcel-Swapping Representation of Turbulent Mixing. Part 1. Formulation and Scaling Properties, Journal of Statistical Physics, 153:142-161 (2013)](https://link.springer.com/content/pdf/10.1007/s10955-013-0811-z.pdf)
 * [A. Kerstein, Hierarchical parcel-swapping representation of turbulent mixing. Part 2. Application to channel flow, Journal of Fluid Mechanics, 750:421-463 (2014)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/hierarchical-parcelswapping-%20%20representation-of-turbulent-mixing-part-2-application-to-channel-flow/19D6D1CAC4D2FAFFC67A67925D7E527B)
 * [A. Kerstein, Hierarchical parcel-swapping representation of turbulent mixing. III. Origins of correlation patterns observed in turbulent boundary layers, Physical Review Fluids, 6, 044611    (2021)](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.6.044611)
## Directory structure
* `cmake_build`: build the code using cmake.
* `data`: contains all data files output during a simulation.
    * The code will generate a subfolder with a name corresponding to case name specified in the run script in the run folder.
        * This case subfolder will contain subfolders `input`, `runtime`, `data`, and `post`, which contain the input data files, runtime output, simulation data files, and post-processed data, respectively.
* `doc`: contains documentation files built using Doxygen. Also includes a simplified, python version of HiPS.
* `input`: contains case input files, notably `input.yaml`, which is the primary input file with code parameters.
    * Other input files are a Cantera mechanism file in `gas_mechanisms` and an optional `restart.yaml` file.
* `post`: contains post-processing scripts and files for given cases. 
   * Output is placed in `data/caseName/post`. These are mostly Python files. Some cases also include experimental data files for comparison and plotting.
* `run`: contains the code executable `hips-run` and several run scripts. These scripts are run by the user to execute the code.
    * The user specifies inputDir as the path to the input file containing the case to run and specifies a case name for variable caseName. Files are created and copied into `data/caseName`, as noted above.
    * The user chooses one of the following run scripts to execute the code: 
      * `runOneRlz.sh` will run a single realization of the code. This is appropriate for some cases, like a statistically stationary channel flow. Many cases require running many realizations to gather turbulent statistics.
      * `runManRlz.sh` will run many realizations on a single processor.
      * `slrmJob.sh` is an example script for running embarrasingly parallel simulations, i.e. one realization for each MPI process.
      * `slrmJob_array.sh` is an example script that runs multiple realizations on a parallel machine using a slurm array.
    * `changeInputParam.py` is a convenience script for changing a value of a variable in the input file. This is convenient when running several cases changing an input parameter and can be used within the run scripts listed above.
* `source`: contains source code (including header files) and `CMakeLists.txt` files.

## Dependencies
### HiPS code
* [Cantera](http://cantera.org): open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport.
* Yaml: input file format. This installation is conveniently built into the HiPS build process. 
* Cmake 3.12 or higher
* (OPTIONAL) Doxygen: builds documentation. 
### Post-processing
Post-processing data produced by HiPS is processed via Python 3 scripts. We recommend Python 3.2 or higher. Scripts may not function properly using Python 2.x. The following packages are required and can be installed via pip3:
* numpy
* scipy
* matplotlib
* glob
* yaml
* sys
* os

## Build
The code is built using CMake. See the README file in the `cmake_build` folder for details.

## Test cases
### HiPS Test
  1. Build the code using CMake.
  2. Navigate to the `run` directory. Open `runOneRlz.sh` and confirm that the input file path and case name are set properly. The defaults in `runOneRlz.sh` are `inputDir="../input/hips"` and `caseName="hips_test"`.
  3. Run `./runOneRlz.sh` to run the case. It should take less than two minutes on an average system. 
  4. Navigate to `post/hips`. 
  5. Run `python3 hips_stats.py [caseName]`. With the default case name, this becomes `python3 hips_stats.py hips_test`. This will generate two plots in `../data/[caseName]/post`. Navigate there to view them. 
  6. In `../data/[caseName]/post`, there should be two newly-generated PDFs. These two plots compare the mean and RMS velocity profiles of the channel flow case just run with HiPS to previous DNS data of the same case. 
