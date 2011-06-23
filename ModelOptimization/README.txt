Project: ModelOptimization
README.txt author: Cory Quammen <cquammen@cs.unc.edu>
Last updated: 2011/06/23

*** BUILD INSTRUCTIONS FOR LINUX AND MAC ***

** PREREQUISITES **

- CMake 2.8 (available at http://www.cmake.org/cmake/resources/software.html)
- A compiler

** CONFIGURING **

1). Create a build directory outside the git repository in which this
    file is located. For example, if the directory containing this
    README.txt file is located in

    /home/user/cquammen/madai/stat/ModelOptimization
    
     a good build directory might be
    
    /home/user/cquammen/madai/build/stat/ModelOptimization-build

2). Navigate to the build directory created in step 1.

3). Type

     cmake .../stat/ModelOptimization

     where .../stat/ModelOptimization is an absolute or relative path
     to the ModelOptimization source directory.

     Using the example directories above, this command is

     cmake /home/user/cquammen/madai/stat/ModelOptimization

4). Type

    make

    The various source files will compile, producing a library
    libModelOptimization.a in a directory called "lib" in the root
    level of the build directory.
