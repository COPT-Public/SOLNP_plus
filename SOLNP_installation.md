## Requirements

CMake  >= 3.5.1

[MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) or [LAPACK](https://netlib.org/lapack/explore-html/index.html)  

[OSQP](https://osqp.org)

## Installation

Currently, we provide C, MATLAB, Python and Julia interface on macOS, Linux and Windows. 

In terms of the linear algebra library, we provide support for both Intel MKL and open-source LAPACK, which should be installed in advance. On macOS, there is no need to install any linear algebra library, as we leverage Apple's proprietary **Accelerate framework**. On Linux and Windows, you can decide whether to use MKL by configuring the parameter **LINK_MKL** in **CMakeLists.txt**. If LAPACK or MKL can not be found automatically, you can specify the path in **CMakeLists.txt** manually. Before installation, you also need to specify the environment variables of OSQP directories:

```bash
export OSQP_HOME=your_path_to_osqp
```

We compile MATLAB interface within MATLAB, you can simply run the following command in MATLAB command line:

```bash
cd SOLNP_PLUS_ROOT
cd interface/Matlab
make_solnp
```

As for C, Python and Julia interface, you need to first compile SOLNP+ source code:

On Linux and macOS :

```bash
cd SOLNP_PLUS_ROOT
mkdir build
cd build
cmake ..
cmake --build .
```

On Windows:

```bash
cd SOLNP_PLUS_ROOT
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
cmake --build .
```

After compiling source code, you can directly use the C interface. We give **CMakeLists.txt** and a simple example, from which you can compile executable file:

```bash
cd SOLNP_PLUS_ROOT
cd interface/C/build
cmake ..
cmake --build .
./solnp_demo
```

You can also Install and test Python interface:

```bash
cd SOLNP_PLUS_ROOT
cd interface/Python
python setup.py install
cd test
python test.py
```

And Julia interface:
```bash
cd SOLNP_PLUS_ROOT
cd interface/Julia
```

- Type in `julia` in command line to start a Julia REPL. Then type in `]` to enter Pkg mode and type in `dev ./SOLNP/`. Type in `[Backspace]` to leave Pkg mode.
- To test whether everything is ready, try to run test file in Julia REPL by typing in
```julia
julia> include("./SOLNP/test/test.jl")
```
