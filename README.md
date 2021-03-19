# GYT: Generalized Young Tableaux

## Table of content

* [Brief description](#brief-description)
* [Compilation](#compilation)
* [Installation](#installation)
* [Invocation](#invocation)

## Brief description

This is a set of C++ programs generalizing Young tableaux to higher
dimensions, designed for fast computation of of solutions for variadic
polynomial equations `p(x_1, ..., x_n) = B` over non-negative
integers. Read the paper `gyt.pdf` for more information.

This distribution consists of the following C++ sources:

    gyt-2d-all.cpp
    gyt-2d-all-gmp.cpp
    gyt-2d.cpp
    gyt-2d-gmp.cpp
    gyt-all.cpp
    gyt-all-gmp.cpp
    gyt-common.cpp
    gyt-common-gmp.cpp
    gyt-common-gmp.hpp
    gyt-common.hpp
    gyt.cpp
    gyt-gmp.cpp
    gyt-pq-common-gmp.hpp
    gyt-pq-common.hpp
    gyt-pq.cpp
    gyt-pq-gmp.cpp
    gyt-proba.cpp
    gyt-proba-gmp.cpp
    gyt-rand.cpp
    gyt-rand-gmp.cpp

The sources `gyt-2d-*.cpp` implemet the original two-dimensional Young
tableaux. The sources `gyt-common*.cpp` and `gyt-common*.hpp` contain
the parts common to all variants of GYT. The headers
`gyt-pq-common*.hpp` describe the priority structure for `gyt-pq*.cpp`
sources. All other sources are described in detail in the paper
`gyt.pdf`.

## Compilation

For the full compilation of these programs, the C++ version of the
[GNU Multiple Precision Arithmetic Library](https://gmplib.org/) (GMP)
must be installed.

To compile all programs, write the command
```Makefile
    make
```
in the root directory. If you do not have GMP installed, write the
command
```Makefile
    make simple
```
in the root directory.

## Installation

The compilation process installs the binaries in the root directory of
the distribution. To install the programs on your computer, write the
command
```Makefile
    make install
```
You need to have superuser privileges to perform this command.

## Invocation

All binaries accept the input from `STDIN` and write the output on
`STDOUT`. See the examples in the subdirectory `data`.
