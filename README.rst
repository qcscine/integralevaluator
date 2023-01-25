SCINE - Integral Evaluator
=============================

Introduction
------------

This repository contains integral evaluation routines, relying on 
the Libint library by Ed Valeev. 
The functionality includes:

- Born--Oppenheimer (BO) and pre-BO one- and two-body integrals.
- Storing of one-body integrals in memory.
- Storing of two-body integrals in memory.
- Direct Fock matrix evaluation routines for BO and pre-BO.
- Cauchy--Schwartz screening for the BO contribution to the Fock matrix.
- First derivatives of one- and two-body integrals.

License and Copyright Information
---------------------------------

For license and copyright information, see the file ``LICENSE.txt`` in this
directory.

Installation and Usage
----------------------

As a user of one of the SCINE modules (such as Kiwi), you do not need
to set up the Integral Evaluator yourself. Instead, this is done as part of the
installation process of the respective SCINE module. Therefore, the following
instructions are only necessary for developers of SCINE, or those interfacing
this library directly.

Dependencies
............

Required software, minimum required versions in brackets, for this SCINE project are:

- A C++ compiler supporting the C++17 standard (at the moment only GCC is supported)
- Libint 2.7.2 (other versions are not tested) [Will be automatically downloaded if not found]
- CMake (3.9)
- Boost (1.65.0)
- Eigen3 (3.3.2)

Installation
............

In order to compile this as a SCINE developer, execute the following
commands::

    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst ..
    make -j 4
    make test
    make install

Known Issues
------------

Intel compilers do not work due to libint at the moment.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE IntegralEvaluator makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
- `libint <https://github.com/evaleev/libint>`_
