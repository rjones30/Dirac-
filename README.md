# Dirac++ - general root/c++ toolkit for computing Feynman amplitudes

## Author

* Richard Jones, University of Connecticut, Storrs, CT

## Description

A general-purpose toolkit for use within the cern/root framework for
computing cross sections and rates with all polarization observables
under the control of the user for both incoming and outgoing particles.
It is decomposed into three lower-level packages that provide the
functionality of 4-vectors (see LorentzPackage.h), Pauli algebra
(see PauliPackage.h), and Dirac algebra (see DiracPackage.h), with
the top-level TCrossSection class integrating them all together.

## History

This package was originally written to provide the compute the polarization
of Compton scattered electrons and photons at GeV-scale energies. From this
beginning, it was fairly simple to generalize into a tool for any QED
process for which the Feynman amplitude can be written down as a sum of
tree-level graphs. Extension of this package to support computation of
higher-order graphs (radiative corrections) is not available at this time.

## Release history

See VERSIONS file in the project directory.

## Usage synopsis

General functions that support plotting of cross sections and particle
generation are provided for specific QED processes in the source files
with names that end with ".C". These functions call methods of the
TCrossSection class to compute differential cross sections and rates.
See the comments within those source files for details.

## Dependencies

You need to have installed cern/root version 6 or greater 
(see http://root.cern.ch) and a working c++ compiler with
support for c++11 language features.

## Building instructions

The core code is compiled into a shared library libDirac.so by the
command:

    $ make

After this, it can be imported into a root session and used by client
code, as illustrated here.

    $ root -l
    root [0] .L libDirac.so
    root [1] .L Brems.C+
    root [3] .L tests.C+
    root [4] TestThreeVectorReal()
    Zero(): TThreeVectorReal(0,0,0)
    Set to 1,2,3: TThreeVectorReal(1,2,3)
    Initialized to 4,5,6: TThreeVectorReal(4,5,6)
    Initialized to (Float_t *) 7,8,9: TThreeVectorReal(7,8,9)
    Initialized to (Double_t *) 10,11,12: TThreeVectorReal(10,11,12)
    Initialized to (TThreeVectorReal) 1,2,3: TThreeVectorReal(1,2,3)
    SpaceInv(): TThreeVectorReal(-1,-2,-3)
    Normalize(1): TThreeVectorReal(-0.267261,-0.534522,-0.801784)
    SetPolar(10,2.5,-1.2) : GetPolar(10,2.5,-1.2)
    z == z ? yes!
    z != w ? yes!
    Rotate(phi,theta,0) gets back z axis? yes!
    (w cross x) != 0 ? yes!
    w dot (w cross x) == 0 ? yes!
    w cross (w - x) == x cross w ? yes!
    root [5]

## Documentation

See comments at the head of specific process implementation sources:

1. Brems.C - bremsstrahlung process
2. Compton.C - Compton scattering process
3. Pairs.C - coherent pair production process
4. Triplets.C - incoherent pair production process

## Troubleshooting

The code is open-source, and the behavior should be apparent from
the names of the methods. For help with bugs and questions about
how to implement additional processes, please contact the author.

## Bugs

None known at this time.

## How to contribute

## Contact the authors

Write to richard.t.jones at uconn.edu.
