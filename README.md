## CONTACT

CONTACT is an advanced simulation program for the study of three-dimensional
frictional contact problems, such as occur between wheel and rail, in roller
bearings and in offset printing devices. It implements the theories for
rolling contact by Dr.ir. E.A.H. Vollebregt and Prof.dr.ir. J.J. Kalker of
Delft University of Technology. 

## Physics model

CONTACT is intended for concentrated contact problems:
* from given (wheel/rail) profiles, determine where contact occurs;
* identifying the size and shape of contact areas and the pressures acting;
* accounting for creepage and friction, determine the tangential surface
   stresses (tractions), and regions with adhesion and (micro-)slip;
* calculation of elastic displacements and subsurface stresses in the bodies'
   interiors;
* for stationary and instationary, rolling and sliding problems;
* with Hertzian and non-Hertzian, "smooth edged" geometries, such that the
   contact is "concentrated" in a small, almost flat part of the bodies' 
   surfaces;
* for massive homogeneous bodies of elastic and viscoelastic materials.

CONTACT aims to provide complete and detailed solutions, with a powerful
interface for integration in multi-body simulation software. The quality of
Finite Elements, at one thousandth  of the computation time needed.

## Requirements

Compilation currently requires use of the Intel Fortran compiler, mainly
because of the use of the Math Kernel Library (Intel MKL) for FFTs and
LAPACK functionality. The compiler can be obtained for free at intel.com.

## Getting started

Introductory information is provided in the [getting started guide](./doc/getting-started.pdf).

## Reporting issues

We welcome your feedback. You can create an issue, start a discussion, or just
contact Edwin Vollebregt via e-mail or his user profile.

## Contribution

We welcome contributions. Please contact Edwin Vollebregt and we'll discuss
how to work together.

The open source version of CONTACT has been supported financially by the US
Federal Railroad Administration. This is greatly acknowledged.

## Binary version

CONTACT is provided by Edwin Vollebregt, through his company Vtech CMCC.
* An "open source version" is provided for the full CONTACT program and
   CONTACT library version.
* A "convenience version" is provided at a service fee: with additional
   quality checking, documentation, with an automated installer, and
   detailed support.
See the section "Software products" at the website of [Vtech CMCC](https://www.cmcc.nl).

## License

Apache License 2.0, see [LICENSE.md](./LICENSE.md).
