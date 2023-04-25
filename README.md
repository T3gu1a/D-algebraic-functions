# D-algebraic functions in Maple

**NLDE** (NonLinear algebra and Differential Equations) is a [Maple](https://www.maplesoft.com/) package to work with D-algebraic functions. These are functions that satisfy algebraic differential equations (ADEs), i.e., differential equations that are polynomial in the independent variable and derivatives of the dependent variables. The package provides:

- **unaryDalg**: for computing ADEs for rational expressions of a single D-algebraic function using elimination with Groebner bases.
- **arithmeticDalg**: for computing ADEs for rational expressions of D-algebraic functions using elimination with Groebner bases.
- **composeDalg**: for computing ADEs for compositions of D-algebraic functions partly using elimination with Groebner bases.
- **diffDalg**: for computing ADEs for derivatives of D-algebraic functions using recursive elimination with by computing resultant.
- **invDalg**: for computing ADEs for functional inverses of D-algebraic functions by explicit construction.
- **AnsatzDalg**: a subpackage with main sub-procedures **unaryDeltak** and **arithmeticDeltak** for doing the same computation (with some extensions) as _unaryDalg_ and _arithmeticDalg_ by an algorithmic search based on linear algebra.
- **DDfiniteToDalg**: for converting a holonomic ODE into an ADE having the same solutions.
- **SysToMinDiffPoly**: for computing input-output equation of dynamical systems.
- **OrderDegreeADE**: for computing the order and the degree of a given ADE. Often useful when the ADE displays on several lines.
- **MultiDalg**: subpackage for operations with multivariate D-algebraic functions. The command for arithmetic operations is almost complete! (April 2023).

## Installation

The easiest way to use **NLDE** in Maple is by putting the file NLDE.mla in your working directory and include the lines
```
  > restart;

  > libname:=currentdir(), libname:

  > with(NLDE) 
```
at the beginning of your Maple worksheet (session). To avoid putting these three lines in all worksheets, one can read the help page of the $\texttt{libname}$ command.

## Requirements and Dependencies

The package can be used with any recent version of Maple (from 2019 onward). 
For Groebner bases computations the package relies on the following Maple packages:
- $\texttt{Groebner}$
- $\texttt{PolynomialIdeal}$

## Author

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), Max Planck Institute for Mathematics in the Sciences (2022-present).
- licence: GNU General Public Licence v3.0.

## Documentation and Examples

Documentation and examples are given at [NLDE documentation](https://T3gu1a.github.io/NLDEdoc/).

## Examples in a Maple worksheet

One can try the examples of the documentation in the Maple worksheet MapleWorksheet-NLDEdoc-examples.mw.

## References

1. [D-algebraic functions](https://arxiv.org/abs/2301.02512). Rida Ait El Manssour, Anna-Laura Sattelberger, Bertrand Teguia Tabuguia. January 2023.

2. [Operations for D-algebraic functions](https://arxiv.org/abs/2304.09675). Bertrand Teguia Tabuguia. April 2023.

3. [Arithmetic of D-algebraic functions](https://) (link to come). ... and Bertrand Teguia Tabuguia. _Work in progress_. We investigate multivariate D-algebraic functions in more details here.

