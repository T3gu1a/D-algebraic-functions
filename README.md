# D-algebraic functions in Maple

**Goal**: The ultimate goal of this project is to develop and implement efficient algorithms for operations involving **D-algebraic functions and sequences**. While elimination methods based on Gröbner bases are often too general for specific problems, our aim is to find highly optimized approaches for D-algebraicity-related elimination challenges that avoid generic Gröbner basis computations.

A prime example is the ``CCfiniteToSimpleRatrec`` command, part of the **DalgSeq** subpackage within **NLDE**. This algorithm makes no explicit use of Gröbner bases, or even explicit linear algebra, yet its output perfectly matches that of the Gröbner bases approach. It consistently returns a rational recursion for the $C^2$-finite sequence under consideration.

**NLDE** (NonLinear algebra and Differential/Difference Equations) is a [Maple](https://www.maplesoft.com/) package to work with D-algebraic functions. These are functions that satisfy algebraic differential equations (ADEs), i.e., differential equations that are polynomial in the independent variable and derivatives of the dependent variables. Some features for the difference case are bieng implemented. The package provides:

- ``unaryDalg``: for computing ADEs for rational expressions of a single D-algebraic function using elimination with Groebner bases.
- ``arithmeticDalg``: for computing ADEs for rational expressions of D-algebraic functions using elimination with Groebner bases.
- ``composeDalg``: for computing ADEs for compositions of D-algebraic functions partly using elimination with Groebner bases.
- ``diffDalg``: for computing ADEs for derivatives of D-algebraic functions using recursive elimination with by computing resultant.
- ``invDalg``: for computing ADEs for functional inverses of D-algebraic functions by explicit construction.
- ``AnsatzDalg``: a subpackage with main sub-procedures ``unaryDeltak`` and ``arithmeticDeltak`` for doing the same computation (with some extensions) as _unaryDalg_ and _arithmeticDalg_ by an algorithmic search based on linear algebra.
- ``DDfiniteToDalg``: for converting a DD-finite ODE into an ADE whose set of solutions contains those of the DD-finite ODE.
- ``SysToMinDiffPoly``: for computing input-output equation of dynamical systems.
- ``OrderDegreeADE``: for computing the order and the degree of a given ADE. Often useful when the ADE displays on several lines.
- ``MultiDalg``: subpackage for operations with multivariate D-algebraic functions. The command, ``arithmeticMDalg``, for arithmetic operations is now available! (May 2023).
- **DalgSeq**: subpackage for operations with difference-algebraic sequences. Its main commands are given below.
  - ``HoloToSimpleRatrec``: to convert a given holonomic equation into a rational recursion satisfied by *most* solutions of that holonomic equation.
  - ``CCfiniteToSimpleRatrec``: to convert a given $C^2$-finite equation into a rational recursion satisfied by *most* solutions of that $C^2$-finite equation.
  - ``DalgGuess``: to search for an algebraic difference equation from finitely many first terms of an unknown sequence.
  - ``arithmeticDalgSeq``, ``unaryDalgSeq``, ``AnsatzDalgSeq``, ``OrderDegreeRec``, ``CCfiniteToDalg``, ``PartialSumDalgSeq``, ``PartialProdDalgSeq``. These commands are either self explanatory or are defined in a similar manner as their differential conterparts.



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

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), Max Planck Institute for Mathematics in the Sciences (2022-2023). University of Oxford (2023 - present)
- licence: GNU General Public Licence v3.0.

## Documentation and Examples

An old documentation with examples are given at [NLDE documentation](https://T3gu1a.github.io/NLDEdoc/).

## Examples in a Maple worksheet

One can try the examples of the documentation in the Maple worksheet MapleWorksheet-NLDEdoc-examples.mw. The file DAlgebraicFunc-Examples-Maple-Worksheet.pdf contains examples from the paper D-Algebraic Functions. The corresponding worksheet is also provided.

## References

1. [D-algebraic functions](https://arxiv.org/abs/2301.02512). Rida Ait El Manssour, Anna-Laura Sattelberger, Bertrand Teguia Tabuguia. January 2023. Journal of Symbolic Computation. DOI: [https://doi.org/10.1016/j.jsc.2024.102377](https://doi.org/10.1016/j.jsc.2024.102377).

2. [Operations for D-algebraic functions](https://arxiv.org/abs/2304.09675). Bertrand Teguia Tabuguia. April 2023. ACM Communications in Computer Algebra, Volume 57, Issue 2. Pages 51--56. June 2023

3. [Arithmetic of D-algebraic functions](https://arxiv.org/abs/2305.00702). Bertrand Teguia Tabuguia. May 2023. Journal of Symbolic Computation. DOI: [https://doi.org/10.1016/j.jsc.2024.102348](https://doi.org/10.1016/j.jsc.2024.102348).
4. [On Rational Recursion for Holonomic Sequences](https://arxiv.org/abs/2404.19136) Bertrand Teguia Tabuguia and James Worrell. April 2024. In: Boulier, F., Mou, C., Sadykov, T.M., Vorozhtsov, E.V. (eds) Computer Algebra in Scientific Computing. CASC 2024. LNCS, vol 14938. Springer, Cham. DOI: [https://doi.org/10.1007/978-3-031-69070-9_18](https://doi.org/10.1007/978-3-031-69070-9_18).

