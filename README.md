[![arXiv][arxiv-shield1]][arxiv-url1]
[![arXiv][arxiv-shield2]][arxiv-url2]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15186632.svg)](https://doi.org/10.5281/zenodo.15186632)

# Continuous-Time Data-Driven Control

This repository contains the code for the activities of the project IMPACT4Mech on data-driven control of continuous-time linear time-invariant systems.

**Main files:** these files contain the MATLAB code for the numerical examples of the paper:
> A. Bosso, M. Borghesi, A. Iannelli, G. Notarstefano, A. R. Teel, "Data-Driven Control of Continuous-Time LTI Systems via Non-Minimal Realizations."

- _batch_reactor_v1.m_: numerical example of data-driven stabilization and control with integral action.

- _surface_vessel_v1.m_: numerical example of data-driven output regulation.

- _canonical.m_: MATLAB function that computes the multivariable controller canonical form of a pair (A, B).

- _find_nu.m_: MATLAB function that estimates the observability index of a plant from a dataset.

- _indices.m_: MATLAB function that computes the controllability indices of a pair (A, B).

- _realization.m_: MATLAB function that computes a canonical non-minimal realization of a plant.

- _stabilize.m_: MATLAB function that computes a stabilizing gain from a dataset.

- _X_solve.m_: MATLAB function that solves the matrix equation XΘ = ΘX, Xβ = φ.

**ecc2025**: folder containing the MATLAB code for the numerical examples of the paper:
> A. Bosso, M. Borghesi, A. Iannelli, G. Notarstefano, A. R. Teel, "Derivative-Free Data-Driven Control of Continuous-Time Linear Time-Invariant Systems." 2025 European Control Conference (ECC).

Content of the folder:

- _ecc2025_algorithm1_v1.m_: numerical example for data-driven control with input-state data

- _ecc2025_algorithm2_v1.m_: numerical example for data-driven control with input-output data



The files require the installation of MOSEK and YALMIP:

MOSEK:  https://docs.mosek.com/10.2/toolbox/index.html

YALMIP: https://yalmip.github.io

## Contact

Alessandro Bosso

alessandro.bosso@unibo.it


[arxiv-shield1]: https://img.shields.io/badge/arxiv-2410.24167-t?style=flat&logo=arxiv&logoColor=white&color=red
[arxiv-shield2]: https://img.shields.io/badge/arxiv-2505.22505-t?style=flat&logo=arxiv&logoColor=white&color=red
[arxiv-url1]: https://arxiv.org/abs/2410.24167
[arxiv-url2]: https://arxiv.org/abs/2505.22505

## Acknowledgments

The research leading to these results has received funding from the European Union's Horizon Europe research and innovation program under the Marie Skłodowska-Curie Grant Agreement No. 101104404 - IMPACT4Mech. https://cordis.europa.eu/project/id/101104404
