# GNC: A Monte Carlo code for dynamics of Galactic Nuclear star Cluster with a central massive black hole

GNC is an open-source **Monte Carlo code** designed to simulate the complex dynamical evolution of stars and compact objects orbiting the supermassive black holes (SMBHs) at the centers of galaxies. Developed by **Fupeng Zhang** and **Pau Amaro Seoane**, GNC provides a powerful tool for researchers in astrophysics and stellar dynamics.

-----

## Key Features

  * **Two-Dimensional Fokker-Planck Method:** Tracks the evolution of particles in energy and angular momentum space, which is highly efficient for simulating long-term relaxation.
  * **Multi-Mass Components:** Natively handles clusters containing different types of stars and compact objects (e.g., stellar-mass black holes, neutron stars).
  * **Core Physics Included:** Accurately models **two-body relaxation** and the effects of the **loss cone**, where objects are consumed by the central black hole.
  * **High Accuracy for Rare Particles:** Implements a weighting method to improve the statistical accuracy for rare but dynamically important objects.
  * **Flexible and Extensible:** The code is designed to be easily extended to include more complex physical processes, such as relativistic effects.

-----

## Getting Started

### Prerequisites

The code is written in **Fortran 2003**. To compile and run it, you'll need a Fortran compiler and an MPI (Message Passing Interface) library.

GNC has been successfully tested with the following common combinations:

  * `gfortran` with `Open MPI`
  * `ifort` with `MPICH`

### Installation and Usage

For detailed instructions on how to compile the code, set up a simulation, and run it, please refer to the complete user guide:

```
Doc/tutorial.pdf
```

-----

## Citing GNC

If you have used GNC in your work, we kindly ask you to please cite the original reference paper:

> Zhang, F., & Amaro-Seoane, P. 2024, The Astrophysical Journal, 961, 232


-----

## Scientific Context

The GNC code was used to investigate the long-term evolution of the structure of stellar clusters around SMBHs. A key finding from the initial study is the development of **tangential anisotropy** (a preference for more circular orbits) in the cluster's central region, while the outer regions remain isotropic. This showcases the code's ability to produce detailed and insightful results for complex astrophysical problems, such as predicting event rates for gravitational wave sources. 

-----

## Questions and Support

If you have any questions, find a bug, or have suggestions for improvements, please feel free to contact the authors:

  * **Fupeng Zhang:** `zhangfupeng@gzhu.edu.cn`
  * **Pau Amaro Seoane:** `amaro@upv.es`

  


