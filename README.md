<hr>

<h1><p align="center"> Cell centred explicit solid dynamics toolkit for OpenFOAM
</p></h1>

<p align="center">
  <a href="https://travis-ci.org/jibranhaider/explicitSolidDynamics" target="_blank">
    <img alt="Travis (.org)" src="https://img.shields.io/travis/jibranhaider/explicitSolidDynamics.svg?label=Build"> &nbsp;
  </a>
  <img alt="OpenFOAM 7" src="https://img.shields.io/badge/OpenFOAM-v7-darkgreen.svg"> &nbsp;
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/blob/master/LICENSE">
    <img alt="GPLv3 license" src="https://img.shields.io/badge/License-GPLv3-orange.svg"> &nbsp;
  </a>
  <img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/jibranhaider/explicitSolidDynamics.svg?color=blue?label=Latest%20commit">
</p>

<br/>

## 1. Introduction

<p align="justify">
This novel toolkit is based on a cell centred Finite Volume Method to predict large deformation in solids. The governing equations comprise of first order hyperbolic conservation laws for linear momentum and deformation gradient tensor of the system. This helps to bridge the gap between Computational Fluid Dynamics and Computational Solid Dynamics. The governing equations are spatially discretised using a low order cell centred Finite Volume Method along with an upwind Riemann solver.
</p>

### Key features
* Second order accuracy for velocities and stresses/strains.
* Parallelised C++ implementation within [OpenFOAM](http://openfoam.org/) code.
* Robust methodology with practical applications.

More details about this work can be found [here](https://jibranhaider.com/research/explicit-solid-dynamics-in-openfoam/).

<br/>
<hr>

## 2. Documentation

* [Installation](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Installation)
* [Tutorials](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Tutorials)
* [Upcoming features](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Upcoming-features)
* [Contributors](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Contributors)
* [References](https://github.com/jibranhaider/explicitSolidDynamics/wiki/References)

<br/>
<hr>

## 3. Applications

### Implosion of a bottle
The simulation can be seen [here](https://youtu.be/K6T6OTjlOzQ).

<img src="/docs/results/implodingBottle/deformation/10.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/20.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/25.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/30.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/35.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/40.png" width="14%">

<br/>

### Crushing of a thin cylinder
The simulation can be seen [here](https://youtu.be/7_0i0FlCHtc).

<img src="/docs/results/crushingCylinder/deformation/10.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/15.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/20.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/25.png" width="23%">

<img src="/docs/results/crushingCylinder/deformation/30.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/35.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/40.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/deformation/50.png" width="23%">


<br/>
<hr>

## 4. How to cite this work?
If you are using this toolkit for your research then please cite it as follows.

<br/>

* [![GITHUB](https://zenodo.org/badge/DOI/10.5281/zenodo.2577032.svg)](https://doi.org/10.5281/zenodo.2577032)

      @misc{haider2019toolkit,
      author    = {Haider, Jibran},
      title     = {{ExplicitSolidDynamics toolkit for OpenFOAM}},
      year      = {2019},
      doi       = {10.5281/zenodo.2577033},
      url       = {https://github.com/jibranhaider/explicitSolidDynamics}}

<br/>
<hr>

## 5. License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file.


<br/>
<hr>

## 6. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.