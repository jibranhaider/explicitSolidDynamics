<hr> 

<h1><p align="center"> Explicit solid dynamics in OpenFOAM 
</p></h1>

<p align="center"> 
  <a href="https://travis-ci.org/jibranhaider/explicitSolidDynamics" target="_blank">
    <img alt="Travis (.org)" src="https://img.shields.io/travis/jibranhaider/explicitSolidDynamics.svg"> &nbsp;
  </a>  
  <img alt="OpenFOAM 6" src="https://img.shields.io/badge/OpenFOAM-6-green.svg"> &nbsp;
  <img alt="OpenFOAM 5" src="https://img.shields.io/badge/OpenFOAM-5-green.svg"> &nbsp;
  <img alt="OpenFOAM 4" src="https://img.shields.io/badge/OpenFOAM-4-green.svg"> &nbsp;
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/blob/master/LICENSE">
    <img alt="GPLv3 license" src="https://img.shields.io/badge/License-GPLv3-coral.svg"> &nbsp;
  </a>
  <a href="https://jibranhaider.weebly.com/research.html">
    <img alt="Website" src="https://img.shields.io/website-up-down-green-red/https/jibranhaider.weebly.com/research.html.svg"> &nbsp;
  </a>  
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/wiki/Manual-for-explicitSolidDynamics-toolkit">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-Wiki-red.svg">
  </a>  
</p>


<p align="center">   
  <img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/jibranhaider/explicitSolidDynamics.svg"> &nbsp;
  <a href="https://www.codacy.com/app/jibranhaider/explicitSolidDynamics?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jibranhaider/explicitSolidDynamics&amp;utm_campaign=Badge_Grad">
    <img alt="Codacy grade" src="https://img.shields.io/codacy/grade/8d44ec596cb64da28c1713e7f88051be.svg?color=darkgreen"> &nbsp; </a>
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/issues">
    <img alt="GitHub open issues" src="https://img.shields.io/github/issues-raw/jibranhaider/explicitSolidDynamics.svg"> &nbsp;
  </a>
  <img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/jibranhaider/explicitSolidDynamics.svg"> &nbsp;
  <img alt="GitHub repo size in bytes" src="https://img.shields.io/github/repo-size/jibranhaider/explicitSolidDynamics.svg?label=Repository%20size"> &nbsp; 
  <a href="https://openfoam.org/dev/coding-style-guide">
  <img alt="Code style" src="https://img.shields.io/badge/Code_style-OpenFOAM-yellow.svg"> 
  </a>  
</p>


<br/>

## 1. Introduction    

<p align="justify">
This novel toolkit is based on a cell centred Finite Volume Method to predict large deformation in solids. The governing equations comprise of first order hyperbolic conservation laws for linear momentum and deformation gradient tensor of the system. This helps to bridge the gap between Computational Fluid Dynamics and Computational Solid Dynamics. The governing equations are spatially discretised using a low order cell centred Finite Volume Method along with an upwind Riemann solver. 
</p> 

#### Key features
* Second order accuracy for velocities and stresses/strains. 
* Parallelised C++ implementation within [OpenFOAM](http://openfoam.org/) code. 
* Robust methodology with practical applications. 

More details about this work can be found [here](http://jibranhaider.weebly.com/research.html).


<br/>
<hr> 

## 2. How to use this toolkit?

Detailed instructions on the usage of the toolkit are provided on this [wiki](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Manual-for-explicitSolidDynamics-toolkit) page.


<br/>
<hr> 

## 3. Applications

#### Implosion of a bottle
The simulation can be seen [here](https://youtu.be/K6T6OTjlOzQ).

<img src="/docs/results/implodingBottle/10.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/20.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/25.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/30.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/35.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/40.png" width="14%">

<br/>

#### Crushing of a thin cylinder
The simulation can be seen [here](https://youtu.be/7_0i0FlCHtc).

<img src="/docs/results/crushingCylinder/10.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/15.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/20.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/25.png" width="23%">

<img src="/docs/results/crushingCylinder/30.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/35.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/40.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/50.png" width="23%">


<br/>
<hr> 

## 4. How to cite this work?

<br/>

* [![CMAME](https://img.shields.io/badge/DOI-10.1016/j.cma.2018.06.010-blue.svg)](https://doi.org/10.1016/j.cma.2018.06.010)

      @article{haider2018upwind,
      title={{An upwind cell centred Total Lagrangian finite volume algorithm for nearly incompressible explicit fast solid dynamic applications}},
      author={Haider, Jibran and Lee, Chun Hean and Gil, Antonio J and Huerta, Antonio and Bonet, Javier},
      journal={Computer Methods in Applied Mechanics and Engineering},
      volume={340},
      pages={684--727},
      year={2018}}

<br/>

* [![IJNME](https://img.shields.io/badge/DOI-10.1002/nme.5293-blue.svg)](https://doi.org/10.1002/nme.5293)

      @article{haider2017first,
      title={{A first-order hyperbolic framework for large strain computational solid dynamics: An upwind cell centred Total Lagrangian scheme}},
      author={Haider, Jibran and Lee, Chun Hean and Gil, Antonio J and Bonet, Javier},
      journal={International Journal for Numerical Methods in Engineering},
      volume={109},
      number={3},
      pages={407--456},
      year={2017}}


<br/>
<hr> 

## 5. Authors and contributors
This toolkit is developed and maintained by [Jibran Haider](http://jibranhaider.weebly.com/) (Swansea University). 

This work would not have been possible without the support and guidance of the following individuals:
* [Dr. Chun Hean Lee](https://www.gla.ac.uk/schools/engineering/staff/chunheanlee/)  (University of Glasgow) 
* [Prof. Antonio J. Gil](https://www.swansea.ac.uk/staff/engineering/a.j.gil/)  (Swansea University)
* [Prof. Javier Bonet](https://www.researchgate.net/profile/Javier_Bonet)  (Greenwich University)
* [Prof. Antonio Huerta](https://www.lacan.upc.edu/huerta/)  (UPC BarcelonaTech)


<br/>
<hr> 

## 6. License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file. 


<br/>
<hr> 

## 7. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.

