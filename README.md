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
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/wiki">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-Wiki-red.svg?label=Documentation"> &nbsp;
  </a>
  <a href="https://doi.org/10.5281/zenodo.2577032">
    <img alt="DOI" src="https://zenodo.org/badge/DOI/10.5281/zenodo.2577032.svg"> &nbsp;
  </a>
  <a href="https://jibranhaider.com/research/explicit-solid-dynamics-in-openfoam">
    <img alt="Website" src="https://img.shields.io/website?down_message=red&label=Website&up_color=green&url=https%3A%2F%2Fimg.shields.io%2Fwebsite%3Furl%3Dhttps%253A%252F%252Fjibranhaider.com%252Fresearch%252Fexplicit-solid-dynamics-in-openfoam">
  </a>   
</p>


<p align="center">   
  <img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/jibranhaider/explicitSolidDynamics.svg?label=Latest%20commit"> &nbsp;
  <a href="https://www.codacy.com/app/jibranhaider/explicitSolidDynamics?utm_source=github.com?&amp;utm_medium=referral&amp;utm_content=jibranhaider/explicitSolidDynamics&amp;utm_campaign=Badge_Grad">
    <img alt="Codacy grade" src="https://img.shields.io/codacy/grade/8d44ec596cb64da28c1713e7f88051be.svg?color=darkgreen&label=Code%20quality"> &nbsp; </a>
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/issues">
    <img alt="GitHub open issues" src="https://img.shields.io/github/issues-raw/jibranhaider/explicitSolidDynamics.svg?label=Open%20issues&color=green"> &nbsp;
  </a>
  <img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/jibranhaider/explicitSolidDynamics.svg?color=lightgrey"> &nbsp;
  <img alt="GitHub repo size in bytes" src="https://img.shields.io/github/repo-size/jibranhaider/explicitSolidDynamics.svg?label=Repository%20size&color=lightgrey"> &nbsp; 
  <a href="https://openfoam.org/dev/coding-style-guide">
  <img alt="Code style" src="https://img.shields.io/badge/Coding_style-OpenFOAM-yellow.svg"> 
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

Detailed instructions on the usage of the toolkit are provided on the [wiki](https://github.com/jibranhaider/explicitSolidDynamics/wiki) page.


<br/>
<hr> 

## 3. Applications

#### Implosion of a bottle
The simulation can be seen [here](https://youtu.be/K6T6OTjlOzQ).

<img src="/docs/results/implodingBottle/deformation/10.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/20.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/25.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/30.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/35.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/deformation/40.png" width="14%">

<br/>

#### Crushing of a thin cylinder
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
If you are using this toolkit for your research then please cite the below mentioned work.

<br/>

* [![GITHUB](https://zenodo.org/badge/DOI/10.5281/zenodo.2577032.svg)](https://doi.org/10.5281/zenodo.2577032)

      @misc{haider2019toolkit,
      author    = {Haider, Jibran},
      title     = {{ExplicitSolidDynamics toolkit for OpenFOAM}},
      year      = {2019},
      doi       = {10.5281/zenodo.2577033},
      url       = {https://github.com/jibranhaider/explicitSolidDynamics}}

<br/>

* [![CMAME](https://img.shields.io/badge/DOI-10.1016/j.cma.2018.06.010-blue.svg)](https://doi.org/10.1016/j.cma.2018.06.010)

      @article{haider2018,
      title     = {{An upwind cell centred Total Lagrangian finite volume algorithm for nearly incompressible explicit fast solid dynamic applications}},
      author    = {Haider, Jibran and Lee, Chun Hean and Gil, Antonio J and Huerta, Antonio and Bonet, Javier},
      journal   = {Computer Methods in Applied Mechanics and Engineering},
      volume    = {340},
      pages     = {684--727},
      year      = {2018},
      doi       = {10.1016/j.cma.2018.06.010}}

<br/>

* [![IJNME](https://img.shields.io/badge/DOI-10.1002/nme.5293-blue.svg)](https://doi.org/10.1002/nme.5293)

      @article{haider2017,
      title     = {{A first-order hyperbolic framework for large strain computational solid dynamics: An upwind cell centred Total Lagrangian scheme}},
      author    = {Haider, Jibran and Lee, Chun Hean and Gil, Antonio J and Bonet, Javier},
      journal   = {International Journal for Numerical Methods in Engineering},
      volume    = {109},
      number    = {3},
      pages     = {407--456},
      year      = {2017},
      doi       = {10.1002/nme.5293}}


<br/>
<hr> 

## 5. Authors and contributors
This toolkit is developed and maintained by [Jibran Haider](http://jibranhaider.weebly.com/). The following individuals are acknowledged for their support:
* [Dr. Chun Hean Lee](https://www.gla.ac.uk/schools/engineering/staff/chunheanlee/)
* [Prof. Antonio J. Gil](https://www.swansea.ac.uk/staff/engineering/a.j.gil/)
* [Prof. Javier Bonet](https://www.researchgate.net/profile/Javier_Bonet)
* [Prof. Antonio Huerta](https://www.lacan.upc.edu/huerta/)


<br/>
<hr> 

## 6. License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file. 


<br/>
<hr> 

## 7. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
