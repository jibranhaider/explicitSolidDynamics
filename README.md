<hr> 

# Explicit solid dynamics in OpenFOAM
This novel toolkit is based on a cell centred Finite Volume Method to predict large deformation in solids. The governing equations comprise of first order hyperbolic conservation laws for linear momentum and deformation gradient tensor of the system. This helps to bridge the gap between Computational Fluid Dynamics and Computational Solid Dynamics. The governing equations are spatially discretised using a low order cell centred Finite Volume Method along with an upwind Riemann solver. 

### Key features
* Second order accuracy for velocities and stresses/strains. 
* Parallelised C++ implementation within [OpenFOAM](http://openfoam.org/) code. 
* Robust methodology with practical applications. 

More details about this work can be found [here](http://jibranhaider.weebly.com/research.html).



<br/><br/>
<hr> 

# How to use this toolkit?

Detailed instructions on the usage of the toolkit are provided on this [wiki](https://github.com/jibranhaider/explicitSolidDynamics.wiki.git) page.



<br/><br/>
<hr> 

# Applications

### Implosion of a bottle
The simulation can be seen [here](https://youtu.be/K6T6OTjlOzQ).

<img src="/docs/results/implodingBottle/10.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/20.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/25.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/30.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/35.png" width="14%"> &nbsp; &nbsp;
<img src="/docs/results/implodingBottle/40.png" width="14%">

<br/>

### Crushing of a thin cylinder
The simulation can be seen [here](https://youtu.be/7_0i0FlCHtc).

<img src="/docs/results/crushingCylinder/10.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/15.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/20.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/25.png" width="23%">

<img src="/docs/results/crushingCylinder/30.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/35.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/40.png" width="23%"> &nbsp; &nbsp;
<img src="/docs/results/crushingCylinder/50.png" width="23%">




<br/><br/>
<hr> 

# Author
This toolkit is developed and maintained by [Jibran Haider](http://jibranhaider.weebly.com/). 



<br/><br/>
<hr> 

# Citation
An article has been submitted for publication in the Journal of Computer Physics Communications.



<br/><br/>
<hr> 

# License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file. 
