# Numerical Simulation of laser annealed 1D stacks of materials
## Overview

This repository contains matlab software to simulate the laser pulse annealing of a stack of materials. Here the simulation is used to determine the maximum temperature achieved during laser annealing inorder to estimate if the stack will be obliterated or deformed by the laser pulse. Using this software, we desgined a stack of materials to locally heat solution processed IZO films on Polyethylene naphthalate substrates that did not damage the substrates.

### Generation of Pulse
Here the pulse is assumed to equalivent to a single freqnecy EM wave spaitally confied by spaiatlay mulitplying the wave by a gusiassn. To simulate the resulting pulse, through Fourier analysis simulation then decomposed the pulse into a linear combination of single freqney EM waves. The resulting interaction of each indepdent wave is solved seperatly.

![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/pulse_generation.png)

### Solving electromagnetic wave equation

The interaction between 

https://latex.codecogs.com/gif.latex?E_%7Bn%2CL%7D%20&plus;E_%7Bn%2CR%7D%20%3D%20E_%7Bn&plus;1%2CL%7D%20&plus;E_%7Bn&plus;1%2CR%7D

https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20E_%7Bn%2CL%7D&plus;E_%7Bn%2CR%7D%7D%7B%5Cpartial%20x%7D%20%3D%20%5Cfrac%7B%5Cpartial%20E_%7Bn&plus;1%2CL%7D&plus;E_%7Bn&plus;1%2CR%7D%7D%7B%5Cpartial%20x%7D

https://latex.codecogs.com/gif.latex?D%20%3D%20%5Cvarepsilon%20E

https://latex.codecogs.com/gif.latex?%5Cbigtriangledown%20E%20%3D%20-%5Cfrac%7B%5Cpartial%20B%7D%7B%5Cpartial%20t%7D

https://latex.codecogs.com/gif.latex?H%20%3D%20%5Cfrac%7BB%7D%7B%5Cmu%20%7D

https://latex.codecogs.com/gif.latex?U%20%3D%20%5Cfrac%7B1%7D%7B2%7D%28E%5Ccdot%20D&plus;B%5Ccdot%20H%29%29

https://latex.codecogs.com/gif.latex?S%20%3D%20E%20%5Ctimes%20H

https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20A%7D%7B%5Cpartial%20t%7D%20%3D-%20%5Cfrac%7B%5Cpartial%20U%7D%7B%5Cpartial%20t%7D%20-%20%5Cbigtriangledown%20%5Ccdot%20S

https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20t%7D%3D%5Cvarrho%20%5Ccdot%20C_%7Bp%7D%5Ccdot%20%5Cfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D-%5Cbigtriangledown%20%28k%5Cbigtriangledown%20T%29

For a solid intro to electromanetic wave dynamics please reference "Introduction to Electrodynamics" by david j. griffiths. 
![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/wave_reflections.png)




### Questions, problems, collaborations?
For questions or problems please create an [issue](https://github.com/OE-FET/numerical_laser_annealing/issues). Please raise issues if sections of code require further explanation. For any other extensions please contact me directly (jwarmitage@gmail.com). Always happy to chat. :D

## How to implement
### Requirements
- standard matlab
Should be a standard installation. Please raise issue if mnecessary. 

## Components
- **Run_example.m** - Main script for running example software.
- **Run_EM.m** - Runs a electromagnetic simulation
- **Run_diffusion_eqn.m** - Runs a thermal simulation
- **quickInt.m** - speeds up intergrations of EM feilds
- **PlotLaser.m** - plots laser 
- **InitConstants.m** - initlaizes physical constants
- **genAR.m** - generate absorption of electromagnetic fields
- **Example_stack.m** - example stack of materials
- **diffusion_eqn.m** - solves thermal diffusion eqn

## Authors:
The software is written by [John Armitage](https://github.com/jwarmitage).
