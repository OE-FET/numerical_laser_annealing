# Numerical simulation of laser pulse annealing of 1D stack
## Overview

This repository contains Matlab software to simulate the laser pulse annealing of a stack of materials. Here the simulation is used to determine the maximum temperature achieved during laser annealing in order to estimate if the stack will be obliterated or deformed by the laser pulse. Using this software, we designed a stack of materials to locally heat solution processed IZO films on polyethylene naphthalate (PEN) substrates without did damaging the PEN.

### Pulse generation
Here the pulse is assumed to equivalent to a single frequency electromagnetic (EM) wave spatially confined by spatially multiplying the wave by a Gaussian. Through Fourier analysis the software decomposes the pulse into a linear combination of single frequency EM waves. The resulting interaction of each independent wave is solved separately and then linearly recombined to approximate the laser pulse.

![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/pulse_generation.png)

### Solving electromagnetic wave equation

At each boundary the amplitude of the EM wave equation must be continuous along with its first derivate. Where n is the nth layer in the stack of material and L/R are the corresponding waves in the nth layer moving left/right respectively. Using these two equations one can solve for the left/right moving waves in each layer given the initial incident wave. For a solid intro to EM wave dynamics please reference "Introduction to Electrodynamics" by David J. Griffiths.

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?E_%7Bn%2CL%7D%20&plus;E_%7Bn%2CR%7D%20%3D%20E_%7Bn&plus;1%2CL%7D%20&plus;E_%7Bn&plus;1%2CR%7D">
</p>

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20x%7D%28E_%7Bn%2CL%7D&plus;E_%7Bn%2CR%7D%29%20%3D%5Cfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20x%7D%28E_%7Bn&plus;1%2CL%7D&plus;E_%7Bn&plus;1%2CR%7D%29">
</p>



The goal of the EM simulation is to solve for the absorption rate (Ar), also known as the rate of work, as a function of time in each layer. This is achieved by substituting the solved EM wave into the following equations. 

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?D%20%3D%20%5Cvarepsilon%20E">
</p>

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?%5Cbigtriangledown%20E%20%3D%20-%5Cfrac%7B%5Cpartial%20B%7D%7B%5Cpartial%20t%7D">
</p>


<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?H%20%3D%20%5Cfrac%7BB%7D%7B%5Cmu%20%7D">
</p>


<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?U%20%3D%20%5Cfrac%7B1%7D%7B2%7D%28E%5Ccdot%20D&plus;B%5Ccdot%20H%29%29">
</p>

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?S%20%3D%20E%20%5Ctimes%20H">
</p>

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?Ar%20%3D%20%5Cfrac%7B%5Cpartial%20W%7D%7B%5Cpartial%20t%7D%20%3D%20J%5Ccdot%20E%20%3D%20-%5Cfrac%7B%5Cpartial%20U%7D%7B%5Cpartial%20t%7D-%5Cbigtriangledown%20%5Ccdot%20S">
</p>



A four layer example stack (including Air, In4ZnO, SiO2 and Si) with the solved EM wave:
![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/wave_reflections.png)
### Thermal diffusion finite element

Explicitly solving the diffusion equation is challenging due to the large number of terms present in the energy absorption rate. To circumvent this obstacle we use finite element to approximate the diffusion equation. The diffusion equation can be set below. Here the boundary conditions of the finite element simulation are held at a fixed temperature (room temperature). For the same stack seen above, the finite element results are presented below:

<p align="center">
  <img width="" height="" src="https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20t%7D%3D%5Cvarrho%20%5Ccdot%20C_%7Bp%7D%5Ccdot%20%5Cfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D-%5Cbigtriangledown%20%28k%5Cbigtriangledown%20T%29">
</p>



![Diffusion Equation](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/Finite_element_results.png)
### Questions, problems, collaborations?
For questions or problems please create an [issue](https://github.com/OE-FET/numerical_laser_annealing/issues). Please raise issues if sections of code require further explanation. For any other extensions please contact me directly (jwarmitage@gmail.com). Always happy to chat. :D

## How to implement
### Requirements
- standard Matlab 2016a
Should be a standard license package installation. 

## Components
- **Run_example.m** - main script for running example software.
- **Run_EM.m** - runs an EM simulation
- **Run_diffusion_eqn.m** - runs a thermal simulation
- **quickInt.m** - speeds up integrations of EM fields
- **InitConstants.m** - initializes physical constants
- **genAR.m** - generate absorption of electromagnetic fields
- **Example_stack.m** - example stack of materials
- **diffusion_eqn.m** - solves thermal diffusion equation 

## Authors:
The software is written by [John Armitage](https://github.com/jwarmitage).
