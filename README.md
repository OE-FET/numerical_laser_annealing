# Numerical Simulation of laser annealed 1D stacks of materials
## Overview

This repository contains matlab software to simulate the maximum temperature reached for a stack of materials. In the power point demonstrate we results from using the software to predict the required stack range 


For example this software solved for the electromanetic by s
abosroed energy (organe)
![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/wave_reflections.png)

### Generation of Pulse
The pulse is generated as picture below. Here we assume that the pulse is guassian distributed and centered at a specific frequency.   

![Laser Pulse](https://github.com/OE-FET/numerical_laser_annealing/blob/master/imgs/pulse_generation.png)



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
