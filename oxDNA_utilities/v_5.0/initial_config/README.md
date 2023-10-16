******************************

Tools to generate initial config and external forces for oxDNA simulations.

******************************

`prep_data.py` is the main file. 
One should select the desired `force_opt` inside the file.
The available forces are:
- repulsive sphere (eventually two distinct spheres, one with smaller radius to confine the resources nearer to the box center, see `rep_sphere.py`)
- mutual traps (usually exploited to allow easier formation of H bonds between predator and resource)
- harmonic traps (used as well to confine filaments)/ N.B.: apparently, harmonic traps and caps are in conflict.

N.B.: Python2 required for external forces
