# Flight_Controller

## Summary
This project was done as part of "ME 354: Automatic Control of Aerospace Vehicles" to design a flight controller for an ongoing research project at Lehigh University. The subject aircraft was called the Gun-Launched Unmanned Aerial System (GLUAS) and was designed to be fired from a mortar shell, unfold a set of 10 feather-shaped wings from the main body, and autonomously survey unknown landscapes.

<p align="center">
  <img width="340" src="Images/printed_gluas.PNG"> 
  <img width="435" src="Images/cad_top_gluas.PNG"> 
</p> 

<p align="center">
  <img width="600" src="Images/cad_front_gluas.PNG"> 
</p> 

## Athena Vortex Lattice
MIT's Athena Vortex Lattice (AVL) was used in approximating the aircraft's stability derivatives using the following files:
* surface_geometry.avl
* flight_conditions.run
* body.mass

<p align="center">
  <img width="500" src="Images/avl.PNG"> 
</p> 

## State Space Model
Longitudinal and lateral state space models were constructed based on the moment stability and force stability derivatives. System behavior was assessed with an elevator doublet response and aileron input response and tuned to reduce oscillation.

## Controller
A pitch-attitude control loop and heading hold control loop were implemented using a combination of PI compensators, lead compensators, and air speed-based transfer functions.
