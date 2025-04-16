# Falling-Sheet-Simulation

This project simulates the motion of a rigid rectangular sheet (such as a lightweight paper-like object) as it falls under the influence of gravity and air resistance using MATLAB. It captures both translational and rotational dynamics to reproduce realistic behavior like flipping, tumbling, and swaying.

##   Objective

To numerically simulate the descent of a rigid, lightweight sheet from a given height and initial angle using physical laws and differential equations. The simulation tracks:

- Translational motion (position and velocity)
- Rotational motion (angular velocity and orientation)
- Effects of gravity, aerodynamic drag, and aerodynamic torque

##  Physics Model

- **Gravitational Force** pulls the sheet downward.
- **Aerodynamic Drag** opposes the motion and depends on the sheet’s velocity and projected area.
- **Aerodynamic Torque** is generated due to uneven air pressure across the sheet's surface, causing it to spin or flip.

##  Simulation Setup

- Mass: 5 g  
- Dimensions: 20 cm × 15 cm  
- Air Density: 1.225 kg/m³  
- Drag Coefficient: 1.2  
- Torque Coefficient: 0.2  
- Gravity: 9.81 m/s²  
- Time Step: 0.005 s  
- Initial conditions (angle, height) are set by the user.

##  Numerical Method

The simulation uses **Euler Integration** to update:

- Linear position and velocity
- Angular velocity
- Rotation matrix (orientation)

At each timestep, the code calculates the forces and torques, updates the current state of the system, and stops when the sheet hits the ground.

##  Output

The simulation generates a plot showing:

- The 3D trajectory of the falling sheet
- Its orientation and motion in real time
- After completion of the motion it generates various plots like z position vs time, velocity vs time, drag force vs time etc.
- It also saves all the calculated data in a csv file at each instance of 0.02 second.
- It also gives the total time taken to hit the ground from initial height.


##  How to Run

1. Open MATLAB.
2. Clone or download this repository.
3. Run `FallingSheetEuler.m`.
4. Enter the desired drop height and angle when prompted.

##  Notes

- The sheet is assumed to be rigid (not flexible or bending).
- The simulation does not account for collisions with other objects or complex turbulence.

---

Feel free to modify the parameters or improve the integration method (e.g., use Runge-Kutta) for higher accuracy because this code uses Euler Integration which is a 
pretty straighforward numerical analysis technique and is not very accurate for extreme precision modelling.


