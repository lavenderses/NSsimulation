# THE SIMULATION OF NAVIER STOKES EQUATION
This is a simulation of Navier Stokes Equation.  
I still keep on modifying these files to opimize, make it faster, and to be adopted 3D simulation.  
Let me know when you recognize something wrong or could be better.    


## SIMULATION MODEL
I adopted "Arakawa C-types Staggared Grid" as a frame and "Fractional Steps Method" to descrete each formula.  
And I also adopted "Jacobi Method" as a iterative method to solve Poisson Eq.    

There are some for loops in functions that solve Advection Eqation because of complicated conditions.    


## HOW TO USE
This simulation is assuming a 2D room. It has Height and Width.    
Variable Lx means its Width and Ly does its Height.    

I put variables below.    

Lx, Ly       : The Width and Height of Imaginary Room [m]  
delt         : The Micro Time  
dell         : The Micro Length == The Length of The Micro Volume  
divLx, divLy : Divided Numbers in x, y Directions  
ux, uy       : Fluid Velocities in x, y Directions at Point (x, y)[numpy array]  
vx, vy       : Temporary Velocities in x, y Directions  
ux_ast       : Calculated Velocity in x Direction at Point (x, y)[numpy array]  
uy_ast       : Calculated Velocity in y Direction at Point (x, y)[numpy array]  
div          : Calculated Divergence at Point (x, y)[numpy array](This Must be Nearly Zero)  
p_rho        : Density(Uniform)  
mu           : The Dinamic Viscosity Coeficient  
H            : Fan Height  
v0           : The First Velocity Condition  
eps          : Tiny Error Constance That Needs in SOR  
w            : Accelation Constance  
cnt_max    : The Number That You Wanna Repeat  
P            : The Presure at Point (x, y)  
fx, fy       : The Forces in x, y Directions at Point (x, y)[numpy array]  


Main script file is "simulate.py", so run this file or import this file to simulate.  
Args in main function is Lx, Ly, H, v0, theta, time_range.  
I put the meanings of each arge above.    

Output format is a animation file(animation.mp4). This will be created in this directory.      


Thank you for seeing my repository.  

