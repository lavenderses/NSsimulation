# THE SIMULATION OF NAVIER STOKES EQUATION
This is a simulation of Navier Stokes Equation.  
I still keep on modifying these files to opimize, make it faster, greater.  
Let me know when you recognize something wrong or could be better.  
<br>
<br>
## SIMULATION MODEL
I adopted "The Arakawa C-types Staggared Grid" as a frame and "The Fractional Steps Method" to descrete each formula.  
And I also adopted "The Jacobi Method" as a iterative method to solve Poisson Eq.  
<br>
There are some for loops in functions that solve Advection Eqation because of complicated conditions.  
<br>
<br>
## ENVIRONMENT
Machine FMVA77GR  
Ubuntu 20.04.2 LTS  
CPU Intel(R) Core(TM) i7-2670QM CPU @ 2.20GHz  
<br>
anaconda 4.10.1  
numpy 1.19.2  
matplotlib 3.3.2  
<br>
<br>
## HOW TO USE
This simulation is assuming a 2D room. It has height and width.  
Variable lx means its Width and ly does its height.  
When you want to do a 3D room, lx means its widt, ly means its depth and lz means its height.
<br>
I put variables below.  
<br>
<br>
lx,&nbsp;ly,&nbsp;lz&nbsp;&nbsp;&nbsp;: The Lengths of Imaginary Room [m].  
DELT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: The Micro Time. Default is 0.01[s].  
DELL&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: The Micro Length == The Length of The Micro Volume. Default is 0.1[m].  
divlx, divly&nbsp;:&nbsp;Divided Numbers in x, y Directions  
ux, uy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Fluid Velocities in x, y Directions at Point (x, y)[numpy array]  
vx, vy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Temporary Velocities in x, y Directions  
ux_ast&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Calculated Velocity in x Direction at Point (x, y)[numpy array]  
uy_ast&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Calculated Velocity in y Direction at Point (x, y)[numpy array]  
div&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Calculated Divergence at Point (x, y)[numpy array](This Must be Nearly Zero)  
RHO&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Density(Uniform). Default is 1.293[kg/m^3](Air).  
MU&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;The Dinamic Viscosity Coeficient. Default is 1.82e-8[Pa s].  
h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Fan Height.  
v0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;The First Velocity Condition.  
EPS&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;Tiny Error Constance That Evaluate Pressure. Default is 1e-8.  
CNT_MAX&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;The Number That You Wanna Repeat. Default is 10000[times].  
p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;The Presure at Point (x, y).  
fx, fy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;The Forces in x, y Directions at Point (x, y)[numpy array].  
save&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:&nbsp;If You Want to Build Animation File in the Script, Set 'False'. If You Want to Build Animation File after Outputing All Img files, Set 'True'. Default is 'False'.
<br>
<br>
Main script file is "simulate.py", so run this file or import this file to simulate.  
Args in main function is lx, ly, (lz when you want to simulate 3D room), h, v0, theta, time_range and save.  
I put the meanings of each arge above.  
<br>
Output format is a animation file(animation_v0_theta.mp4) written by ffmpeg.  
This will be created in this current directory.  
<br>
<br>
## Thank you for seeing my repository!
<br>
