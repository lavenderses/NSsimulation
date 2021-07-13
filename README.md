# THE SIMULATION OF NAVIER STOKES EQUATION
This is a simulation of Navier Stokes Equation.
I still keep on modifying these files to opimize, make it faster, greater.
Let me know when you find something wrong or could be better.

## SIMULATION MODEL
I adopted "The Arakawa C-types Staggared Grid" as a frame and "The Fractional Steps Method" to descrete each formula.
And I also adopted "The Jacobi Method" as a iterative method to solve Poisson Eq.
There are some for loops in functions that solve Advection Eqation because of complicated conditions.

![CodeCogsEqn-2](https://user-images.githubusercontent.com/72122101/125432614-6cf5165c-8475-4483-aea9-38ffe35d4e5e.png)

## ENVIRONMENTS
||Version etc|
----|----
|Machine|FMVA77GR|
|OS|Ubuntu 20.04.2 LTS|
|CPU|Intel(R) Core(TM) i7-2670QM CPU @ 2.20GHz|
|Python|3.8.10|
|conda|4.10.3|
|numpy|1.19.2|
|matplotlib|3.3.2|
## USAGE
This simulation is assuming a 2D room. It has height and width.
Variable `lx` means its Width and `ly` does its height.
When you want to do a 3D room, `lx` means its widt, `ly` means its depth and `lz` means its height.
I put variables below.
|`Variable`|Explanation|
----|----
|`lx, ly, lz`|The Lengths of Imaginary Room [m].|
|`DELT`| The Micro Time. Default is 0.01[s].|
|`DELL`| The Micro Length of The Micro Volume. Default is 0.1[m].|
|`divlx, divly`|Divided Numbers in x, y Directions.|
|`ux, uy`| Fluid Velocities in x, y Directions at Point (x, y)`ndarray`.|
|`vx, vy`|Temporary Velocities in x, y Directions.|
|`ux_ast`|Calculated Velocity in x Direction at Point (x, y)`ndarray`.|
|`uy_ast`|Calculated Velocity in y Direction at Point (x, y)`ndarray`.|  
|`div`|Calculated Divergence at Point (x, y)`ndarray`(This Must be Nearly Zero).|
|`RHO`|Density(Uniform). Default is 1.293[kg/m^3](Air).|
|`MU`|The Dinamic Viscosity Coeficient. Default is 1.82e-8[Pa s].|
|`h`| Fan Height.|
|`v0`|The First Velocity Condition.|
|`EPS`|Tiny Error Constance That Evaluate Pressure. Default is 1e-8. |
|`CNT_MAX`|The Number That You Wanna Repeat. Default is 10000[times]. |
|`p`|The Presure at Point (x, y).|
|`fx, fy`|The Forces in x, y Directions at Point (x, y)`ndarray`.|
|`save`| If You Want to Build Animation File in the Script, Set 'False'. If You Want to Build Animation File after Outputing All Img files, Set 'True'. Default is 'False'.|
Main script file is simulater.py, so run this file or import this file to simulate.  
Args in main function is `lx`, `ly`, (`lz` when you want to simulate 3D room), `h`, `v0`, `theta`, `time_range` and `save`.  
I put the meanings of each arge above.  
Output format is a animation file(animation_v0_theta.mp4) written by ffmpeg.  
This will be created in this current directory.  
## EXSAMPLE
```sh
$ python simulater.py 2D
```
#### Example Result
- - -
![0000011999](https://user-images.githubusercontent.com/72122101/125426234-06c142ab-96f0-415d-8452-c9d174bda9cd.png "Example Output Image.")
If you want to watch Example movie, visit my [Google Drive](https://drive.google.com/file/d/1ElK341AROLEM2WMz7PFNz2mwZGhHhUIx/view?usp=sharing). This is a 30s movie and might be a little bit heavy(9.3MB).

## REFERENCE
https://warp.da.ndl.go.jp/info:ndljp/pid/3481431/www.eto.titech.ac.jp/contents/sub03/chapter01.html
https://www.imi.kyushu-u.ac.jp/PDF/Hirose_20160612.pdf

## Thank you for seeing my repository!
