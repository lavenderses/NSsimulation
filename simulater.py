import numpy as np
from decimal import Decimal, ROUND_HALF_UP
from functions import advection_viscosity as adv
from functions import initialization as ini
from functions import poisson
from functions import plots


'''
-------variables-----------------------------------------------------
Lx, Ly       : The Lengths of Imaginary Room [m]
delt         : The Micro Time
dell         : The Micro Length == The Length of The Micro Volume
divLx, divLy : Divided Numbers in x, y Directions
ux, uy       : Fluid Velocities in x, y Directions at Point (x, y)[numpy array]
vx, vy       : Temporary Velocities in x, y Directions
ux_calced    : Calculated Velocity in x Direction at Point (x, y)[numpy array]
uy_calced    : Calculated Velocity in y Direction at Point (x, y)[numpy array]
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
---------------------------------------------------------------------
'''


def simulate(Lx, Ly, H, v0, theta, time_range, delt=0.01, dell=0.1, p_rho=1.2, mu=1.82e-5, eps=1e-8, w=1.8, cnt_max=10000):
    #Set Constance Values
    divLx = int(Decimal(str(Lx / dell)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divLy = int(Decimal(str(Ly / dell)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    ux = np.zeros((divLx + 3, divLy + 2))
    uy = np.zeros((divLx + 2, divLy + 3))
    div = np.zeros((divLx, divLy))
    P = np.zeros((divLx, divLy))
    plts = []
    num = 0
    passed_time = 0

    #First Plot (must be white)
    plts = plots.plot(ux, uy, Lx, Ly, delt, dell, dell, plts, v0)

    while passed_time < time_range:
        print(passed_time)
        num += 1

        #First Condition
        ux, uy = ini.force_and_first(ux, uy, H, v0, theta)

        #Advection
        ux = adv.advectionX(ux, uy, delt, dell, dell)
        uy = adv.advectionY(ux, uy, delt, dell, dell)

        #Viscosity
        ux, uy = adv.viscosity(ux, uy, mu, delt, dell, dell, p_rho)

        #Solce Poisson Eq to Align ux&uy (div u = 0)
        P = poisson.Poisson(ux, uy, div, delt, dell, dell, P, eps, w, cnt_max)

        #FIX THE VELOCITY
        ux, uy = poisson.fix_u(ux, uy, P, delt, dell, dell, p_rho)
        ux, uy = ini.force_and_first(ux, uy, H, v0, theta)

        #PLOT
        plts = plots.plot(ux, uy, Lx, Ly, delt, dell, dell, plts, v0)

        passed_time += delt


if __name__ == '__main__':
    #SET EACH CONSTANCE
    Lx = 4.0
    Ly = 2.5
    H = [[1, 3],
         [1, 4]]
    v0 = 5.0
    theta = 30
    time_range = 60.0

    simulate(Lx, Ly, H, v0, theta, time_range)
