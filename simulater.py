import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime
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
---------------------------------------------------------------------
'''


def simulate(Lx, Ly, H, v0, theta, time_range, delt=0.01, dell=0.1, p_rho=1.293, mu=1.82e-5, eps=1e-8, w=1.8, cnt_max=10000):
    #Set Constance Values
    divLx = int(Decimal(str(Lx / dell)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divLy = int(Decimal(str(Ly / dell)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    ux = np.zeros((divLx + 3, divLy + 2))
    uy = np.zeros((divLx + 2, divLy + 3))
    P = np.zeros((divLx + 2, divLy + 2))
    passed_time = 0
    theta = np.deg2rad(theta)
    plt.xlim(-0.1, Lx + 0.1)
    plt.ylim(-0.1, Ly + 0.1)

    #First Condition
    ux, uy = ini.force_and_first(ux, uy, delt, H, v0, theta)

    #First Plot (must be white)
    plts = []
    fig = plt.figure()

    #First Plot (must be white)
    plt_img = plots.plot(ux, uy, divLx, divLy, delt, dell, dell, v0, plt, passed_time)
    plts.append([plt_img])

    while passed_time < time_range:
        print(passed_time)
        #Advection
        ux_ast = adv.advectionX(ux, uy, delt, dell, dell)
        uy_ast = adv.advectionY(ux, uy, delt, dell, dell)

        #Viscosity
        ux_ast, uy_ast = adv.viscosity(ux, uy, ux_ast, uy_ast, delt, dell, dell, mu, p_rho)

        #Solce Poisson Eq to Align ux&uy (div u = 0)
        P = poisson.Poisson(ux_ast, uy_ast, delt, dell, dell, P, eps, w, cnt_max)

        #Fix each Velocity
        ux, uy = poisson.fix_u(ux_ast, uy_ast, P, delt, dell, dell, p_rho)
        ux, uy = ini.force_and_first(ux, uy, delt, H, v0, theta)

        #Plot
        passed_time += delt
        plt_img = plots.plot(ux, uy, divLx, divLy, delt, dell, dell, v0, plt, passed_time)
        plts.append([plt_img])

    ani = anime.ArtistAnimation(fig, plts, interval=100, blit=True)
    ani.save('./animation.mp4', writer="ffmpeg")
#    ani.show()


if __name__ == '__main__':
    #Set Each Constance
    Lx = 4.0
    Ly = 2.5
    H = [[1, 3],
         [1, 4]]
    v0 = 5.0
    theta = 0
    time_range = 20.0

    simulate(Lx, Ly, H, v0, theta, time_range)
