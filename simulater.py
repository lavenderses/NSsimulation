import numpy as np
import matplotlib.pyplot as plt
'''import matplotlib.animation as anime'''
from decimal import Decimal, ROUND_HALF_UP
from functions import advection_viscosity as adv
from functions import initialization as ini
from functions import poisson
from functions import plots
from functions3d import advection_viscosity3d as adv3d
from functions3d import initialization3d as ini3d
from functions3d import poisson3d
from functions3d import plots3d

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
cnt_max      : The Number That You Wanna Repeat
P            : The Presure at Point (x, y)
fx, fy       : The Forces in x, y Directions at Point (x, y)[numpy array]
---------------------------------------------------------------------
'''


def simulate(lx, ly, hs, v0, theta, time_range, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8, W=1.8, CNT_MAX=10000):
    #Set Constance Values
    divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    ux = np.zeros((divlx + 3, divly + 2))
    uy = np.zeros((divlx + 2, divly + 3))
    P = np.zeros((divlx + 2, divly + 2))
    passed_time = 0
    num = 0
    theta = np.deg2rad(theta)

    #First Condition
    ux, uy = ini.force_and_first(ux, uy, DELT, hs, v0, theta)

    #First Plot (must be white)
    plots.plot(ux, uy, divlx, divly, DELT, DELL, DELL, v0, passed_time, num)
    '''plt.colorbar(plt_img)
    plts.append([plt_img])'''

    while passed_time < time_range:
        print(passed_time)
        num += 1
        #Advection
        ux_ast = adv.advectionX(ux, uy, DELT, DELL, DELL)
        uy_ast = adv.advectionY(ux, uy, DELT, DELL, DELL)

        #Viscosity
        ux_ast, uy_ast = adv.viscosity(ux, uy, ux_ast, uy_ast, DELT, DELL, DELL, MU, RHO)

        #Solce Poisson Eq to Align ux&uy (div u = 0)
        P = poisson.Poisson(ux_ast, uy_ast, DELT, DELL, DELL, P, EPS, W, CNT_MAX)

        #Fix each Velocity
        ux, uy = poisson.fix_u(ux_ast, uy_ast, P, DELT, DELL, DELL, RHO)
        ux, uy = ini.force_and_first(ux, uy, DELT, hs, v0, theta)

        #Plot
        passed_time += DELT
        plots.plot(ux, uy, divlx, divly, DELT, DELL, DELL, v0, passed_time, num)

    '''ani = anime.ArtistAnimation(fig, plts, interval=100, blit=True)
    ani.save('./animation.mp4', writer="ffmpeg")'''



def simulate3d(lx, ly, lz, hs, v0, theta, phai, time_range, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8, W=1.8, CNT_MAX=10000):
    #Set Constance Values
    divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divlz = int(Decimal(str(lz / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    ux = np.zeros((divlx + 3, divly + 2, divlz + 2))
    uy = np.zeros((divlx + 2, divly + 3, divlz + 2))
    uz = np.zeros((divlx + 2, divly + 2, divlz + 3))
    p = np.zeros((divlx + 2, divly + 2, divlz + 2))
    passed_time = 0
    num = 0
    theta = np.deg2rad(theta)
    phai = np.deg2rad(phai)

    #First Condition
    ux, uy, uz = ini3d.force_and_first(ux, uy, uz, DELT, hs, v0, theta, phai)

    #First Plot (must be white)
    fig = plt.figure(figsize=(18, 13), dpi=100)
    #plots3d.plot(fig, ux, uy, uz, lx, ly, lz, DELT, DELL, DELL, DELL, v0, passed_time, num)

    while passed_time < time_range:
        print('time : ', passed_time)
        num += 1
        #Advection
        ux_ast = adv3d.advection_x(ux, uy, uz, DELT, DELL, DELL, DELL)
        uy_ast = adv3d.advection_y(ux, uy, uz, DELT, DELL, DELL, DELL)
        uz_ast = adv3d.advection_z(ux, uy, uz, DELT, DELL, DELL, DELL)

        #print('A', np.max(ux_ast))
        #Viscosity
        ux_ast, uy_ast, uz_ast = adv3d.viscosity(ux, uy, uz, ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, MU, RHO)

        #print('V', np.max(ux_ast))
        #Solce Poisson Eq to Align ux&uy (div u = 0)
        p = poisson3d.poisson(ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, p, EPS, W, CNT_MAX)
        #print('P', np.max(p))

        #Fix each Velocity
        ux, uy, uz = poisson3d.fix_u(ux_ast, uy_ast, uz_ast, p, DELT, DELL, DELL, DELL, RHO)
        #print('F', np.max(ux), np.max(p))
        ux, uy, uz = ini3d.force_and_first(ux, uy, uz, DELT, hs, v0, theta, phai)

        #print('L', np.max(ux), np.max(uy), np.max(uz))
        #Plot
        passed_time += DELT
        plots3d.plot(fig, ux, uy, uz, lx, ly, lz, DELT, DELL, DELL, DELL, v0, passed_time, num)
        plt.savefig('../imgs/plt{}.png'.format(num))
        plt.show()



if __name__ == '__main__':
    #Set Each Constance
    room_x = 4.0
    room_y = 2.0
    room_z = 2.5
    h = [[1, 3],
         [1, 4]]
    h_3d = [[1, 1, 3],
         [1, 1, 4]]
    v0 = 5.0
    the = 60
    pha = 30
    time_range = 60.0

#    simulate(room_x, room_y, h, v0, theta, time_range)
    simulate3d(room_x, room_y, room_z, h_3d, v0, the, pha, time_range)
