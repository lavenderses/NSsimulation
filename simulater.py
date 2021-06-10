import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime
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
lx, ly, lz   : The Lengths of Imaginary Room [m]
DELT         : The Micro Time
DELL         : The Micro Length == The Length of The Micro Volume
divlx, divly : Divided Numbers in x, y Directions
ux, uy       : Fluid Velocities in x, y Directions at Point (x, y)[numpy array]
vx, vy       : Temporary Velocities in x, y Directions
ux_ast       : Calculated Velocity in x Direction at Point (x, y)[numpy array]
uy_ast       : Calculated Velocity in y Direction at Point (x, y)[numpy array]
div          : Calculated Divergence at Point (x, y)[numpy array](This Must be Nearly Zero)
RHO        : Density(Uniform)
MU           : The Dinamic Viscosity Coeficient
h            : Fan Height
v0           : The First Velocity Condition
EPS          : Tiny Error Constance That Needs in SOR
CNT_MAX      : The Number That You Wanna Repeat
p            : The Presure at Point (x, y)
fx, fy       : The Forces in x, y Directions at Point (x, y)[numpy array]
---------------------------------------------------------------------
'''


def simulate(lx, ly, hs, v0, theta, time_range, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8, CNT_MAX=10000, save=False):
    #Set Constance Values
    divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'), rounding=ROUND_HALF_UP))
    ux = np.zeros((divlx + 3, divly + 2))
    uy = np.zeros((divlx + 2, divly + 3))
    p = np.zeros((divlx + 2, divly + 2))
    passed_time = 0
    num = 0
    imgs = []
    dirname = './imgs2d_{}_{}'.format(v0, theta)
    theta = np.deg2rad(theta)

    #points
    xx = np.array([np.full(divly + 2, i) for i in range(divlx + 2)])
    yy = np.tile(np.arange(divly + 2), divlx + 2)


    if save:
        os.makedirs(dirname, exist_ok=True)
        ax = None
    else:
        fig = plt.figure(figsize=(20, 14), dpi=100)
        ax = fig.add_subplot(111,
                             projection='3d',
                             xlim=(-0.1, lx * 1.1),
                             ylim=(-0.1, ly * 1.1))
        ax.set(xlabel='x', ylabel='y')

    while passed_time < time_range:
        #First Condition
        ux, uy = ini.force_and_first(ux, uy, DELT, hs, v0, theta)

        #Plot
        im = plots.plot(ux, uy, lx, ly, DELT, DELL, DELL, xx, yy, v0, passed_time, num, ax, dirname=dirname)
        if not save and num == 0:
            fig.colorbar(im)
        imgs.append([im])

        print('time : ', passed_time)
        #Advection
        ux_ast = adv.advection_x(ux, uy, DELT, DELL, DELL)
        uy_ast = adv.advection_y(ux, uy, DELT, DELL, DELL)

        #Viscosity
        ux_ast, uy_ast = adv.viscosity(ux, uy, ux_ast, uy_ast, DELT, DELL, DELL, MU, RHO)

        #Solce Poisson Eq to Align ux&uy (div u = 0)
        p = poisson.poisson(ux_ast, uy_ast, DELT, DELL, DELL, p, EPS, CNT_MAX)

        #Fix each Velocity
        ux, uy = poisson.fix_u(ux_ast, uy_ast, p, DELT, DELL, DELL, RHO)

        passed_time += DELT
        num += 1

    if not save:
        anim = anime.ArtistAnimation(fig, imgs, interval=100)
        anim.save('./animation_{}_{}.mp4'.format(v0, theta), writer='ffmpeg')



def simulate3d(lx, ly, lz, hs, v0, theta, phai, time_range, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8, CNT_MAX=10000, save=False):
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
    dirname = './imgs3d_{}_{}_{}'.format(v0, theta, phai)
    theta = np.deg2rad(theta)
    phai = np.deg2rad(phai)
    imgs = []

    xx = np.array([np.full((divlz + 2) * (divly + 2), i) for i in range(divlx + 2)])
    yy = np.tile(np.array([np.full(divlz + 2, j) for j in range(divly + 2)]).flatten(), divlx + 2)
    zz = np.tile(np.arange(divlz + 2), (divly + 2) * (divlx + 2))

    if save:
        os.makedirs(dirname, exist_ok=True)
        ax = None
    else:
        fig = plt.figure(figsize=(20, 14), dpi=100)
        ax = fig.add_subplot(111,
                             projection='3d',
                             xlim=(-0.1, lx * 1.1),
                             ylim=(-0.1, ly * 1.1),
                             zlim=(-0.1, lz * 1.1))
        ax.set(xlabel='x', ylabel='y', zlabel='z')

    while passed_time < time_range:
        #First Condition
        ux, uy, uz = ini3d.force_and_first(ux, uy, uz, DELT, hs, v0, theta, phai)

        #Plot
        img = plots3d.plot(ux, uy, uz, lx, ly, lz, DELT, DELL, DELL, DELL, xx, yy, zz, v0, passed_time, num, ax, dirname=dirname)
        if not save and num == 0:
            fig.colorbar(img)
            imgs.append([img])

        print('time : ', passed_time)
        #Advection
        ux_ast = adv3d.advection_x(ux, uy, uz, DELT, DELL, DELL, DELL)
        uy_ast = adv3d.advection_y(ux, uy, uz, DELT, DELL, DELL, DELL)
        uz_ast = adv3d.advection_z(ux, uy, uz, DELT, DELL, DELL, DELL)

        #Viscosity
        ux_ast, uy_ast, uz_ast = adv3d.viscosity(ux, uy, uz, ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, MU, RHO)

        #Solce Poisson Eq to Align ux&uy (div u = 0)
        p = poisson3d.poisson(ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, p, EPS, CNT_MAX)

        #Fix each Velocity
        ux, uy, uz = poisson3d.fix_u(ux_ast, uy_ast, uz_ast, p, DELT, DELL, DELL, DELL, RHO)

        passed_time += DELT
        num += 1

    if not save:
        anim = anime.ArtistAnimation(fig, imgs, interval=100)
        anim.save('./animation_{}_{}_{}.mp4'.format(v0, theta, phai), writer='ffmpeg')


if __name__ == '__main__':
    #Set Each Constance
    room_x = 4.0
    room_y = 2.0
    room_z = 2.5
    h = [[1, 3],
         [1, 4]]
    h_3d = [[1, 3, 3],
            [1, 3, 4],
            [1, 4, 3],
            [1, 4, 4]]
    v0 = 5.0
    the = 0
    pha = 0
    time_range = 60.

    simulate(room_x, room_y, h, v0, the, time_range, save=True)
    #simulate3d(room_x, room_y, room_z, h_3d, v0, the, pha, time_range, save=True)
