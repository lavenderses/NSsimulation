import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as anime
from decimal import Decimal, ROUND_HALF_UP
import functions as F2d
import functions3d as F3d


class Simulate2D:
    def __init__(self, lx, ly, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8):
        divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'),
                                                          rounding=ROUND_HALF_UP))
        divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'),
                                                          rounding=ROUND_HALF_UP))
        self.ux = np.zeros((divlx + 3, divly + 2))
        self.uy = np.zeros((divlx + 2, divly + 3))
        self.p = np.zeros((divlx + 2, divly + 2))
        self.lx, self.ly = lx, ly
        self.DELT, self.DELL = DELT, DELL
        self.MU, self.RHO, self.EPS = MU, RHO, EPS
        self.imgs = []

        #points
        self.xx = np.array([np.full(divly + 2, i) for i in range(divlx + 2)])
        self.yy = np.tile(np.arange(divly + 2), divlx + 2)

        plt.rcParams['font.size'] = 26
        self.fig = plt.figure(figsize=(20, 14), dpi=100)

    def update(self, hs, v0, theta, time_range, w=None, rad_range=None, CNT_MAX=10000, save=False):
        lx, ly  = self.lx, self.ly
        DELT, DELL = self.DELT, self.DELL
        theta = np.deg2rad(theta)

        dirname = './imgs_ex/2d_{}_{}'.format(v0, theta)
        os.makedirs(dirname, exist_ok=True)
        passed_time, num = 0, 0

        while passed_time < time_range:
            self.ux, self.uy, theta, w = F2d.ff(self.ux, self.uy, self.DELT, hs, v0, theta, w, rad_range)

            #Plot
            ax = self.fig.add_subplot(111,
                                      xlim=(-0.1, self.lx * 1.1),
                                      ylim=(-0.1, self.ly * 1.1),
                                      aspect='equal',)
            ax.set(xlabel='Room Length x/m', ylabel='Room Height y/m')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', '5%', pad='3%')
            im = F2d.plot(self.ux, self.uy, lx, ly, DELT, DELL, DELL, self.xx, self.yy, v0, num, ax)
            self.fig.colorbar(im, cax=cax)
            self.fig.savefig('{}/{:0=10}.png'.format(dirname, num))
            if not save:
                self.imgs.append([im])
            plt.cla()
            plt.clf()

            print('time : {:.2f}'.format(num * DELT))
            #Advection
            ux_ast = F2d.adv_x(self.ux, self.uy, DELT, DELL, DELL)
            uy_ast = F2d.adv_y(self.ux, self.uy, DELT, DELL, DELL)
    
            #Viscosity
            ux_ast, uy_ast = F2d.vis(self.ux, self.uy, ux_ast, uy_ast, DELT, DELL, DELL, self.MU, self.RHO)
    
            #Solce Poisson Eq to Align ux&uy (div u = 0)
            self.p = F2d.poisson(ux_ast, uy_ast, DELT, DELL, DELL, self.p, self.EPS, CNT_MAX)

            #Fix each Velocity
            self.ux, self.uy = F2d.fix_u(ux_ast, uy_ast, self.p, DELT, DELL, DELL, self.RHO)
    
            passed_time = num * DELT
            num += 1

        if not save:
            anim = anime.ArtistAnimation(self.fig, self.imgs, interval=DELT * 1000)
            anim.save('./animation_{}_{}.mp4'.format(v0, theta), writer='ffmpeg')
        

class Simulate3D:
    def __init__(self, lx, ly, lz, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8):
        #Set Constance Values
        divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        divlz = int(Decimal(str(lz / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        self.ux = np.zeros((divlx + 3, divly + 2, divlz + 2))
        self.uy = np.zeros((divlx + 2, divly + 3, divlz + 2))
        self.uz = np.zeros((divlx + 2, divly + 2, divlz + 3))
        self.p = np.zeros((divlx + 2, divly + 2, divlz + 2))
        self.lx, self.ly, self.lz = lx, ly, lz
        self.DELT, self.DELL = DELT, DELL
        self.MU, self.RHO, self.EPS = MU, RHO, EPS
        self.imgs = []
    
        self.xx = np.array([np.full((divlz + 2) * (divly + 2), i) for i in range(divlx + 2)])
        self.yy = np.tile(np.array([np.full(divlz + 2, j) for j in range(divly + 2)]).flatten(), divlx + 2)
        self.zz = np.tile(np.arange(divlz + 2), (divly + 2) * (divlx + 2))

        plt.rcParams['font.size'] = 26
        self.fig = plt.figure(figsize=(24, 14), dpi=100)

    def update(self, hs, v0, theta, phai, time_range, w=None, rad_range=None, CNT_MAX=10000, save=False):
        lx, ly, lz  = self.lx, self.ly, self.lz
        DELT, DELL = self.DELT, self.DELL
        MU, RHO, EPS = self.MU, self.RHO, self.EPS
        phai = np.deg2rad(phai)
        passed_time, num = 0, 0

        theta = np.deg2rad(theta)
        phai = np.deg2rad(phai)

        dirname = './imgs_ex/3d_{}_{}_{}'.format(v0, theta, phai)
        os.makedirs(dirname, exist_ok=True)

        # Particles Array
        ptclsx = np.random.rand(2000) * lx
        ptclsy = np.random.rand(2000) * ly
        ptclsz = np.random.rand(2000) * lz
        particles = np.stack([ptclsx, ptclsy, ptclsz]).T

        while passed_time < time_range:
            #First Condition
            self.ux, self.uy, self.uz = F3d.ff(self.ux, self.uy, self.uz, DELT, hs, v0, theta, phai)

            #Plot
            ax = self.fig.add_subplot(111,
                                      projection='3d',
                                      xlim=(-0.1, lx * 1.1),
                                      ylim=(-0.1, ly * 1.1),
                                      zlim=(-0.1, lz * 1.1),)
            ax.set_xlabel('Room Length x/m', labelpad=20)
            ax.set_ylabel('Room Depth y/m', labelpad=20)
            ax.set_zlabel('Room Height z/m', labelpad=20)
            ax.set_title('t = {:.2f} [s]'.format(num * DELT))

            #img = F3d.plot(self.ux, self.uy, self.uz, lx, ly, lz, DELT, DELL, DELL, DELL, self.xx, self.yy, self.zz, v0, ax)
            img = F3d.p_plot(self.ux, self.uy, self.uz, particles, lx, ly, lz, DELT, DELL, DELL, DELL, v0, ax)
            #self.fig.colorbar(img)
            self.fig.savefig('{}/{:0=10}.png'.format(dirname, num))
            if not save:
                self.imgs.append([img])
            plt.cla()
            plt.clf()

            print('time : {:.2f}'.format(num * DELT))
            #Advection
            ux_ast = F3d.adv_x(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)
            uy_ast = F3d.adv_y(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)
            uz_ast = F3d.adv_z(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)

            #Viscosity
            ux_ast, uy_ast, uz_ast = F3d.vis(self.ux, self.uy, self.uz, ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, MU, RHO)

            #Solce Poisson Eq to Align ux&uy (div u = 0)
            self.p = F3d.poisson(ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, self.p, EPS, CNT_MAX)

            #Fix each Velocity
            self.ux, self.uy, self.uz = F3d.fix_u(ux_ast, uy_ast, uz_ast, self.p, DELT, DELL, DELL, DELL, RHO)

            passed_time = num * DELT
            num += 1

        if not save:
            anim = anime.ArtistAnimation(self.fig, self.imgs, interval=100)
            anim.save('./animation_{}_{}_{}.mp4'.format(v0, theta, phai), writer='ffmpeg')


if __name__ == '__main__':
    # Set Each Constance
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
    time_range = 30.

    rng = np.deg2rad(60)
    rad_range = [-rng, rng]
    T = 18.0
    w = 2 * np.pi / T

    '''
    model = Simulate2D(room_x, room_y)
    model.update(h, v0, the, time_range, w=w, rad_range=rad_range)
    '''
    model = Simulate3D(room_x, room_y, room_z)
    model.update(h_3d, v0, the, pha, time_range)
    #'''
