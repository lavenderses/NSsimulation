import numpy as np
import matplotlib.pyplot as plt
import time
DELTAT=0.5#Δt
OMEGA=1.8#SORの加速係数

#移流。壁のところは速度を更新してないことに注意
def Advection(x,y,vx,vy):
    vx_after = np.zeros((x + 3, y + 2), dtype=np.float64)  # x速度
    vy_after = np.zeros((x + 2, y + 3), dtype=np.float64)  # y速度
    #x方向の移流
    for i in range(2,x+1):
        for j in range(1,y+1):
            u = vx[i, j]
            v = (vy[i - 1,j] + vy[i,j] + vy[i - 1,j + 1] + vy[i,j + 1]) / 4
            if (u >= 0.0) & (v >= 0.0):
                vx_after[i,j] = vx[i,j] - u * (vx[i,j] - vx[i-1,j]) * DELTAT - v * (vx[i,j]-vx[i,j-1]) * DELTAT
            if (u >= 0.0) & (v < 0.0):
                vx_after[i,j] = vx[i,j] - u * (vx[i,j] - vx[i-1,j]) * DELTAT - v * (vx[i,j+1]-vx[i,j]) * DELTAT
            if (u < 0.0) & (v >= 0.0):
                vx_after[i,j] = vx[i,j] - u * (vx[i+1,j] - vx[i,j]) * DELTAT - v * (vx[i,j]-vx[i,j-1]) * DELTAT
            if (u < 0.0) & (v < 0.0):
                vx_after[i,j] = vx[i,j] - u * (vx[i+1,j] - vx[i,j]) * DELTAT - v * (vx[i,j+1]-vx[i,j]) * DELTAT
    #y方向の移流
    for i in range(1,x+1):
        for j in range(2,y+1):
            u = (vx[i,j - 1] + vx[i + 1,j - 1] + vx[i,j] + vx[i+1,j]) / 4
            v = vy[i,j]
            if (u >= 0.0) & (v >= 0.0):
                vy_after[i,j] = vy[i,j] - u * (vy[i,j] - vy[i-1,j]) * DELTAT - v * (vy[i,j]-vy[i,j-1]) * DELTAT
            if (u >= 0.0) & (v < 0.0):
                vy_after[i,j] = vy[i,j] - u * (vy[i,j] - vy[i-1,j]) * DELTAT - v * (vy[i,j+1]-vy[i,j]) * DELTAT
            if (u < 0.0) & (v >= 0.0):
                vy_after[i,j] = vy[i,j] - u * (vy[i+1,j] - vy[i,j]) * DELTAT - v * (vy[i,j]-vy[i,j-1]) * DELTAT
            if (u < 0.0) & (v < 0.0):
                vy_after[i,j] = vy[i,j] - u * (vy[i+1,j] - vy[i,j]) * DELTAT - v * (vy[i,j+1]-vy[i,j]) * DELTAT
    vx[:,:] = vx_after.copy()
    vy[:,:] = vy_after.copy()
    return
def Div(x,y,vx,vy,s):
    for i in range(1,x+1):
        for j in range(1,y+1):
            s[i,j]=( -vx[i,j] -vy[i,j] +vx[i+1,j] +vy[i,j+1] )/DELTAT
    return
def Poisson(x,y,p,s):
    for loop in range(30):#適当に30回ループで切り上げることにする
        for i in range(1,x+1):
            for j in range(1,y+1):
                #もし壁なら、p[i,j]の圧力を代入。壁情報をbool型配列で管理しておくといろんな壁が再現できる
                if i==1:#左の壁
                    p[i-1,j]=p[i,j]
                if i==x:#右の壁
                    p[i+1,j]=p[i,j]
                if j==1:#上の壁
                    p[i,j-1]=p[i,j]
                if j==y:#下の壁
                    p[i,j+1]=p[i,j]
                p[i,j] = (1.0 - OMEGA) * p[i,j] + OMEGA / 4 * (p[i - 1,j] + p[i + 1,j] + p[i,j - 1] + p[i,j + 1] - s[i,j])
    return

def Rhs(x,y,vx,vy,p):
    for i in range(1,x+1):
        for j in range(1,y+1):
            vx[i,j] -= (p[i,j] - p[i - 1,j]) * DELTAT
            vy[i,j] -= (p[i,j] - p[i,j - 1]) * DELTAT
    return

#粒子座標の速度を抽出して座標更新
#スタガード格子なのでxとy速度場の参照がずれる
def Flowparticles(vx,vy,prt):
    for i in range(prt.shape[0]):
        xx=np.clip(prt[i,0],0.0,vx.shape[0]-2)
        yy=np.clip(prt[i,1],0.0,vy.shape[1]-2)
        ixx = np.int32(xx)
        iyy = np.int32(yy-0.5)#スタガード格子なのでいろいろずれる
        sxx = xx - ixx
        syy = (yy-0.5) - iyy
        spdx = (((1.0 - sxx) * vx[ixx,iyy] + sxx * vx[ixx+1,iyy]) * (1.0 - syy) + (
                    (1.0 - sxx) * vx[ixx,iyy+1] + sxx * vx[ixx+1,iyy+1]) * syy) * DELTAT
        ixx = np.int32(xx-0.5)#スタガード格子なのでいろいろずれる
        iyy = np.int32(yy)
        sxx = (xx-0.5) - ixx
        syy = yy - iyy
        spdy = (((1.0 - sxx) * vy[ixx,iyy] + sxx * vy[ixx+1,iyy]) * (1.0 - syy) + (
                    (1.0 - sxx) * vy[ixx,iyy+1] + sxx * vy[ixx+1,iyy+1]) * syy) * DELTAT
        prt[i, 0] += spdx
        prt[i, 1] += spdy
    return







def CFD_plot(x,y):
    vx = np.zeros((x + 3, y + 2), dtype=np.float64)  # x速度
    vy = np.zeros((x + 2, y + 3), dtype=np.float64)  # y速度
    p = np.zeros((x + 2, y + 2), dtype=np.float64)  # 圧力
    s = np.zeros((x + 2, y + 2), dtype=np.float64)  # ダイバージェンス
    prt = np.random.rand(1024, 2) * np.array((x,y),dtype=np.float64)+1  # 粒子の初期座標、x,y
    scat = plt.scatter(prt[:,0],prt[:,1])#初期化的に一度plotしなければならない.そのときplotしたオブジェクトを受け取る受け取る必要がある．
    frame = 0
    while True:
        #移流
        Advection(x, y, vx, vy)
        #粘性は今回は考えない・・・
        #外力
        vx[1, 4]= .5#適当に速度固定
        vx[1, 3]= .5#適当に速度固定
        #ダイバージェンス計算
        Div(x, y, vx, vy, s)
        #圧力計算
        Poisson(x, y, p, s)
        #修正
        Rhs(x, y, vx, vy, p)
        #ここで全部の速度、圧力が次のステップに更新されたことになるので可視化
        Flowparticles(vx, vy, prt)
        # plot 逐次画面更新
        scat.set_offsets(prt)
        plt.pause(0.015) #sec
        # print
        frame += 1
        if frame%16==0:
            print(frame)


if __name__ == "__main__":
    CFD_plot(6,8)
