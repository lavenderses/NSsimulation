def force_and_first(ux, uy, H, v0, delt):
    #GRAVITY
    #uy += 9.8 * delt

    #BOUNDARY CONDITION (v = 0 ON THE WALL)
    #x = 0 and LAST x (LEFT and RIGHT WALL)
    for y in range(ux.shape[1]):
        ux[0][y] = 0
        ux[1][y] = 0
        ux[ux.shape[0] - 2][y] = 0
        ux[ux.shape[0] - 1][y] = 0
    
    for y in range(uy.shape[1]):
        uy[0][y] = 0
        uy[1][y] = 0
        uy[uy.shape[0] - 2][y] = 0
        uy[uy.shape[0] - 1][y] = 0
        
    #y = 0 and LAST y (FLOOR and CEILING)
    for x in range(ux.shape[0]):
        ux[x][0] = 0
        ux[x][1] = 0
        ux[x][ux.shape[1] - 2] = 0
        ux[x][ux.shape[1] - 1] = 0
    
    for x in range(uy.shape[0]):
        uy[x][0] = 0
        uy[x][1] = 0
        uy[x][uy.shape[1] - 2] = 0
        uy[x][uy.shape[1] - 1] = 0

    #FAN
    for h in H:
        ux[h[0]][h[1]] = v0

    return ux, uy