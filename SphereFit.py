def SphereFit(coords):
    #from http://stackoverflow.com/questions/15785428/how-do-i-fit-3d-data
    import numpy as np
    from scipy.optimize import leastsq

    p0 = [0, 0, 0, 1]

    def fitfunc(p, coords):
        x0, y0, z0, R = p
        x, y, z = coords.T
        return np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)

    errfunc = lambda p, x: fitfunc(p, x) - p[3]

    p1, flag = leastsq(errfunc, p0, args=(coords,))
    
    #p1: (x, y, z, radius)
    return p1[-1]