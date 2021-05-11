import numpy as np

def partial_x_cd(f,dx,nx):

    """This function computes the first order central difference of f(x,y) wrt x
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell in x-direction
    nx: Number of points in the domain in x-direction
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """
    
    dfdx = np.zeros_like(f)
    
    dfdx[1:nx-1,:] = 1/(2*dx) * (f[2:nx,:]-f[0:nx-2,:])
    dfdx[-1,:] = 1/(2*dx) *(f[0,:] - f[-2,:])
    dfdx[0,:] = 1/(2*dx) *(f[1,:] - f[-1,:])
    
    return dfdx

def partial_x_bd(f,dx,nx):
    
    """This function computes the first order backward derivative of f(x,y) wrt x
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell in x-direction
    nx: Number of points in the domain in x-direction
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """

    dfdx = np.zeros_like(f)
    
    dfdx[1:nx,:] = 1/(dx) *  (f[1:nx,:] - f[0:nx-1,:])
    dfdx[0,:] = 1/(dx) * (f[0,:] - f[-1,:])
    
    return dfdx

def partial_x_fd(f,dx,nx):
    
    """This function computes the first order forward derivative of f(x,y) wrt x

    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell in x-direction
    nx: Number of points in the domain in x-direction
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """

    dfdx = np.zeros_like(f)
    
    dfdx[0:nx-1,:] = 1/(dx) * (f[1:nx,:] - f[0:nx-1,:])
    dfdx[-1,:] = 1/(dx) * (f[0,:] - f[-1,:])

    return dfdx

def partial_y_cd(f,dy,ny):

    """This function computes the first order central difference of f(x,y) wrt y
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dy: Width of each cell in y-direction
    ny: Number of points in the domain in y-direction
    -------------------------------------------------------------------------------------
    Returns:
    fdy: Partial derivative of f
    """
    
    dfdy = np.zeros_like(f)
    
    dfdy[:,1:ny-1] = 1/(2*dy) * (f[:,2:ny] - f[:,0:ny-2])
    dfdy[:,-1] = 1/(2*dy) *(f[:,0] - f[:,-2])
    dfdy[:,0] = 1/(2*dy) *(f[:,1] - f[:,-1])
    
    return dfdy

def partial_y_bd(f,dy,ny):
    
    """This function computes the first order backward derivative of f(x,y) wrt y
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dy: Width of each cell in y-direction
    ny: Number of points in the domain in y-direction
    -------------------------------------------------------------------------------------
    Returns:
    dfdy: Partial derivative of f
    """

    dfdy = np.zeros_like(f)
    
    dfdy[:,1:ny] = 1/(dy) *  (f[:,1:ny] - f[:,0:ny-1])
    dfdy[:,0] = 1/(dy) * (f[:,0] - f[:,-1])
    
    return dfdy

def partial_y_fd(f,dy,ny):
    
    """This function computes the first order forward derivative of f(x,y) wrt y

    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dy: Width of each cell in y-direction
    ny: Number of points in the domain in y-direction
    -------------------------------------------------------------------------------------
    Returns:
    dfdy: Partial derivative of f
    """

    dfdy = np.zeros_like(f)
    
    dfdy[:,0:ny-1] = 1/(dy) * (f[:,1:ny] - f[:,0:ny-1])
    dfdy[:,-1] = 1/(dy) * (f[:,0] - f[:,-1])

    return dfdx

def adv_x(adv_speed,tracer,dx,nx,ID_adv = 0):

    """Evaluates advection term in x-direction for a 2D geometry

    -------------------------------------------------------------------------------------
    Arguments:
    adv_speed: x-directed advection velocity
    tracer: The tracer that needs to be advected
    dx: Width of each cell in x-direction
    nx: Number of points in the domain in x-direction
    ID_adv: ID to run forward (>0), backward(<0) or central(=0) difference scheme
    -------------------------------------------------------------------------------------
    Returns:
    f3: Advected value of the tracer (in x-direction)
    """

    if (ID_adv == 0): 
        f3 = adv_speed*partial_x_cd(tracer,dx,nx)
        return f3
    elif (ID_adv < 0.0):
        f3 = adv_speed*partial_x_bd(tracer,dx,nx)
        return f3
    elif (ID_adv > 0.0):
        f3 = adv_speed*partial_x_fd(tracer,dx,nx)
        return f3
    else:
        print('Please provide right value for ID_diff_type')


def adv_y(adv_speed,tracer,dy,ny,ID_adv = 0):

    """Evaluates advection term in y-direction for a 2D geometry

    -------------------------------------------------------------------------------------
    Arguments:
    adv_speed: y-directed advection velocity
    tracer: The tracer that needs to be advected
    dy: Width of each cell in y-direction
    ny: Number of points in the domain in y-direction
    ID_adv: ID to run forward (>0), backward(<0) or central(=0) difference scheme
    -------------------------------------------------------------------------------------
    Returns:
    f3: Advected value of the tracer (in x-direction)
    """

    if (ID_adv == 0): 
        f3 = adv_speed*partial_y_cd(tracer,dy,ny)
        return f3
    elif (ID_adv < 0.0):
        f3 = adv_speed*partial_y_bd(tracer,dy,ny)
        return f3
    elif (ID_adv > 0.0):
        f3 = adv_speed*partial_y_fd(tracer,dy,ny)
        return f3
    else:
        print('Please provide right value for ID_diff_type')
