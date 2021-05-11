import numpy as np

def partial_x_cd(f,dx,nx):

    """This function computes the first order central difference of f(x)
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell
    nx: Number of points in the domain
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """
    
    dfdx = np.zeros_like(f)
    for j in np.arange(1,f.size-1,1):
        dfdx[j] = 1/(2*dx) * (f[j+1]-f[j-1])    #Backward difference
    
    dfdx[-1] = 1/(2*dx) *(f[0] - f[-2])
    dfdx[0] = 1/(2*dx) *(f[1] - f[-1])
    return dfdx

def partial_x_bd(f,dx,nx):
    
    """This function computes the first order backward derivative of x
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell
    nx: Number of points in the domain
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """

    dfdx = np.zeros_like(f)
    for j in range(f.size):
        dfdx[j] = 1/(dx) * (f[j]-f[j-1])    #Backward difference

    return dfdx

def partial_x_fd(f,dx,nx):
    
    """This function computes the first order forward derivative of x
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell
    nx: Number of points in the domain
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Partial derivative of f
    """
    
    dfdx = np.zeros_like(f)
    for j in range(f.size - 1):
        dfdx[j] = 1/(dx) * (f[j+1]-f[j])    #Forward difference
    
    dfdx[f.size] = 1/(dx) * (f[0] - f[f.size])

    return dfdx

def partial_x2_cd(f,dx,nx):
    
    """This function computes the second order central derivative of x
    
    -------------------------------------------------------------------------------------
    Arguments:
    f: Function which needs to be differentiated
    dx: Width of each cell
    nx: Number of points in the domain
    -------------------------------------------------------------------------------------
    Returns:
    dfdx: Second derivative of f(x)
    """

    dfdx = np.zeros_like(f)
    #Interior points
    dfdx[1:nx-1] = 1/(dx**2) * (f[2:nx] - 2*f[1:nx-1] + f[0:nx-2])
    
    #Boundary points
    dfdx[0] = 1/(dx**2) * (f[1] - 2*f[0] + f[-1])
    dfdx[-1] = 1/(dx**2) * (f[0] - 2*f[-1] + f[-2])
    
    return dfdx


def adv_x(adv_speed,tracer,dx,nx,ID_diff_type = 0):
	
	if (ID_diff_type == 0): 
		f3 = adv_speed*partial_x_cd(tracer,dx,nx)
		return f3
	elif (ID_diff_type < 0.0):
		f3 = adv_speed*partial_x_bd(tracer,dx,nx)
		return f3
	elif (ID_diff_type > 0.0):
		f3 = adv_speed*partial_x_fd(tracer,dx,nx)
		return f3
	else:
		print('Please provide right value for ID_diff_type')

def diff_x(kappa,tracer,dx,nx,ID_diff_type = 0):
    if (ID_diff_type == 0): 
        f3 = kappa*partial_x2_cd(tracer,dx,nx)
        return f3
    else:
        print('Please provide right value for ID_diff_type')