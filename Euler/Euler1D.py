def partial_x_cd(f,dx,nx):
	dfdx = 1/(2*dx) *(f[2:nx] - f[0:nx-2])
	return dfdx

def partial_x_bd(f,dx,nx):
	dfdx = 1/(dx) * (f[1:nx-1] - f[0:nx-2])  #Stable for positive advection velocity
	return dfdx

def partial_x_fd(f,dx,nx):
	dfdx = 1/(dx) * (f[2:nx] - f[1:nx-1])  #Stable for negative advection velocity
	return dfdx

def partial_x2_cd(f,dx,nx):
	dfdx = 1/(dx**2) * (f[2:nx] - 2*f[1:nx-1] + f[0:nx-2])
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
	#elif (ID_diff_type < 0.0):
	#	f3 = adv_speed*partial_x2_bd(tracer,dx,nx)
	#	return f3
	#elif (ID_diff_type > 0.0):
	#	f3 = adv_speed*partial_x2_fd(tracer,dx,nx)
	# 	return f3
	else:
		print('Please provide right value for ID_diff_type')