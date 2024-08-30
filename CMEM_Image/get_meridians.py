#This will work out the XZ and XY planes from a 3D dataset. 
import numpy as np 

def calculate_meridian_planes(x_3d, y_3d, z_3d, var_3d):
	'''This will actually work out the XZ and XY plane data properly by taking means between the nearest planes
		
	Parameters
	----------
	x_3d - 3D array of x positions
	y_3d - 3D array of y positions
	z_3d - 3D array of z positions 
	var_3d - 3D array of data (e.g emissivity)
		
	Returns
	-------
	xp_y - 2D array of x values in XZ plane. 
	yp_y - 2D array of y values in XZ plane. 
	zp_y - 2D array of z values in XZ plane. 
	var_y - 2D array of data values in XZ plane. 
	xp_z - 2D array of x values in XY plane. 
	yp_z - 2D array of y values in XY plane. 
	zp_z - 2D array of z values in XY plane. 
	var_z - 2D array of data values in XY plane. 
		
	''' 
			
	#XZ plane. 
	#Get lowest positive y value and highest negative y value. 
	i_yl = np.where(y_3d[0,:,0] < 0, y_3d[0,:,0], -np.inf).argmax()
	i_yu = np.where(y_3d[0,:,0] > 0, y_3d[0,:,0], np.inf).argmin()
	
	xp_yl = x_3d[:,i_yl]
	xp_yu = x_3d[:,i_yu]
	xp_y = (xp_yl+xp_yu)/2.
		
	yp_yl = y_3d[:,i_yl]
	yp_yu = y_3d[:,i_yu]
	yp_y = (yp_yl+yp_yu)/2.
		
	zp_yl = z_3d[:,i_yl]
	zp_yu = z_3d[:,i_yu]
	zp_y = (zp_yl+zp_yu)/2.
		
	var_yl = var_3d[:,i_yl]
	var_yu = var_3d[:,i_yu]
	var_y = (var_yl+var_yu)/2.
		
	#XY plane. 
	#Get lowest positive z value and highest negative z value. 
	i_zl = np.where(z_3d[:,0,0] < 0, z_3d[:,0,0], -np.inf).argmax()
	i_zu = np.where(z_3d[:,0,0] > 0, z_3d[:,0,0], np.inf).argmin()
		
	xp_zl = x_3d[i_zl]
	xp_zu = x_3d[i_zu]
	xp_z = (xp_zl+xp_zu)/2. 
		
	yp_zl = y_3d[i_zl]
	yp_zu = y_3d[i_zu]
	yp_z = (yp_zl+yp_zu)/2. 
		
	zp_zl = z_3d[i_zl]
	zp_zu = z_3d[i_zu]
	zp_z = (zp_zl+zp_zu)/2. 
		
	var_zl = var_3d[i_zl]
	var_zu = var_3d[i_zu]
	var_z = (var_zl+var_zu)/2. 
		
	return xp_y, yp_y, zp_y, var_y, xp_z, yp_z, zp_z, var_z 
		
