#This will resample the PPMLR data into the SMILE FOV and work out an image. 

import numpy as np
import matplotlib.pyplot as plt

class ppmlr_image():
	'''This class takes in the ppmlr simulation object and the smile fov object and calculates an image through the simulation.'''
	
	def __init__(self, ppmlr, smile):
		'''This takes in the ppmlr and smile objects.'''
		
		self.ppmlr = ppmlr
		self.smile = smile 
		
	def get_weighted_eta(self, px, py, pz):
		'''This will get the nearest x, y and z values for a given point. 
		
		Parameters
		----------
		px - x coordinate of point in Smile Fov
		py - y coordinate of point in Smile Fov
		pz - z coordinate of point in Smile Fov 
		
		'''
		
		#For now, just do one pixel and do the first point. 
		px = self.smile.xpos[0][0][50]
		py = self.smile.ypos[0][0][50]
		pz = self.smile.zpos[0][0][50]
		print (px, py, pz)
		
		
		#Make empty eta array to fill in. 
		peta = np.zeros((self.smile.xpos[0][0][0].shape))
		
		#If px, py or pz are outside the boundaries of the simulation, eta is zero. 
		if (px < self.ppmlr.x[0]) or (px > self.ppmlr.x[-1]) or (py < self.ppmlr.y[0]) or (py > self.ppmlr.y[-1]) or (pz < self.ppmlr.z[0]) or (pz > self.ppmlr.z[-1]):
			peta = 0 
			return peta
			
		else: 
			#It must be inside the cube.
			
			#Get the indices of x0, y0 and z0.  
			ix = self.get_x0_x1(px) 
			iy = self.get_y0_y1(py)
			iz = self.get_z0_z1(pz) 
			print (ix, iy, iz)
			#Get the x, y, z and eta values for the vertices. 
			self.vert_x = self.ppmlr.x_3d[iz:iz+2, iy:iy+2, ix:ix+2]
			self.vert_y = self.ppmlr.y_3d[iz:iz+2, iy:iy+2, ix:ix+2]
			self.vert_z = self.ppmlr.z_3d[iz:iz+2, iy:iy+2, ix:ix+2]
			self.etav = self.ppmlr.eta_3d[iz:iz+2, iy:iy+2, ix:ix+2]
			
			#Get the radial distances to each vertex. 
			r_vertices = self.get_r_to_vertex(px, py, pz, self.vert_x, self.vert_y, self.vert_z)
		
			#Get weights for each vertex.
			weights = 1/r_vertices 
			
			#Calculate the weighted eta value for the coordinates (px, py, pz) 
			peta = (weights*self.etav).sum()/weights.sum()
			
			return peta 
	
	def get_r_to_vertex(self, px, py, pz, vx, vy, vz):
		'''This calculates the radial distance to a vertex. '''
		
		#r = np.sqrt((vertex[0]-px)**2 + (vertex[1]-py)**2 + (vertex[2]-pz)**2)
		r = np.sqrt((vx-px)**2 + (vy-py)**2 + (vz-pz)**2)
		return r 
		
		
	def get_x0_x1(self, px):
		'''This will get the index of x0'''
		
		diff = self.ppmlr.x - px
		ix = diff[diff <= 0].argmax()	
		return ix
		
		
	def get_y0_y1(self, py):
		'''This will get the index of y0'''
			
		diff = self.ppmlr.y - py
		iy = diff[diff <= 0].argmax()	
		return iy
		
		
	def get_z0_z1(self, pz):
		'''This will get the index of z0'''
		
		diff = self.ppmlr.z - pz
		iz = diff[diff <= 0].argmax()	
		return iz
