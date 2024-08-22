#This code will try and make an ellipse to represent an orbit with classical satellite state representation (see Vallado). 
#It calculates the radius vectors as functions of time. 
#It solves Kepler's equation to work out the eccentric and true anomalies,
#from which it can construct the orbit. 
#You need to give an array of times in hours to the class. 

import numpy as np
import matplotlib.pyplot as plt 
import os

from spacepy import coordinates as coord
from spacepy.time import Ticktock 
import datetime as dt
from matplotlib.patches import Wedge, Polygon, Circle

class ellipse():
	'''This will make an object to describe an ellipse. '''
	
	def __init__(self, t, ra=19, rp=2, inc=70, raan=80, omega=300, ptime=dt.datetime(2024,8,10)):
		'''This takes in the six parameters to describe the orbit. 
		
		Parameters
		----------
		t - time in hours. Should be an array from 0 to expected period.  
		rp - radius of perigee (RE)
		ra - radius of apogee (RE) 
		p - semi parameter (RE) - Calculated from rp and ra.
		e - eccentricity (0-1) - Calculated from rp and ra. 
		inc - inclination (deg)
		raan - right angle of ascending node (RAAN) (deg)
		omega - argument of periapsis (deg)
		ptime - periapsis time as a datetime object. 
		
		'''
		
		self.RE = 6370000
		self.ptime = ptime
		
		#Calculate eccentricity and semiparameter in metres. 
		self.rp = rp*self.RE
		self.ra = ra*self.RE
		
		self.e = (self.ra-self.rp)/(self.ra+self.rp) 
		
		self.p = (2*self.rp*self.ra)/(self.ra+self.rp) 
		print ('Eccentricity = ', self.e)
		print ('Semiparameter = ', self.p/self.RE) 
		
		
		self.plot_path = os.environ.get("PLOT_PATH")+"orbits/"
		
		self.inc = np.deg2rad(inc) 
		self.raan = np.deg2rad(raan)
		self.omega = np.deg2rad(omega)
		self.t = t*3600 
		
		#Store mu value for Earth. 
		self.mu = 6.67e-11 * 6e24 
		
		#Calculate semi-major axis. 
		self.a = self.p/(1 - self.e**2)
		
		#Calculate semi-minor axis. 
		self.b = self.a*((1 - self.e**2)**0.5) 
		
		#Calculate specific angular momentum of satellite. 
		self.h = np.sqrt(self.mu*self.p) 
		
		#Calculate period of satellite. 
		self.period = 2*np.pi*self.a*self.b/self.h 
		
		#Calculate the mean anomaly. 
		self.M = self.t*(self.mu/(self.a**3))**0.5
		
		#Next you need to get the eccentric anomaly by solving Kepler's equation
		#using the Newton-Raphson method outlined in Vallado. 
		self.solve_kepler()
		
		#Next calculate the true anomaly. 
		self.nu = 2*np.arctan2(((1+self.e)**0.5)*np.tan(self.EA/2),((1-self.e)**0.5))
		inward = np.where(self.nu < 0 ) 
		self.nu[inward] = self.nu[inward] + 2*np.pi
		
		#Calculate the radius and velocity vectors. 
		self.get_radius_vectors() 
		
		#Run the checks on the angles. 
		self.check_angles() 
		
		#Get datetimes for each time from periapsis time. 
		self.get_datetime_from_periapsis()
		
		#Convert GSE to GSM coordinates. 
		self.gse_to_gsm()
		
	def solve_kepler(self, tol=1e-5):
		'''This will solve Kepler's equation using the 
		Newton-Raphson method outlined in Vallado'''
		
		print ("Solving Kepler's equation...")
		self.EA = np.zeros(self.M.size)
		
		for m in range(len(self.M)):
			#Select initial guess depending on value of M. 
			if self.M[m] < np.pi:
				E0 = self.M[m] + self.e
			else:
				E0 = self.M[m] - self.e 
			
			#Set incorrect diff just to make it run. 
			diff = 0.1 
			while diff > tol: 
				E1 = E0 + (self.M[m] - E0 + self.e*np.sin(E0))/(1 - self.e*np.cos(E0))
				diff = abs(E1-E0)
				E0=E1
			
			self.EA[m] = E1
	
	def get_radius_vectors(self):
		'''This will contain all the calculations for the radius and velocity vectors, 
		and the rotations to put them in the correct orientation.'''
		
		#Calculate radius of satellite in m. 
		self.r = self.p/(1 + self.e*np.cos(self.nu)) 
		
		#Calculate speed of satellite in m/s. 
		self.v = np.sqrt(self.mu*(2/self.r - 1/self.a))
		
		#Calculate the flight path angle. 
		self.phi = np.zeros((self.r.size))
		self.phi[self.nu < np.pi] = np.arccos(self.h/(self.r[self.nu < np.pi]*self.v[self.nu < np.pi]))
		self.phi[self.nu >= np.pi] = -np.arccos(self.h/(self.r[self.nu >= np.pi]*self.v[self.nu >= np.pi]))
		
		
		#Get initial x, y and z values. 
		self.x0 = self.r*np.cos(self.nu)
		self.y0 = self.r*np.sin(self.nu)
		self.z0 = self.r*np.zeros(self.x0.size)
		
		#Get initial periapsis vector. 
		self.xp0 = self.rp
		self.yp0 = 0
		self.zp0 = 0 
		
		#Get initial semiparameter vector. 
		self.xs0 = 0
		self.ys0 = self.p
		self.zs0 = 0
		
		#Rotate ellipse around z axis by argument of periapsis. 
		self.x1 = self.x0*np.cos(self.omega) - self.y0*np.sin(self.omega)
		self.y1 = self.x0*np.sin(self.omega) + self.y0*np.cos(self.omega) 
		self.z1 = self.z0 
		
		#Rotate periapsis vector around z axis by argument of periapsis. 
		self.xp1 = self.xp0*np.cos(self.omega) - self.yp0*np.sin(self.omega)
		self.yp1 = self.xp0*np.sin(self.omega) + self.yp0*np.cos(self.omega)
		self.zp1 = self.zp0 
		
		#Rotate periapsis vector around z axis by argument of periapsis. 
		self.xs1 = self.xs0*np.cos(self.omega) - self.ys0*np.sin(self.omega)
		self.ys1 = self.xs0*np.sin(self.omega) + self.ys0*np.cos(self.omega)
		self.zs1 = self.zs0
		
		#Rotate ellipse around x axis by inclination. 
		self.x2 = self.x1
		self.y2 = self.y1*np.cos(self.inc) - self.z1*np.sin(self.inc)
		self.z2 = self.y1*np.sin(self.inc) + self.z1*np.cos(self.inc) 
		
		#Rotate periapsis vector around x by inclination. 
		self.xp2 = self.xp1
		self.yp2 = self.yp1*np.cos(self.inc) - self.zp1*np.sin(self.inc)
		self.zp2 = self.yp1*np.sin(self.inc) + self.zp1*np.cos(self.inc)
		
		#Rotate periapsis vector around x by inclination. 
		self.xs2 = self.xs1
		self.ys2 = self.ys1*np.cos(self.inc) - self.zs1*np.sin(self.inc)
		self.zs2 = self.ys1*np.sin(self.inc) + self.zs1*np.cos(self.inc)
		
		#Rotate ellipse around z axis by RAAN. 
		self.x_gse = self.x2*np.cos(self.raan) - self.y2*np.sin(self.raan)
		self.y_gse = self.x2*np.sin(self.raan) + self.y2*np.cos(self.raan) 
		self.z_gse = self.z2 
		
		#Rotate periapsis vector around z by RAAN. 
		self.xp_gse = self.xp2*np.cos(self.raan) - self.yp2*np.sin(self.raan)
		self.yp_gse = self.xp2*np.sin(self.raan) + self.yp2*np.cos(self.raan) 
		self.zp_gse = self.zp2
		
		#Rotate periapsis vector around z by RAAN. 
		self.xs_gse = self.xs2*np.cos(self.raan) - self.ys2*np.sin(self.raan)
		self.ys_gse = self.xs2*np.sin(self.raan) + self.ys2*np.cos(self.raan) 
		self.zs_gse = self.zs2
		
		#Create radius vectors. 
		self.r_vector = np.array((self.x_gse, self.y_gse, self.z_gse)).T 
		
		#Create periapsis vector and unit vector. 
		self.periapsis_vector = np.array([self.xp_gse, self.yp_gse, self.zp_gse])
		self.periapsis_unit_vector = self.periapsis_vector/self.rp
		
		#Create periapsis vector and unit vector. 
		self.semip_vector = np.array([self.xs_gse, self.ys_gse, self.zs_gse])
		self.semip_unit_vector = self.semip_vector/self.p
		
		#Get magnitude of new radius vectors. Should be same as r. 
		self.r_mag = np.sqrt(self.x_gse**2 + self.y_gse**2 + self.z_gse**2)
		
		#Get r unit vector. 
		self.r_unit_vector = np.array([self.r_vector[i]/self.r_mag[i] for i in range(self.r_mag.size)])
		
		#Get normal vector. 
		self.n_unit_vector = np.cross(self.periapsis_unit_vector, self.semip_unit_vector)
		
		#Get h vector. 
		self.h_vector = self.h*self.n_unit_vector 
		
		#Get tangential vector for each r vector. 
		self.tangent_unit = np.array([np.cross(self.n_unit_vector, self.r_unit_vector[i]) for i in range(self.r_mag.size)]) 
		
		#Get velocity unit vector. 
		self.v_unit_vector = np.array([np.cos(self.phi[i])*self.tangent_unit[i] + np.sin(self.phi[i])*self.r_unit_vector[i] for i in range(self.r_mag.size)]) 
		
		#Get velocity vector. 
		self.v_vector = np.array([self.v[i]*self.v_unit_vector[i] for i in range(self.r_mag.size)])
		
	def check_angles(self): 
		'''This will check that the angles of the final ellipse match the inclination, 
		argument of perigee and right ascension of the ascending node. '''
		
		print ('Checks on final ellipse...')
		
		#Input i, j and k vectors. 
		self.I = [1,0,0]
		self.J = [0,1,0]
		self.K = [0,0,1] 
		
		#Check inclination angle. Correct.  
		cosi = np.dot(self.K, self.n_unit_vector)
		inclination = np.rad2deg(np.arccos(cosi))
		print ('Inclination = {:.2f}'.format(inclination)) 
		
		#Get line of nodes vector and get RAAN. 
		self.nodes = np.cross(self.K, self.n_unit_vector) 
		self.nodes_mag = np.sqrt(self.nodes[0]**2 + self.nodes[1]**2 + self.nodes[2]**2) 
		self.nodes_unit = self.nodes/self.nodes_mag 
		#print ('Nodes = ', self.nodes_unit)
		cossigma = np.dot(self.I, self.nodes_unit)
		sigma = np.rad2deg(np.arccos(cossigma))
		print ('RAAN = {:.2f}'.format(sigma)) 
		
		#Check argument of perigee. 
		cosomega = np.dot(self.nodes_unit, self.r_unit_vector[0]) 
		omega = np.rad2deg(np.arccos(cosomega))
		print ('Arg. Perigee = {:.2f}'.format(omega)) 
		
	def get_datetime_from_periapsis(self):
		'''This will work out a datetime object for each point around the orbit by adding the time to the datetime object of the periapsis point.'''
    	
		self.dt_list = []
    	
		for t in range(len(self.t)):
			deltat = dt.timedelta(seconds=float(self.t[t])) 
			self.dt_list.append(self.ptime+deltat) 
    		
    
	def gse_to_gsm(self):
		'''This will use spacepy to convert from GSE to GSM'''
		
		#Need coordinates in appropriate form for conversion. 
		self.coords_gse = np.zeros((self.x_gse.size,3))
		self.coords_gse[:,0] = self.x_gse
		self.coords_gse[:,1] = self.y_gse
		self.coords_gse[:,2] = self.z_gse 
    	
		#Needs to take data in a a list of xyz points. 
		coord_obj = coord.Coords(self.coords_gse, 'GSE', 'car')
    	
		#Add time information. 
		coord_obj.ticks = Ticktock(self.dt_list, 'UTC') 
    	
		#To convert to GSM. 
		self.coords_gsm = coord_obj.convert('GSM', 'car') 
		
		
	def plot_ellipse_3d(self, lims=(-20,20), elev=30, azim=30, scatter=True):
		'''This plots the ellipse in 3D space.'''
		
		#Set conditions for inbound and outbound. 
		outward = np.where(self.phi > 0)
		inward = np.where(self.phi < 0)
		
		fig = plt.figure(figsize=(8,6))
		ax1 = fig.add_subplot(121, projection='3d')
	
		#FINAL GSE. 
		if scatter:
			ax1.scatter(self.x_gse[outward]/self.RE, self.y_gse[outward]/self.RE, self.z_gse[outward]/self.RE, c='cyan', marker='x', s=5)
			ax1.scatter(self.x_gse[inward]/self.RE, self.y_gse[inward]/self.RE, self.z_gse[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax1.plot(self.x_gse[outward]/self.RE, self.y_gse[outward]/self.RE, self.z_gse[outward]/self.RE, 'cyan')
			ax1.plot(self.x_gse[inward]/self.RE, self.y_gse[inward]/self.RE, self.z_gse[inward]/self.RE, 'darkblue')
		
		#Add periapsis vector after RAAN. 
		ax1.plot([0,self.periapsis_vector[0]/self.RE], [0, self.periapsis_vector[1]/self.RE], [0, self.periapsis_vector[2]/self.RE], 'k')
		
		
		ax1.set_xlabel(r'$x_{GSE}$')
		ax1.set_ylabel(r'$y_{GSE}$')
		ax1.set_zlabel(r'$z_{GSE}$')
		ax1.set_title('Orbit in GSE')
		
		#Add axes into plot.  
		ax1.plot(lims, (0,0), (0,0), 'k')
		ax1.plot((0,0), lims, (0,0), 'k')
		ax1.plot((0,0), (0,0), lims, 'k')
		
		ax1.set_xlim(lims)
		ax1.set_ylim(lims)
		ax1.set_zlim(lims)
		ax1.set_aspect('equal')
		
		self.add_earth(ax1)
		
		ax1.view_init(elev,azim) 
		
		#Add GSM position to main graph. 
		ax2 = fig.add_subplot(122, projection='3d')
		
		if scatter:
			ax2.scatter(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.y[outward]/self.RE, self.coords_gsm.z[outward]/self.RE, c='cyan', marker='x', s=5)
			ax2.scatter(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, self.coords_gsm.z[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax2.plot(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.y[outward]/self.RE, self.coords_gsm.z[outward]/self.RE, 'cyan')
			ax2.plot(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, self.coords_gsm.z[inward]/self.RE, 'darkblue')
		
		ax2.set_xlabel(r'$x_{GSM}$')
		ax2.set_ylabel(r'$y_{GSM}$')
		ax2.set_zlabel(r'$z_{GSM}$')
		ax2.set_title('Orbit in GSM')
		
		#Add axes into plot.  
		ax2.plot(lims, (0,0), (0,0), 'k')
		ax2.plot((0,0), lims, (0,0), 'k')
		ax2.plot((0,0), (0,0), lims, 'k')
		
		ax2.set_xlim(lims)
		ax2.set_ylim(lims)
		ax2.set_zlim(lims)
		ax2.set_aspect('equal')
		
		self.add_earth(ax2)
		
		ax2.view_init(elev,azim) 
		
		fig.text(0.95, 0.05, 'Outbound', color='cyan', fontsize=8, ha='right')
		fig.text(0.95, 0.03, 'Inbound', color='darkblue', fontsize=8, ha='right')
		
		fig.savefig(self.plot_path+"time_ellipses_3d.png") 
		
	def plot_ellipse_2d(self, scatter=False):
		'''This will plot the xy and xz planes of the orbit in GSE and GSM. 
		
		Parameters
		----------
		scatter - Boolean to plot a scatter plot or a line plot. 
		'''
		
		#Set conditions for inbound and outbound. 
		outward = np.where(self.phi > 0)
		inward = np.where(self.phi < 0)
		
		fig = plt.figure(figsize=(6,6))
		fig.subplots_adjust(wspace=0.3)
		#GSE
		ax1 = fig.add_subplot(221)
		if scatter:
			ax1.scatter(self.x_gse[outward]/self.RE, self.z_gse[outward]/self.RE, c='cyan', marker='x', s=5)
			ax1.scatter(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax1.plot(self.x_gse[outward]/self.RE, self.z_gse[outward]/self.RE, 'cyan')
			ax1.plot(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, 'darkblue')
		
		ax1.set_xlabel(r'$x_{GSE}$')
		ax1.set_ylabel(r'$z_{GSE}$') 
		ax1.set_title('Orbit in GSE')
		self.make_earth_2d(ax1, rotation=-90)
		ax1.set_aspect('equal')
		
		ax2 = fig.add_subplot(223)
		if scatter:
			ax2.scatter(self.x_gse[outward]/self.RE, self.y_gse[outward]/self.RE, c='cyan', marker='x', s=5)
			ax2.scatter(self.x_gse[inward]/self.RE, self.y_gse[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax2.plot(self.x_gse[outward]/self.RE, self.y_gse[outward]/self.RE, 'cyan')
			ax2.plot(self.x_gse[inward]/self.RE, self.y_gse[inward]/self.RE, 'darkblue')
		ax2.set_xlabel(r'$x_{GSE}$')
		ax2.set_ylabel(r'$y_{GSE}$') 
		#ax2.set_title('Orbit in GSE')
		self.make_earth_2d(ax2, rotation=-90)
		ax2.set_aspect('equal')
		
		ax3 = fig.add_subplot(222)
		if scatter:
			ax3.scatter(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.z[outward]/self.RE, c='cyan', marker='x', s=5)
			ax3.scatter(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.z[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax3.plot(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.z[outward]/self.RE, 'cyan')
			ax3.plot(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.z[inward]/self.RE, 'darkblue')
		ax3.set_xlabel(r'$x_{GSM}$')
		ax3.set_ylabel(r'$z_{GSM}$') 
		ax3.set_title('Orbit in GSM')
		self.make_earth_2d(ax3, rotation=-90)
		ax3.set_aspect('equal')
		
		ax4 = fig.add_subplot(224)
		if scatter:
			ax4.scatter(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.y[outward]/self.RE, c='cyan', marker='x', s=5)
			ax4.scatter(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, c='darkblue', marker='x', s=5)
		else:
			ax4.plot(self.coords_gsm.x[outward]/self.RE, self.coords_gsm.y[outward]/self.RE, 'cyan')
			ax4.plot(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, 'darkblue')
		ax4.set_xlabel(r'$x_{GSM}$')
		ax4.set_ylabel(r'$y_{GSM}$') 
		#ax4.set_title('Orbit in GSM')
		self.make_earth_2d(ax4, rotation=-90)
		ax4.set_aspect('equal')
		
		fig.text(0.95, 0.05, 'Outbound', color='cyan', fontsize=8, ha='right')
		fig.text(0.95, 0.03, 'Inbound', color='darkblue', fontsize=8, ha='right')
		
		fig.savefig(self.plot_path+'time_ellipses_2d.png')
			
	def plot_orbital_parameters(self):
		'''This will plot how various parameters change with true anomaly.
		
		Parameters
		----------
		t - Boolean to plot as a function of time instead of true anomaly. '''
		
		#nu_ticks = np.linspace(0,360,13)
		#t_ticks = np.linspace(0,self.period/3600,13)
		
		fig = plt.figure(figsize=(6,8))
		fig.subplots_adjust(hspace=0.5, left=0.15)
		
		#Radial position (RE)
		ax1 = fig.add_subplot(411)
		ax1.plot(self.t/3600, (self.r/self.RE), 'k')
		ax1.set_xlim(0, self.period/3600)
		ax1.set_ylabel('Radius of Orbit (RE)') 
		ax1.grid(which='both')
		ax1.set_title('Orbital Parameters') 
		ax1.set_ylim(0,)
		ax1.minorticks_on()
		
		#Speed (km/s)
		ax2 = fig.add_subplot(412)
		ax2.plot(self.t/3600, self.v/1000, 'k')
		ax2.set_xlim(0, self.period/3600)
		ax2.set_ylabel('Speed (km/s)')
		ax2.grid(which='both')
		ax2.set_ylim(0,)
		ax2.minorticks_on()
		
		#Flight path angle (deg). 
		ax3 = fig.add_subplot(413)
		ax3.plot(self.t/3600, np.rad2deg(self.phi), 'k')
		ax3.set_xlim(0, self.period/3600)
		ax3.set_ylabel('Flight Path \nAngle (deg)') 
		ax3.grid(which='both')
		ax3.minorticks_on()
		
		#Eccentric Anomaly  
		ax4 = fig.add_subplot(414)
		ax4.plot(self.t/3600, np.rad2deg(self.EA), 'k')
		ax4.set_xlim(0, self.period/3600)
		ax4.set_xlabel('Time (hrs)')
		ax4.set_ylabel('Eccentric\n Anomaly (deg)') 
		ax4.grid(which='both')
		ax4.minorticks_on()
		ax4.set_ylim(0,360)
		ax4.set_yticks([0,90,180,270,360])
		fig.savefig(self.plot_path+"time_orbital_params.png") 
		
		
	def add_earth(self, ax):
		'''This will add a sphere for the Earth. '''
		
		#Create a spherical surface. 
		radius = 1
		u = np.linspace(0, 2*np.pi, 100) 
		v = np.linspace(0, np.pi, 100) 
		x = radius* np.outer(np.cos(u), np.sin(v))
		y = radius* np.outer(np.sin(u), np.sin(v))
		z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

		ax.plot_surface(x, y, z, color='k', lw=0, alpha=1)
		
	def make_earth_2d(self, ax, rotation=0):
		'''This will add a little plot of the Earth on top for reference. '''

		# Add white circle first. 
		r=1
		circle = Circle((0,0), r, facecolor='w', edgecolor='navy')
		ax.add_patch(circle)

		# Add nightside. 
		theta2 = np.arange(181)-180+rotation
		xval2 = np.append(r*np.cos(theta2*(np.pi/180)),0)
		yval2 = np.append(r*np.sin(theta2*(np.pi/180)),0)
		verts2 = [[xval2[i],yval2[i]] for i in range(len(xval2))]
		
		polygon2 = Polygon(verts2, closed=True, edgecolor='navy', facecolor='navy', alpha=1) 
		ax.add_patch(polygon2)	
