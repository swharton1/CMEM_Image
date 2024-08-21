#This code will try and make an ellipse to represent an orbit with classical satellite state representation (see Vallado). 

import numpy as np
import matplotlib.pyplot as plt 
import os

from spacepy import coordinates as coord
from spacepy.time import Ticktock 
import datetime as dt

class ellipse():
	'''This will make an object to describe an ellipse. '''
	
	def __init__(self, nu, ra=19, rp=2, inc=70, raan=80, omega=300, ptime=dt.datetime(2024,8,10)):
		'''This takes in the six parameters to describe the orbit. 
		
		Parameters
		----------
		rp - radius of perigee (RE)
		ra - radius of apogee (RE) 
		p - semi parameter (RE) - Calculated from rp and ra.
		e - eccentricity (0-1) - Calculated from rp and ra. 
		inc - inclination (deg)
		raan - right angle of ascending node (RAAN) (deg)
		omega - argument of periapsis (deg)
		nu - true anomaly (deg) - Should be an array from 0 to 359 degrees. 
		ptime - periapsis time as a datetime object. 
		
		'''
		
		self.RE = 6370000
		self.ptime = ptime
		
		#Calculate eccentricity and semiparameter. 
		self.rp = rp*self.RE
		self.ra = ra*self.RE
		
		self.e = (self.ra-self.rp)/(self.ra+self.rp) 
		
		self.p = (2*self.rp*self.ra)/(self.ra+self.rp) 
		print ('Eccentricity = ', self.e)
		print ('Semiparameter = ', self.p/self.RE) 
		
		
		self.plot_path = os.environ.get("PLOT_PATH")+"orbits/"
		
		#self.p = p*self.RE
		#self.e = e
		self.inc = np.deg2rad(inc) 
		self.raan = np.deg2rad(raan)
		self.omega = np.deg2rad(omega)
		self.nu = np.deg2rad(nu) 
		
		#Store mu value for Earth. 
		mu = 6.67e-11 * 6e24 
		
		#Calculate semi-major axis. 
		self.a = self.p/(1 - self.e**2)
		
		#Calculate semi-minor axis. 
		self.b = self.a*((1 - self.e**2)**0.5) 
		
		#Calculate radius of satellite in m. 
		self.r = self.p/(1 + self.e*np.cos(self.nu)) 
		
		#Calculate specific angular momentum of satellite. 
		self.h = np.sqrt(mu*self.p) 
		
		#Calculate period of satellite. 
		self.period = 2*np.pi*self.a*self.b/self.h 
		
		#Calculate speed of satellite in m/s. 
		self.v = np.sqrt(mu*(2/self.r - 1/self.a))
		
		#Calculate the flight path angle. 
		#cos_fpa = self.h/(self.r*self.v)
		#self.phi = [np.arccos(cos_fpa[i]) if self.nu[i] >= np.pi else -np.arccos(cos_fpa[i]) for i in range(len(self.nu))]
		
		self.phi = np.zeros((self.r.size))
		self.phi[self.nu < np.pi] = np.arccos(self.h/(self.r[self.nu < np.pi]*self.v[self.nu < np.pi]))
		self.phi[self.nu >= np.pi] = -np.arccos(self.h/(self.r[self.nu >= np.pi]*self.v[self.nu >= np.pi]))
		
		#Input i, j and k vectors. 
		self.I = [1,0,0]
		self.J = [0,1,0]
		self.K = [0,0,1] 
		
		#Get initial x, y and z values. 
		self.x0 = self.r*np.cos(self.nu)
		self.y0 = self.r*np.sin(self.nu)
		self.z0 = self.r*np.zeros(self.x0.size)
		
		
		#Rotate ellipse around z axis by argument of periapsis. 
		self.x1 = self.x0*np.cos(self.omega) - self.y0*np.sin(self.omega)
		self.y1 = self.x0*np.sin(self.omega) + self.y0*np.cos(self.omega) 
		self.z1 = self.z0 
		
		
		#Rotate ellipse around x axis by inclination. 
		self.x2 = self.x1
		self.y2 = self.y1*np.cos(self.inc) - self.z1*np.sin(self.inc)
		self.z2 = self.y1*np.sin(self.inc) + self.z1*np.cos(self.inc) 
		
		#Rotate ellipse around z axis by RAAN. 
		self.x_gse = self.x2*np.cos(self.raan) - self.y2*np.sin(self.raan)
		self.y_gse = self.x2*np.sin(self.raan) + self.y2*np.cos(self.raan) 
		self.z_gse = self.z2 
		
		#Create radius vectors. 
		self.r_vector = np.array((self.x_gse, self.y_gse, self.z_gse)).T 
		
		#Get magnitude of new radius vectors. Should be same as r. 
		self.r_mag = np.sqrt(self.x_gse**2 + self.y_gse**2 + self.z_gse**2)
		
		#Get r unit vector. 
		self.r_unit_vector = np.array([self.r_vector[i]/self.r_mag[i] for i in range(self.r_mag.size)])
		
		#Get normal vector. 
		self.n_unit_vector = np.cross(self.r_unit_vector[0], self.r_unit_vector[90])
	
		#Get h vector. 
		self.h_vector = self.h*self.n_unit_vector 
		
		#Get tangential vector for each r vector. 
		self.tangent_unit = np.array([np.cross(self.n_unit_vector, self.r_unit_vector[i]) for i in range(self.r_mag.size)]) 
		
		#Get velocity unit vector. 
		self.v_unit_vector = np.array([np.cos(self.phi[i])*self.tangent_unit[i] + np.sin(self.phi[i])*self.r_unit_vector[i] for i in range(self.r_mag.size)]) 
		
		#Get velocity vector. 
		self.v_vector = np.array([self.v[i]*self.v_unit_vector[i] for i in range(self.r_mag.size)])
		
		#Check inclination angle. Correct.  
		cosi = np.dot(self.K, self.n_unit_vector)
		inclination = np.rad2deg(np.arccos(cosi))
		print ('Inclination = ', inclination) 
		
		#Get line of nodes vector and get RAAN. 
		self.nodes = np.cross(self.K, self.n_unit_vector) 
		self.nodes_mag = np.sqrt(self.nodes[0]**2 + self.nodes[1]**2 + self.nodes[2]**2) 
		self.nodes_unit = self.nodes/self.nodes_mag 
		print ('Nodes = ', self.nodes_unit)
		cossigma = np.dot(self.I, self.nodes_unit)
		sigma = np.rad2deg(np.arccos(cossigma))
		print ('RAAN = ', sigma) 
		
		#Check argument of perigee. 
		cosomega = np.dot(self.nodes_unit, self.r_unit_vector[0]) 
		print (cosomega)
		omega = np.rad2deg(np.arccos(cosomega))
		print ('Arg. Perigee = ', omega) 
		
		#Calculate the time since perigee. 
		self.calc_time() 
		
	def calc_time(self):
		'''This will calculate the time since periapsis for a given true anomaly.
		
		It will use a simple numerical method to perform the integration. Currently, 
		it just uses the spacing of the true anomaly array.'''
		
		#Split the main bracket into two terms for ease of calculation. 
		#self.term1 = (self.e*np.sin(self.nu))/((self.e**2-1)*(self.e*np.cos(self.nu)+1))	
		
		#self.arg_arctanh = ((self.e-1)*np.tan(self.nu/2))/((self.e**2 -1)**0.5)
		#self.term2 = (2*np.arctanh(self.arg_arctanh))/((self.e**2 - 1)**1.5)
		
		#Calculate p^2/h
		self.coeff = self.p**2/self.h
		
		self.t = np.zeros(self.nu.size) 
		
		#Calculate 1/(1+ecos(nu))^2
		self.nu_func = 1/(1 + self.e*np.cos(self.nu))**2
		
		for t in range(len(self.t)):
			if t == 0:
				print (self.nu_func[0])
				self.t[t] = 0
			
			else:
				self.integral = self.trapezium_rule(self.nu[1]-self.nu[0], self.nu_func[0:t+1])
				self.t[t] = self.coeff*self.integral 
		
		#Now calculate the time in seconds. 
		#self.t = self.coeff*self.integral 
	
	def trapezium_rule(self, nu_spacing, nu_func):
		'''This will integrate a function using the trapezium rule. '''
		
		return (nu_spacing/2)*(nu_func[0] + nu_func[-1] + 2*sum(nu_func[1:-1]))
		
	def get_datetime_from_periapsis(self):
		'''This will work out a datetime object for each point around the orbit by adding the time to the datetime object of the periapsis point.'''
    	
		self.dt_list = []
    	
		for t in range(len(self.t)):
			deltat = dt.timedelta(seconds=self.t[t]) 
			self.dt_list.append(self.ptime+deltat) 
    	
		print (self.dt_list)	
    
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
		#t = dt.datetime(self.dt_list)
		coord_obj.ticks = Ticktock(self.dt_list, 'UTC') 
    	
		#To convert to GSM. 
		self.coords_gsm = coord_obj.convert('GSM', 'car') 
		self.x_gsm = self.coords_gsm.x
		self.y_gsm = self.coords_gsm.y
		self.z_gsm = self.coords_gsm.z 
    	
	def plot_ellipse(self, lims=(-10,10), elev=45, azim=45):
		'''This plots the ellipse in 3D space.'''
		
		fig = plt.figure(figsize=(6,8))
		ax1 = fig.add_subplot(311, projection='3d')
		
		#Original ellipse
		#ax1.plot(self.x0/self.RE, self.y0/self.RE, self.z0/self.RE, 'b')
		
		#Add original periapsis vector. 
		#ax1.plot([0,self.x0[0]/self.RE], [0, self.y0[0]/self.RE], [0, self.z0[0]/self.RE], 'b')
		
		#Rotated by argument of periapsis. 
		#ax1.plot(self.x1/self.RE, self.y1/self.RE, self.z1/self.RE, 'r') 
		
		#Add periapsis vector after argument of periapsis. 
		#ax1.plot([0,self.x1[0]/self.RE], [0, self.y1[0]/self.RE], [0, self.z1[0]/self.RE], 'r')
		
		#Rotated by inclination. 
		#ax1.plot(self.x2/self.RE, self.y2/self.RE, self.z2/self.RE, 'g')
		
		#Add periapsis vector after inclination. 
		#ax1.plot([0,self.x2[0]/self.RE], [0, self.y2[0]/self.RE], [0, self.z2[0]/self.RE], 'g')
		
		#Rotated by RAAN. FINAL GSE. 
		zpos = np.where(self.z_gse >=0)
		zneg = np.where(self.z_gse < 0)
		ax1.plot(self.x_gse[zpos]/self.RE, self.y_gse[zpos]/self.RE, self.z_gse[zpos]/self.RE, 'k')
		ax1.plot(self.x_gse[zneg]/self.RE, self.y_gse[zneg]/self.RE, self.z_gse[zneg]/self.RE, 'gray')
		
		#Add periapsis vector after RAAN. 
		ax1.plot([0,self.x_gse[0]/self.RE], [0, self.y_gse[0]/self.RE], [0, self.z_gse[0]/self.RE], 'k')
		#Add other r vector. 
		#ax1.plot([0,self.x3[90]/self.RE], [0, self.y3[90]/self.RE], [0, self.z3[90]/self.RE], 'k')
		#Add normal unit vector. 
		ax1.plot([0,self.n_unit_vector[0]*5], [0, self.n_unit_vector[1]*5], [0, self.n_unit_vector[2]*5], 'b')
		
		#Add line of nodes. 
		ax1.plot([0, 5*self.nodes[0]], [0, 5*self.nodes[1]], [0, 5*self.nodes[2]], 'g')
		
		#Add GSM position to main graph. 
		ax1.plot(self.x_gsm/self.RE, self.y_gsm/self.RE, self.z_gsm/self.RE, 'b')
		
		#Add I vector. 
		#ax1.plot([0, 20*self.I[0]], [0, 20*self.I[1]], [0, 20*self.I[2]], 'r')
		
		#Add tangent vectors. 
		#for i in range(self.r_mag.size):
		#	if i%5==0:
		#		ax1.plot([self.x3[i]/self.RE, (self.x3[i]/self.RE+5*self.v_unit_vector[i][0])], [self.y3[i]/self.RE, (self.y3[i]/self.RE+5*self.v_unit_vector[i][1])], [self.z3[i]/self.RE, (self.z3[i]/self.RE+5*self.v_unit_vector[i][2])], 'r')
		
		ax1.set_xlabel('x')
		ax1.set_ylabel('y')
		ax1.set_zlabel('z')
		
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
		
		#Add another plot to show the variation of radial distance and velocity with true anomaly. 
		ax2 = fig.add_subplot(312)
		ax2.plot(np.rad2deg(self.nu), (self.r/self.RE), 'k')
		ax2.set_ylabel('Radius of Orbit (RE)') 
		ax2.set_xlabel('True Anomaly (deg)')
		
		ax3 = fig.add_subplot(313)
		ax3.plot(np.rad2deg(self.nu), self.t/3600, 'k')
		ax3.set_ylabel('Time (s)') 
		ax3.set_xlabel('True Anomaly (deg)')
		
		#Add period on. 
		ax3.plot([0,360], [self.period/3600, self.period/3600], 'k--', label="T={:.1f}hrs".format(self.period/3600))
		ax3.legend(loc='best')
		#ax3 = fig.add_subplot(414)
		#ax3.plot(np.rad2deg(self.nu), np.rad2deg(self.phi), 'k')
		#ax3.set_ylabel('Flight Path Angle (deg)') 
		#ax3.set_xlabel('True Anomaly (deg)')
		
		
		fig.savefig(self.plot_path+"example.png") 
		
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
		
		
