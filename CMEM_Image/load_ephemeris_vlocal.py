import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from datetime import datetime
import os
import pytz
#import re
from spacepy import coordinates as coord
from spacepy.time import Ticktock 
#import warnings
#warnings.filterwarnings('error')
from matplotlib.patches import Wedge, Polygon, Circle

UTC = pytz.utc
EARTH_RADIUS = 6371.009
RE = EARTH_RADIUS

class orbit():
    '''This class contains Yasir's code to read the orbit file and my code to read Yasir's code and plot what it reads!'''
    
    def __init__(self, stime=(2025,10,1), etime=(2025,10,3,3), calc_gsm=True):
        
        #Get plot path. 
        self.plot_path = os.environ.get('PLOT_PATH')+'smile_orbits/'
        
        self.stime = stime
        self.etime = etime
        self.calc_gsm = calc_gsm
        
        #Read in the data. 
        self.read_orbit_data(stime=self.stime, etime=self.etime, calc_gsm=self.calc_gsm)
        
        
    def format_date(self, t,rdate=""):
    
        t_dt = datetime.strptime(t, "%d %b %Y %H:%M:%S.%f")

        if rdate:
            replacement = datetime.strptime(rdate, "%Y-%m-%d")
            t_dt = t_dt.replace(year=replacement.year, month=replacement.month, day=replacement.day)

        date = t_dt.strftime("%Y%m%d")
        date = int(date)
        ut = t_dt.hour + t_dt.minute/60 + t_dt.second/3600

        return UTC.localize(t_dt), date, ut

    def format_re(self, x):
        return x/EARTH_RADIUS

    def load_ephemeris_vlocal(self, ephemeris_file, opt_skiprows=7, opt_date=None, opt_enddate=None):

        with open(ephemeris_file, 'r') as readfile:
            readlines = readfile.readlines()
            readlines = readlines[opt_skiprows:]
            time = []
            x = []
            y = []
            z = []
            x_der = []
            y_der = []
            z_der = []
    
            for l in readlines:
                line = l.split()
                line_time = line[0]+' '+line[1]+' '+line[2]+' '+line[3]
        
                time.append(line_time)
                x.append(float(line[4]))
                y.append(float(line[5]))
                z.append(float(line[6]))
                x_der.append(float(line[7]))
                y_der.append(float(line[8]))
                z_der.append(float(line[9]))

            time_len = len(time)
            date = range(0,time_len)
            ut = range(0,time_len)
            time = np.array(time, dtype='str')
            date = np.array(date, dtype='int')
            ut = np.array(ut, dtype='float')
            x = np.array(x, dtype='float')
            y = np.array(y, dtype='float')
            z = np.array(z, dtype='float')
    
            x = self.format_re(x)
            y = self.format_re(y)
            z = self.format_re(z)

            x_der = np.array(x_der, dtype='float')
            y_der = np.array(y_der, dtype='float')
            z_der = np.array(z_der, dtype='float')

            #print('time_len = ', time_len)
            for i in range(0, time_len):
                if opt_date:
                    rdate = opt_date
                    time[i], date[i], ut[i] = self.format_date(time[i], rdate)
                else:
                    time[i], date[i], ut[i] = self.format_date(time[i]) 
    
            readout = {'time' :time, 
               'date' : date,
               'ut' : ut,
               'x_gse' :x,
               'y_gse' :y,
               'z_gse' :z,
               'x_der':x_der,
               'y_der':y_der,
               'z_der':z_der
               }
    
        return readout
    
    def read_orbit_data(self, stime=(2025,10,1), etime=(2025,10,2), calc_gsm=False):
        '''This will read Yasir's function to just get the times I want.
    
        Parameters
        ----------
        stime - start time as a tuple (yyyy,mm,dd)
        etime - end time as a tuple (yyyy,mm,dd) 
        calc_gsm - boolean to convert to GSM. If after 2025, it will use the previous year's date. 
    
        Returns
        -------
        data - dictionary containing data in GSE coordinates. 
    
        '''    

        #Convert parameters to datetime objects. 
        stime = dt.datetime(*stime)
        etime = dt.datetime(*etime)

        #Read the whole file in. 
        orbit_file = '/data/sol-ionosphere/ys378/SMILE/smile_updated_ephemeris_gse_2025_2028.txt'
        data = self.load_ephemeris_vlocal(orbit_file) 
    
        #Convert time column to datetime objects. 
        data['dtime'] = np.array([dt.datetime.strptime(t, "%Y-%m-%d %H:%M:%S+00:0") for t in data['time']])
    
        #Now filter by time. 
        i = np.where((data['dtime'] >= stime) & (data['dtime'] < etime))
    
        new_data = {}
    
        #Create new dictionary with just the filtered data. 
        for k in data.keys():
            new_data[k] = data[k][i]
    
        if calc_gsm: 
        
            #Convert these coordinates to GSM. 
            #Need coordinates in appropriate form for conversion.
            print ('Convert coordinates to GSM...')
            coords_gse = np.zeros((new_data['x_gse'].size,3))
            coords_gse[:,0] = new_data['x_gse']
            coords_gse[:,1] = new_data['y_gse']
            coords_gse[:,2] = new_data['z_gse']
        
            #Needs to take data in a a list of xyz points. 
            coord_obj = coord.Coords(coords_gse, 'GSE', 'car')
        
            #Check how many years to subtract based on last time. 
            n = int(new_data['dtime'][-1].year-2024)
        
            #Add time information. 
            coord_obj.ticks = Ticktock(new_data['dtime']-dt.timedelta(days=365*n), 'UTC') 
        
            #To convert to GSM. 
            coords_gsm = coord_obj.convert('GSM', 'car') 
    
            new_data['x_gsm'] = coords_gsm.x
            new_data['y_gsm'] = coords_gsm.y
            new_data['z_gsm'] = coords_gsm.z 
        
            
        self.new_data = new_data  
    
    def plot_ellipse_2d(self, scatter=False, gse_col='darkblue', gsm_col='r'):
        '''This will plot the xy and xz planes of the orbit in GSE and GSM. 
        
        Parameters
        ----------
        scatter - Boolean to plot a scatter plot or a line plot. 
        '''
        
        self.dtime = self.new_data['dtime']
        self.x_gse = self.new_data['x_gse']
        self.y_gse = self.new_data['y_gse']
        self.z_gse = self.new_data['z_gse']
        self.x_gsm = self.new_data['x_gsm']
        self.y_gsm = self.new_data['y_gsm']
        self.z_gsm = self.new_data['z_gsm']
        
        
        #Set conditions for inbound and outbound. 
        #outward = np.where(self.phi > 0)
        #inward = np.where(self.phi < 0)
        
        fig = plt.figure(figsize=(6,6))
        fig.subplots_adjust(wspace=0.3)
        
        #GSE
        ax1 = fig.add_subplot(221)
        if scatter:
            ax1.scatter(self.x_gse, self.z_gse, c=gse_col, marker='x', s=5)
            #ax1.scatter(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, c='darkblue', marker='x', s=5)
        else:
            ax1.plot(self.x_gse, self.z_gse, gse_col)
            #ax1.plot(self.x_gse[inward]/self.RE, self.z_gse[inward]/self.RE, 'darkblue')
        
        ax1.set_xlabel(r'$x_{GSE}$')
        ax1.set_ylabel(r'$z_{GSE}$') 
        ax1.set_title('Orbit in GSE')
        self.make_earth_2d(ax1, rotation=-90)
        ax1.set_aspect('equal')
        
        ax2 = fig.add_subplot(223)
        if scatter:
            ax2.scatter(self.x_gse, self.y_gse, c=gse_col, marker='x', s=5)
            #ax2.scatter(self.x_gse, self.y_gse, c='darkblue', marker='x', s=5)
        else:
            ax2.plot(self.x_gse, self.y_gse, gse_col)
            #ax2.plot(self.x_gse, self.y_gse, 'darkblue')
        ax2.set_xlabel(r'$x_{GSE}$')
        ax2.set_ylabel(r'$y_{GSE}$') 
        #ax2.set_title('Orbit in GSE')
        self.make_earth_2d(ax2, rotation=-90)
        ax2.set_aspect('equal')
        
        ax3 = fig.add_subplot(222)
        if scatter:
            ax3.scatter(self.x_gsm, self.y_gsm, c=gsm_col, marker='x', s=5)
            #ax3.scatter(self.x_gsm, self.z_gsm, c='darkblue', marker='x', s=5)
        else:
            ax3.plot(self.x_gsm, self.z_gsm, gsm_col)
            #ax3.plot(self.x_gsm, self.z_gsm, 'darkblue')
        ax3.set_xlabel(r'$x_{GSM}$')
        ax3.set_ylabel(r'$z_{GSM}$') 
        ax3.set_title('Orbit in GSM')
        self.make_earth_2d(ax3, rotation=-90)
        ax3.set_aspect('equal')
        
        ax4 = fig.add_subplot(224)
        if scatter:
            ax4.scatter(self.x_gsm, self.y_gsm, c=gsm_col, marker='x', s=5)
            #ax4.scatter(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, c='darkblue', marker='x', s=5)
        else:
            ax4.plot(self.x_gsm, self.y_gsm, gsm_col)
            #ax4.plot(self.coords_gsm.x[inward]/self.RE, self.coords_gsm.y[inward]/self.RE, 'darkblue')
        ax4.set_xlabel(r'$x_{GSM}$')
        ax4.set_ylabel(r'$y_{GSM}$') 
        #ax4.set_title('Orbit in GSM')
        self.make_earth_2d(ax4, rotation=-90)
        ax4.set_aspect('equal')
        
        #fig.text(0.95, 0.05, 'Outbound', color='cyan', fontsize=8, ha='right')
        #fig.text(0.95, 0.03, 'Inbound', color='darkblue', fontsize=8, ha='right')
        #fig.text(0.05, 0.06, 'Inclination = {:.1f} deg'.format(np.rad2deg(self.inc)), fontsize=8, ha='left')
        #fig.text(0.05, 0.04, 'Arg. Periapsis = {:.1f} deg'.format(np.rad2deg(self.omega)), fontsize=8, ha='left')
        #fig.text(0.05, 0.02, 'RAAN* = {:.1f} deg'.format(np.rad2deg(self.raan)), fontsize=8, ha='left') 
        
        fig.text(0.5,0.95,'{} - {}'.format(self.dtime[0].strftime("%Y%m%d %H:%M"),self.dtime[-1].strftime("%Y%m%d %H:%M")), ha='center', fontsize=12)
        
        fig.savefig(self.plot_path+'orbit_plots_2d_{}_{}.png'.format(self.dtime[0].strftime("%Y%m%d %H:%M"),self.dtime[-1].strftime("%Y%m%d %H:%M")))
        
        
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
