import numpy as np
import datetime as dt
from datetime import datetime
import pytz
import re


UTC = pytz.utc
EARTH_RADIUS = 6371.009
RE = EARTH_RADIUS

def format_date(t,rdate=""):
    
    t_dt = datetime.strptime(t, "%d %b %Y %H:%M:%S.%f")

    if rdate:
       replacement = datetime.strptime(rdate, "%Y-%m-%d")
       t_dt = t_dt.replace(year=replacement.year, month=replacement.month, day=replacement.day)

    date = t_dt.strftime("%Y%m%d")
    date = int(date)
    ut = t_dt.hour + t_dt.minute/60 + t_dt.second/3600

    return UTC.localize(t_dt), date, ut

def format_re(x):
    return x/EARTH_RADIUS

def load_ephemeris_vlocal(ephemeris_file, opt_skiprows=7, opt_date=None, opt_enddate=None):

    readfile = open(ephemeris_file, 'r')
    readlines = readfile.readlines()
    readlines = readlines[opt_skiprows:]

    readlines_trunc = []

    time = []
    x = []
    y = []
    z = []
    x_der = []
    y_der = []
    z_der = []
    
    for line in readlines:
        readlines_trunc.append(line.strip())
        line = re.split("\s+", line.strip())
        line_time = line[0] + ' ' + line[1] + ' ' + line[2] + ' ' + line[3]
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

    x = format_re(x)
    y = format_re(y)
    z = format_re(z)

    x_der = np.array(x_der, dtype='float')
    y_der = np.array(y_der, dtype='float')
    z_der = np.array(z_der, dtype='float')

    #print('time_len = ', time_len)
    for i in range(0, time_len):
        if opt_date:
            rdate = opt_date
            time[i], date[i], ut[i] = format_date(time[i], rdate)
        else:
            time[i], date[i], ut[i] = format_date(time[i]) 
    
    readout = {'time' :time, 
               'date' : date,
               'ut' : ut,
               'x_re' :x,
               'y_re' :y,
               'z_re' :z,
               'x_der':x_der,
               'y_der':y_der,
               'z_der':z_der
               }
    
    readfile.close()
    
    return readout
    
def read_orbit_data(stime=(2025,10,1), etime=(2025,10,2)):
	'''This will read Yasir's function to just get the times I want.
	
	Parameters
	----------
	stime - start time as a tuple (yyyy,mm,dd)
	etime - end time as a tuple (yyyy,mm,dd) 
	
	Returns
	-------
	data - dictionary containing data in GSE coordinates. 
	
	'''    

	#Convert parameters to datetime objects. 
	stime = dt.datetime(*stime)
	etime = dt.datetime(*etime)

	#Read the whole file in. 
	orbit_file = '/data/sol-ionosphere/ys378/SMILE/smile_updated_ephemeris_gse_2025_2028.txt'
	data = load_ephemeris_vlocal(orbit_file) 
	
	#Convert time column to datetime objects. 
	data['dtime'] = np.array([dt.datetime.strptime(t, "%Y-%m-%d %H:%M:%S+00:0") for t in data['time']])
	
	#Now filter by time. 
	i = np.where((data['dtime'] >= stime) & (data['dtime'] < etime))
	
	new_data = {}
	
	#Create new dictionary with just the filtered data. 
	for k in data.keys():
		new_data[k] = data[k][i]
	
	return new_data 
	
