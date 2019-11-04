#!/Users/you/bin/mypython
#----------------------------------------------------------------------------
# Name:         winsim7.py
# Purpose:      A fixed solar panel simulator returning instantaneous energies
#		in W/m^2 and energies over time in W*Hrs/m^2/time unit.
#		The user must be able to call "yearly_energy()" and "instant_energy(date_time)" 
#		from the interpretor.
#
#               Built for simulating Morgan Solar Window as an efficiency matrix of efficiency vs.
#               altitude and azimuth is implemented which can reflect the window's optical efficiencies.
#
#		Possibility exists for adapting simulator to include tracking
#		solar panels
#
# Author:       Peter Morgan
#
# Created:      2009, Toronto, Canada
#----------------------------------------------------------------------------

def helper():
        print """def instant_energy(self,date_time, accessed_from_interpretor = True)
def time_energy(self,start_time = 0.0,end_time = 365.0,list = False,output=None,step=1e-2)
winsim:
def __init__(self,lat=[43,40],win_vector=array([0.0,-1.0,0.0]),area = 1,CF = 1,
="perfect_eff.csv")
winsim_USA:
def __init__(self,city,win_vector=array([0.0,-1.0,0.0]),area = 1,eff_file="perfect_eff.csv",create_horz_inst=True)"""

from numpy import cos,sin,arcsin,arccos,pi,array,dot,cross,arange,size
from math import sqrt
from pylab import find,plot

import sys
sys.path.append("/Desktop/Window\ Simulator/Efficiency \Matrices/")

class winsim:

	def __init__(self,lat=39.833,N=array([0.0,-1.0,0.0]),B=array([0.0,0.0,-1.0]),area = 1,eff_file="perfect_eff.csv",eff_silicon = 0.19): 

                if type(lat) == float:  self.dec_lat = lat
                elif type(lat) == list: self.dec_lat = lat[0] + float(lat[1])/60	# dec_lat == decimal latitude (float)
		self.N = N; self.B = B; self.C = cross(N,B)
		self.mag_N = self.mag(N)
		self.area = area                                
		self.eff_matrix = self.read_array(eff_file)
		self.alt_array,self.azm_array,self.alt_increment,self.azm_increment \
		= self.create_altazm_arrays(self.eff_matrix)
		self.eff_silicon = 1.0

	def instant_energy(self, date_time, accessed_from_interpretor = True,is_horz_instance = False):
		alt,azm = self.calc_altazm(date_time)                           # alt,azm == altitude,azimuth
		if alt < 0.0:
			return 0.0
		else:                                           
                        S = self.calc_sun_vector(alt,azm)                       # S == solar position vector.
                        print S
                        win_alt,win_azm = self.calc_win_altazm(S)               # win_alt,win_azm == altitude & azimuth
			airmass = self.calc_airmass(alt)                        # from perspective of window orientation.
			fluence = self.calc_fluence(S)
			opt_eff = self.calc_efficiency(win_alt,win_azm)
                        energy = self.area*self.eff_silicon*opt_eff*fluence*1380*0.7**(airmass**0.678)
			if accessed_from_interpretor == False:
				return energy
			elif accessed_from_interpretor:
				print "%g W/%g m^2"%(energy,self.area)
                

	def time_energy(self,start_time = 0.0,end_time = 365.0,list=False,step=1e-2,eff_file=None):
                if eff_file != None:
                        self.change_eff_file(eff_file)                    
		h = step
		f = self.instant_energy		        # 	x == date_time. 
		x = start_time				# ie: 	x = 1.5 <--> Jan 1st, 12:00
		energy = 0.0				#	x = 365.0 <--> Dec 31st, 00:00
		if list == False:
                        while x <= end_time:
                                energy += h * 24 * ((f(x+h,False)+f(x,False))/2.0)
                                x += h
                        return "%g W*hrs/%g m^2/year" % (energy,self.area)
		elif list:
			return array([[time,self.instant_energy(time,False)] for time in arange(start_time,end_time,h)])

#------------------------------------------------------------------------------

	def create_altazm_arrays(self,matrix):	# This assumes that the altitude is on the
		alt_size = size(matrix,1)	# X-axis of the matrix.
		azm_size = size(matrix,0)
		alt_array = arange(-(pi/2.0),(pi/2.0)+(pi/2.0)/(alt_size-1),pi/(alt_size-1))
		azm_array = arange(-(pi/2.0),(pi/2.0)+(pi/2.0)/(azm_size-1),pi/(azm_size-1))
		return alt_array,azm_array,pi/(alt_size-1),pi/(azm_size-1)
		
	def read_array(self,filename, separator = ','):					# Returns a matrix from csv file
		print filename
		myfile = open(filename,'r')
		row_num = 0
		while myfile.readline() != '':
			row_num += 1
		myfile.seek(0)
		matrix = []
		for row in range(row_num):
			line = myfile.readline()
			fields = line.strip().split(separator)
			matrix.append(fields)
		for i in range(row_num):
			for j in range(len(fields)):
				matrix[i][j] = float(matrix[i][j])
		return array(matrix)

        def change_eff_file(self,eff_file):
                self.eff_matrix = self.read_array(eff_file)
		self.alt_array,self.azm_array,self.alt_increment,self.azm_increment \
		= self.create_altazm_arrays(self.eff_matrix)

        def change_eff_silicon(self,eff):
                self.eff_silicon = eff

#------------------------------------------------------------------------------

	def calc_altazm(self,date_time):			# positive azimuth == sun is in the East (morning)
		date_time = float(date_time)			# negative azimuth == sun is in the West (afternoon)
		lat = self.DegreesToRads(self.dec_lat)
		SDY = int(date_time)								# SDY == days since Jan.1st
		d = self.DegreesToRads(-23.45)*cos(((2*pi)/365)*(SDY+10))               # checked
		h = self.HoursToRads((date_time - SDY)*24 - 12.0)			
		alt = arcsin( cos(h)*cos(d)*cos(lat) + sin(d)*sin(lat) )
		azm = arcsin( sin(h)*cos(d)/cos(alt) )
		return alt,azm					

        def calc_win_altazm(self,S):
                proj_NC_S = S - self.proj(self.B,S)
                win_azm = arccos(dot(self.N,proj_NC_S)/(self.mag(self.N)*self.mag(proj_NC_S)))
                win_alt = arccos(dot(S,proj_NC_S)/(self.mag(S)*self.mag(proj_NC_S)))
  
                if (self.proj(self.B,S)/self.B)[0] > 0 or\
                   (self.proj(self.B,S)/self.B)[1] > 0 or\
                   (self.proj(self.B,S)/self.B)[2] > 0:
                    win_alt = - win_alt
                if (self.proj(self.C,S)/self.C)[0] < 0 or\
                   (self.proj(self.C,S)/self.C)[1] < 0 or\
                   (self.proj(self.C,S)/self.C)[2] < 0:
                    win_azm = - win_azm
                return win_alt,win_azm
		
	def calc_efficiency(self,win_alt,win_azm):
                f = self.RadsToDegrees
                alta = self.alt_array
                azma = self.azm_array
                if win_azm > (pi/2.0) or win_azm < (-pi/2.0):
                        return 0.0
		near_alt = size(find(self.alt_array < win_alt))-1		# near_alt == index of discrete point		#checked
		near_azm = size(find(self.azm_array < win_azm))-1		# in the array alt_array nearest win_alt	#checked
##                print "alt =      %-10g%-10g\nnear_alt = %-10g%-10g\nazm =      %-10g%-10g\nnear_azm = %-10g%-10g"\
##                      %(win_alt,f(win_alt),alta[near_alt],f(alta[near_alt]),
##                        win_azm,f(win_azm),azma[near_azm],f(azma[near_azm]))
##
		d_alt = (self.eff_matrix[near_azm][near_alt+1]	                # checked													
				- self.eff_matrix[near_azm][near_alt]) / (self.alt_increment)
		d_azm = (self.eff_matrix[near_azm+1][near_alt]
				- self.eff_matrix[near_azm][near_alt]) / (self.azm_increment)
		efficiency = (self.eff_matrix[near_azm][near_alt] + d_alt* (win_alt - 
						self.alt_array[near_alt]) + d_azm*(win_azm - self.azm_array[near_azm]))
		return efficiency
		
	def calc_airmass(self,alt):
		z = pi/2 - alt										# z == zenith angle
		a = 1.002432; 	b = 0.148386;	c = 0.0096467		# Constants used in
		d = 0.149864;	e = 0.0102963;	f = 0.000303978		# calculating airmass.  All checked.
		return (a*cos(z)**2 + b*cos(z) + c) / (cos(z)**3 + d*cos(z)**2 + e*cos(z) + f)	 	

	def calc_fluence(self,S):
		return dot(S,self.N)/(self.mag(S)*self.mag_N) 
		
	def calc_sun_vector(self,alt,azm):
		azm -= pi/2.0
		alt = pi/2.0 - alt									# Here alt == zenith angle
		return array([cos(azm)*sin(alt),sin(azm)*sin(alt),cos(alt)])
		
		# Using polar coordinates assumes that 0 is at x axis and grows positive towards
		# y-axis s.th. y-axis is at pi/2.
		# The convention we have adopted for this simulator is that the y-axis points north,
		# x-axis points north.  0 rads azimuth is due south (negative y axis), and grows
		# positive towards x-axis.
		# If we were to use the existing azm in these calculations, they would assume that 0
		# is at the x axis.  We have subtracted pi/2 to place 0 at the negative y-axis.
		
#------------------------------------------------------------------------------

        def mag(self,vector):
            return sqrt(dot(vector,vector))
        
        def proj(self,V,S):
                return (dot(S,V)/dot(V,V))*V
                
	def DegreesToRads(self,angle):
		return angle*pi/180.0

	def RadsToDegrees(self,angle):
                return angle*180.0/pi

	def HoursToRads(self,hourangle):
		return hourangle*pi/12

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class winsim_USA(winsim):
        
#-------- Overloaded Functions from Winsim ------------------------------------
        
	def __init__(self,city="Denver",N=array([0.0,-1.0,0.0]),B=array([0.0,0.0,-1.0]),
                     area = 1,eff_file="perfect_eff.csv",eff_silicon = 0.19,
                     create_horz_inst=True): 

                city_dict = {"Denver":("Denver_TMY3.csv",39.833,104.650,105),"Denver-Cent":("Denver2_TMY3.csv",39.742,105.179,105)}
                self.dec_lat = city_dict[city][1];
                self.dec_long = city_dict[city][2];
                long_meridian = city_dict[city][3]

##		self.dec_lat = lat[0] + float(lat[1])/60.0			# dec_lat == decimal latitude (float)
##		self.dec_long = long[0] + float(long[1])/60.0
##		self.dec_long_meridian = long[0] + float(long[1])/60.0

                self.CF_hist = []

		self.N = N; self.B = B
		self.C = cross(N,B)
		self.mag_N = self.mag(N)		
		self.area = area                                
		self.eff_matrix = self.read_array(eff_file)
		self.alt_array,self.azm_array,self.alt_increment,self.azm_increment \
		= self.create_altazm_arrays(self.eff_matrix)
                self.eff_silicon = eff_silicon
        
                self.month_dict = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
    		self.TMY3array = self.read_TMY3(city_dict[city][0])

    		if create_horz_inst:
                    self.horz_instance = winsim_USA(city,array([0.0,0.0,1.0]),array([0.0,1.0,0.0]),1,"perfect_eff.csv",0.19,False)

	def instant_energy(self, date_time, accessed_from_interpretor = True,is_horz_instance = False):
		alt,azm = self.calc_altazm(date_time)
		if alt < 0.0:

                        self.CF_hist.append([date_time,0.0])

			return 0.0

		
		else:                                   # CF == Instantaneous Correction Factor
                        S = self.calc_sun_vector(alt,azm)
                        win_alt,win_azm = self.calc_win_altazm(S)
			airmass = self.calc_airmass(alt)		# Checked!
			fluence = self.calc_fluence(S)
			opt_eff = self.calc_efficiency(win_alt,win_azm)
                        ideal_energy = self.area*fluence*1380*0.7**(airmass**0.678)
			if is_horz_instance == False:
                                CF = float(self.calc_TMY3_energy(date_time))/float(self.horz_instance.instant_energy(date_time,False,True))
                                if CF > 2.0:
                                      CF = 2.0
                                self.CF_hist.append([date_time,CF])

                                energy = CF * opt_eff * self.eff_silicon * ideal_energy
                                if accessed_from_interpretor == False:
                                        return energy
                                else:
                                        print "%g W/%g m^2"%(energy,self.area)       
			else:
				return ideal_energy

#-------- New Functions only existing in Winsim_USA ---------------------------

        def calc_TMY3_energy(self,date_time):           # TMY3 files have energy for every hour.
                xi = date_time                          # 'calc_TMY3_energy()' interpolates this data
                                                        # for any given point in the year.
                x1_index = size(find(self.TMY3array[:,0]<=date_time))-1
                if self.TMY3array[x1_index,0] == date_time:
                        return self.TMY3array[x1_index,1]
                else:
                        x2_index = x1_index + 1
                        x1,x2 = float(self.TMY3array[x1_index,0]),float(self.TMY3array[x2_index,0])
                        a = (x2-xi)/(x2-x1); b = (xi-x1)/(x2-x1)
                        interpolated_energy = a*float(self.TMY3array[x1_index,1])\
                        + b*float(self.TMY3array[x2_index,1])
                        return interpolated_energy
        
	def Date_to_DecTime(self,date,time,separator='/'):      # This function reads the time as
                dec_time=0.0                                    # recorded in TMY3 files.
                mmddyy = date.strip().split(separator)
                hour = float(time[:2])
                for month in arange(1.0,float(mmddyy[0])):
                        dec_time += self.month_dict[month]
                dec_time = float(dec_time)+float(mmddyy[1])-1.0+hour/24.0

                # "-1.0" is subtracted from the 'dec_time' float as the first day
                # in this simulator is 0.0, and thus noon on Jan.1 == 0.5 as opposed
                # to 1.5
                
                return dec_time
        
##                return self.local_to_solar_time(dec_time,self.dec_long,self.dec_long_meridian)
##
##        def local_to_solar_time(self,date_time,long,long_meridian):
##            N = int(date_time)
##            x = ((360*(N-1))/365.242)*(pi/180)
##            LST = (date_time - N) * 24
##            LC = (float(long) - float(long_meridian))/15.0
##            ET = 0.258*cos(x) - 7.416*sin(x) - 3.648*cos(2*x) - 9.228* sin(2*x)
##
##            ST = LST + (ET/60.0) - LC
##            return N+ST/24

        def read_TMY3(self,filename,output=None,separator = ','):

            # returns an array with the decimal date_time in column 1, and the GHI
            # (Global Horizontal Irradiance) in column 2.

            TMY3array = []
            file = open(filename,'r')
            row_num = 0
            while file.readline() != '':
                row_num += 1
                
            file.seek(0); file.readline(); file.readline()      # Jump to 3rd Line
            for row in range(1,row_num - 1):
                line = file.readline()
                fields = line.strip().split(separator)
                TMY3array.append((self.Date_to_DecTime(fields[0],fields[1])-(1.0/48.0),float(fields[4])))

            if output == None:
                return array(TMY3array)
            else:
                newfile = open(output,'w')
                for x in range(row_num - 2):
                        newfile.writeline(str(TMY3array[x])+'\n')
                newfile.close()
                print "Saved in %s" % (output)

#------------------------------------------------------------------------------

##                                END OF MODULE
                














