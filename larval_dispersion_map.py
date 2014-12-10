
"""
Created on Tue Nov 25 09:31:17 2014

@author: af26

Uses velocity fields from the model of Gallego et al 2013 (or Logemann 20??)
to advect and diffuse larvae.

Written in Python 2.7 (though I don't think there is anything here which isn't
Python 3) because I couldn't find the basemap package for Python 3 on Windows.

Only tested on Windows using Anaconda python installation.

The data file path will need to be modified. Data file is an ascii text file.
(It may be possible to run off the compressed file but I haven't figured out 
how yet) 

10/12/2014 Now includes larval behaviour. Starting at bed, rising through water 
column, then sinking. Unfortunately the model fields at depth don't appear to 
be up to the job. Need to implement running with POLCOMM model fields.


"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d

inputfilename = "C:/Users/af26/Documents/MSModelData/clim.txt"

STARTDAY = 182        
SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = 3600.0

KM = np.array([1.0, 1.0, 0.0002]) #constant diffusion coefficient m2/s

# constants for larval behaviour
# larvae released at bed, head upwards with increasing swim speeds up to 
# average age SWIMMAX. Then swimming gradually directed more downwards from 
# age DESCENDAGE up to DEADAGE.

SWIMSLOW = 0.0          #initial swimming speed
SWIMFAST = 3.0 / 1000.0 # max swimming speed
SWIMSTART = 0.0         #age in days at which start swimming
SWIMMAX = 14.0          #average age in days at which max swimming speed is reached
DESCENDAGE = 21.0       # average age at which probability of heading down starts
                        # to increase
FULLDESCENDAGE = 42.0    # now fully heading down
MINSETTLEAGE = 30.0     # minimum age at which can settle given suitable 
                        # habitat
DEADAGE = 56.0          # Average age at which dead

def setupGrid(filename):
    '''
    sets up the model output grid with information from file filename.
    But what does it return? Or do I need a grid class?
    '''
    
    infile = open(filename,"r")
    header = infile.readline()
    words = header.split()

    minlon = float(words[0])
    maxlon = float(words[1])
    minlat = float(words[2])
    maxlat = float(words[3])

    header = infile.readline()
    words = header.split()
    
    infile.close()
    
    nx = int(words[0])
    ny = int(words[1])
    nz = int(words[2])
    
    latstep = (maxlat - minlat) / float(ny - 1)
    lonstep = (maxlon - minlon) / float(nx - 1)
    kstep = [10.0,15.0,25.0,50.0,50.0,50.0,200.0,200.0,400.0,500.0,3500.0]

    gridlon = np.zeros((nx),dtype = float)
    gridlat = np.zeros((ny),dtype = float)
    gridk   = np.zeros((3,nz),dtype = float)
    
    # grid longitudes from model represent the centre of a grid box
    for i in range(nx):
        gridlon[i] = minlon + (i * lonstep)
      
    # grid latitudes again are the centre of a box
    for i in range(ny):
        gridlat[i] = minlat + (i * latstep)
        
    # input data file gives the top of each vertical grid box, heights of boxes
    # input in kstep above. gridk[0] stores the upper surface, grid[1] the lower 
    # surface, and gridk[2] the centre point
        
    gridk[0, 0] = 0.0
    gridk[1, 0] = kstep[0]
    gridk[2, 0] = kstep[0] / 2.0
    for i in range(1,nz):
        gridk[0, i] = gridk[0, i-1] + kstep[i-1]
        gridk[1, i] = gridk[0, i] + kstep[i]
        gridk[2, i] = (gridk[1, i] + gridk[0, i]) / 2.0
    
    return gridlon, gridlat, gridk, minlon, maxlon, minlat, maxlat, lonstep, latstep


def readHeader(filename):
    infile = open(filename,"r")
    header = infile.readline()
    header = infile.readline()
    words = header.split()
        
    nx = int(words[0])
    ny = int(words[1])
    nz = int(words[2])
    
    return nx, ny, nz, infile

def readHydrographicData(infile, dayno, nx, ny, nz, first):
    '''
    str -> 2 * array[i,j,k]
    Returns the 3-d  u and v velocity fields read from file 'filename'
    
    '''

    #file has no vertical velocity w, but I put it in for future use

    u = np.zeros((nx,ny,nz),dtype = float)
    v = np.zeros((nx,ny,nz),dtype = float)
    w = np.zeros((nx,ny,nz),dtype = float)
    
    # scroll through input file to the required day on first call
    
    if first:
        data_lines_per_day = nx * ny * 5
    
        for iday in range(data_lines_per_day * (dayno - 1)):
            infile.readline() 
        
    # read the data
    i = 0
    j = 0
    k = 0
    for i in range(nx):
        for j in range(ny):
            line = ''
            for m in range(5):
                line = line + ' ' + infile.readline()
            words = line.split()   
            for k in range(nz):
                u[i,j,k] = float(words[2 * k])
                v[i,j,k] = float(words[2 * k + 1])
#    infile.close()

    return u, v, w
    
def landmask(u,v,w):
    # mask is True when it is a water box, False when it is land.
    # these are not properly distinguished in the data file except by 
    # velocities being equal to zero.
    return np.logical_or(u != 0.0,v != 0.0,w != 0.0)
    
class Larva:
    
    # each larva will be an instance of this class

    def __init__(self, pos, vel, source):
        self.rundays = 0.0
        self.pos = np.array([pos[0], pos[1], pos[2]])
        self.newpos = np.array([pos[0], pos[1], pos[2]])
        self.vel = np.array([vel[0], vel[1], vel[2]])
        self.xlon_history = []
        self.xlat_history = []
        self.depth_history = []
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        self.turb = np.zeros((3),dtype = float)
        self.source = source
        self.i = 10000
        self.j = 10000
        self.fu2 = []
        self.fv2 = []
        self.fw2 = []
        
        #larval behaviour constants - vary slightly for each larva
        self.swimslow = SWIMSLOW
        self.swimfast = SWIMFAST 
        self.swimstart = SWIMSTART
        self.swimmax = SWIMMAX + np.random.normal(0.0,1.0)
        self.descendage = DESCENDAGE+ np.random.normal(0.0,1.0)
        self.fulldescendage = FULLDESCENDAGE+ 2 * np.random.normal(0.0,1.0)

# uncomment to activate settling and dying
#        self.minsettleage = MINSETTLEAGE
#        self.deadage = DEADAGE
        
        
    def get_position(self):
        return self.pos

    def get_source(self):
        return self.source
    
    def get_track(self):
        return self.xlon_history, self.xlat_history
        
    def get_depth_history(self):
        return self.depth_history        
    
    def get_velocity(self):
        return self.vel
    
    def vertical_interpolation(self,i,j):
        # extend arrays above surface and below bed for interpolation
        # only does this when the larva has moved horizontally to a new 
        # box. Otherwise uses stored values. To reduce recalculation.
    
        if ((i != self.i) or (j != self.j)):
            self.i = i
            self.j = j
            nz = len(gridk[2])
            zlevs = np.insert(np.append(gridk[2],6000.0),0,-5.0)
            u1 = u[i,j]
            v1 = v[i,j]
            w1 = w[i,j]
            u2 = np.insert(np.append(u1,u1[nz-1]),0,u1[0])
            v2 = np.insert(np.append(v1,v1[nz-1]),0,v1[0])
            w2 = np.insert(np.append(w1,w1[nz-1]),0,w1[0])
        
            # cubic spline interpolation 
            self.fu2 = interp1d(zlevs, u2, kind='cubic')
            self.fv2 = interp1d(zlevs, v2, kind='cubic')
            self.fw2 = interp1d(zlevs, w2, kind='cubic')

    def advection(self):
        
        # presently no interpolation in the horizontal direction
        # interpolation in the vertical
        # this was to save me coding time, reduce run-time
        
        # update velocity
        # i coordinate nearest larva
        i = round((self.pos[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((self.pos[1] - minlat) / latstep)
        
        # interpolate velocities in the vertical
        self.vertical_interpolation(i,j)
              
        self.vel[0] = self.fu2(self.pos[2])
        self.vel[1] = self.fv2(self.pos[2])
        self.vel[2] = self.fw2(self.pos[2])
        
    def diffusion(self):
        # 2-d horizontal diffusion with constant k. Gaussian randon walk.
        # vertical diffusion, Gaussian random walk, different constant.
    
        rnorm = np.random.normal(0.0,1.0,3)
        self.turb = rnorm * np.sqrt(2.0 * KM * DT)

    def vertical_behaviour(self):
        
        # implements swimming (or could be bouyancy) in the vertical.
        # days for each phase are taken from Larrson et al.
        # basically swims up at increasing rates to surface layer (5 m), 
        # stays there for a couple of weeks then heads back down to the bed
        
        # set percentage chance of heading upwards. If it doesn't go up it goes down.
        
        # set random swimming if in surface layer
        swim = np.zeros((3),dtype = float)                          
        if (self.pos[2] < 10.0):
            percent_up = 0.5
        else:
            if (self.rundays < self.swimstart):
                percent_up = 0.5
            elif (self.rundays < self.descendage):
                percent_up = 0.75
            elif (self.rundays > self.fulldescendage):
                percent_up = 0.25
            else:
                percent_up = 0.75 - 0.5 * ((self.rundays - self.descendage) 
                                      / (self.fulldescendage - self.descendage))
                                      
               
                
        if (self.rundays < self.swimstart):
            swimspeed = self.swimslow
        elif (self.rundays > self.swimmax):
            swimspeed = self.swimfast
        else:
            swimspeed = self.swimslow + ((self.swimfast - self.swimslow)
                                          * (self.rundays - self.swimstart)
                                          / (self.swimmax - self.swimstart))
                
#        print self.rundays, percent_up
        
        # set swim speed based on percent_up and swimspeed
        
        swimspeed = (swimspeed * 
                    float(np.random.choice([1,-1], 1, 
                                        p=[percent_up,1.0 - percent_up ])))
        swim[2] = swimspeed                                
                                        
        return swim              
                       
    def isoutofarea(self, position):

        # i coordinate nearest larva
        i = round((position[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((position[1] - minlat) / latstep)
        
        if i < 0 or i >= nx:
            return True
        if j < 0 or j >= ny:
            return True
            
        return False
        
    def isonland(self, position):
        
        # i coordinate nearest larva
        i = round((position[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((position[1] - minlat) / latstep)
          
        # determines k box with larva in
        zlevs = gridk[1]
        nz = len(zlevs)
        z = position[2]

        for layer in range(nz):
            if z <= zlevs[layer]:
                k = layer
                break
        
        return not mask[i,j,k]
                        
                
    def update(self, rundays):
        
        self.rundays = rundays
        
        # updates the larva position, returns a boolean value - True if 
        # the larva has hit the bed or left the area, otherwise False
        
        # update position
        m_to_degree_lon = M_TO_DEGREE / np.cos(np.radians(self.pos[0]))
        m_to_degree = np.array([m_to_degree_lon, M_TO_DEGREE, 1.0])
              
#        print self.pos

        # advection
        self.advection()       
        self.newpos = self.pos + self.vel * m_to_degree * DT
        
#        print 'advection', self.newpos
        
        # test if advected out of area, if so end
        
        if self.isoutofarea(self.newpos):
            self.pos = self.newpos
            self.xlon_history.append(self.pos[0])
            self.xlat_history.append(self.pos[1])
            self.depth_history.append(self.pos[2])
            return True
            
        # test if advected on to land
        # if it is don't advect
            
        if self.isonland(self.newpos):
            self.newpos = self.pos
#            print 'on land', self.newpos
        
        # lock in advection update
        self.pos = self.newpos
#        print 'after advection', self.pos

        # vertical swimming
            
        self.newpos = self.pos - DT * self.vertical_behaviour()
        
#        print 'swimming', self.newpos, self.pos
        
        # test if swims on to bed
        # if it is don't swim

        if self.isonland(self.newpos):
            self.newpos = self.pos
#            print 'on land', self.newpos 
        
        # test if swims out of surface
        # if it is don't swim

        if (self.newpos[2] < 0.0):
            self.newpos = self.pos
            
        # lock in swimming update
        self.pos = self.newpos
#        print 'after swimming', self.pos
            
                
        # diffusion
        self.diffusion()
        self.newpos = self.newpos + self.turb * m_to_degree
        
#        print 'diffusion', self.newpos
        # test if diffused out of area, if so end
        
        if self.isoutofarea(self.newpos):
            self.pos = self.newpos
            self.xlon_history.append(self.pos[0])
            self.xlat_history.append(self.pos[1])
            self.depth_history.append(self.pos[2])
            return True
            
        # test if diffused on to land
        # if it is don't diffuse
            
        if self.isonland(self.newpos):
            self.newpos = self.pos
#            print 'on land', self.newpos 
            
        # test if diffused out of surface
        # if it is, no vertical diffusion

        if (self.newpos[2] < 0.0):
            self.newpos[2] = self.pos[2]
                        
        # lock in diffusion update
        self.pos = self.newpos
#        print 'after diffusion', self.pos
                
        # store the new position in track history
        
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        
        return False
        


np.random.seed(1)



# read in and calculate the model grid variables

gridlon, gridlat, gridk, minlon, maxlon, minlat, maxlat, lonstep, latstep = setupGrid(inputfilename) 
nx, ny, nz, infile = readHeader(inputfilename)

# initialise larvae. 
# Using grids of larvae at the same depth around a central point.

larvae_group = set([])
larvae_outofarea = set([])

# Rockall Bank 1
lon = -15.0
lat = 55.56
for i in range(10):
    for j in range(10):
        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
                                lat - 0.125 + j * 0.025, 2000.0], 
                                [0.0,0.0,0.0],'rockall1'))
                                
# Rockall Bank 2
lon = -14.23
lat = 57.84
for i in range(10):
    for j in range(10):
        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
                                lat - 0.125 + j * 0.025, 300.0], 
                                [0.0,0.0,0.0],'rockall2'))

#Hebrides Terrace 
lon = -10.3
lat = 56.5
for i in range(10):
    for j in range(10):
        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
                                lat - 0.125 + j * 0.025, 2000.0], 
                                [0.0,0.0,0.0],'hebrides'))

# East Mingulay
lon = -7.405
lat = 56.78889
for i in range(10):
    for j in range(10):
        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
                                lat - 0.125 + j * 0.025, 75.0], 
                                [0.0,0.0,0.0],'mingulay'))
#

# read in the opening day's data

u, v, w = readHydrographicData(infile,STARTDAY, nx, ny, nz, first = True)
mask = landmask(u, v, w)

# the main program loop

day = STARTDAY
nsteps = int(SECONDS_IN_DAY * 56.0 / DT)
endday = False
runtime = 0.0
for t in range(nsteps):
    
    rundays = runtime / SECONDS_IN_DAY
    
#    print runtime, rundays
    
    for larva in set(larvae_group):
        left = larva.update(rundays)
        if left:
            larvae_outofarea.add(larva)
            larvae_group.remove(larva)
    runtime = (t + 1) * DT
    # read in a new current data field at the start of each day (daily mean fields)
    if ((runtime % SECONDS_IN_DAY) < DT/2.0):
        day = day + 1
        print day
        u, v, w = readHydrographicData(infile,day, nx, ny, nz, first = False)
            
infile.close()

# plot the results
# initialise map

m = Basemap(projection='lcc',llcrnrlat=54.,llcrnrlon=-16.,urcrnrlat=62.,\
            urcrnrlon=2.,lat_1=50.,lon_0 = -7.0,resolution='c')
            
# Draw the larvae tracks on the map. End each track with a dot.
# Ugly if loop code allows different colour tracks for different sources.
           
for larva in larvae_group:
    x, y = larva.get_track()
    x1, y1, z1 = larva.get_position()
    source = larva.get_source()
    if source == 'hebrides':
        col = '#444444'
    elif source == 'mingulay':
        col = '#444444'
    elif source == 'rockall1':
        col = '#444444'
    else: #source = 'rockall2'
        col = '#444444'
    m.plot(x,y, latlon = True, color = col)
    m.scatter(x1,y1,latlon = True, color = col)
    
for larva in larvae_outofarea:
    x, y = larva.get_track()
    x1, y1, z1 = larva.get_position()
    source = larva.get_source()
    if source == 'hebrides':
        col = '#000000'
    elif source == 'mingulay':
        col = '#000000'
    elif source == 'rockall1':
        col = '#000000'
    else: #source = 'rockall2'
        col = '#000000'
    m.plot(x,y, latlon = True, color = col)
    m.scatter(x1,y1,latlon = True, color = col)
        
# etopo is the topography/ bathymetry map background
m.etopo()
m.drawcoastlines()

# alternative plain map background
#m.fillcontinents(color='coral',lake_color='skyblue')

# draw parallels and meridians.
m.drawparallels(np.arange(50.,66.,1.),labels = [1,1,0,0])
m.drawmeridians(np.arange(-16.,17.,2.),labels = [0,0,0,1])
m.drawmapboundary(fill_color='skyblue')
plt.title("Larval dispersal July release at 5 m")
#plt.savefig('foo.pdf')



# this draws a t,z plot of the position of the larvae in the water column

plt.figure()
tstep = np.array(range(nsteps+1))
t = tstep / 24.0
for larva in larvae_group:
    z = larva.get_depth_history()
    plt.plot(t,z)

plt.show()