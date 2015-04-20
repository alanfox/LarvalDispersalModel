
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
import matplotlib.path as mplPath
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import shapefile
from bngtolatlon import OSGB36toWGS84

inputfilename = "C:/Users/af26/MSModelData/clim.txt"
MPA_SOURCE = 'Geikie Slide and Hebridean Slope'

NUM_LARVAE = 100

STARTDAY = 1        
SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = 3600.0

KM = np.array([1.0, 1.0, 0.0002]) #constant diffusion coefficient m2/s

VERTICAL_INTERP = True
ANIMATE = False
SETTLING = True
DEATH = False

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
    
def plot_animate():
    
    for larva in larvae_group:
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
        m.scatter(x1,y1,latlon = True, color = col)

    for larva in larvae_outofarea:
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
        m.scatter(x1,y1,latlon = True, color = col)

    m.etopo()
    m.drawcoastlines()
    
    #m.fillcontinents(color='coral',lake_color='skyblue')
    # draw parallels and meridians.
    m.drawparallels(np.arange(50.,66.,1.),labels = [1,1,0,0])
    m.drawmeridians(np.arange(-16.,17.,2.),labels = [0,0,0,1])
    m.drawmapboundary(fill_color='skyblue')
    plt.title("Larval dispersal. 'Realistic larval behaviour (Larsson et al.)")

    plot_file_name = '_temp%05d.png' % (day-STARTDAY)
    
    plt.savefig(plot_file_name)
    plt.clf()
    
def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
    
class Mpa:
    
    def __init__(self, shape, record, area_type):
        
        
        self.shape = shape
        self.record = record
        self.area_type = area_type
        self.bbox = shape.bbox
        self.bbox_points = [[self.bbox[0],self.bbox[1]],
                            [self.bbox[2],self.bbox[1]],
                            [self.bbox[2],self.bbox[3]],
                            [self.bbox[0],self.bbox[3]]]
        self.points = shape.points
# convert from bng to latlon
        if self.area_type == 'MPA' or self.area_type == 'MAR_SAC' or self.area_type == 'SPA':
            self.bbox = self.bng2lonlat_bbox(self.bbox)
            self.bbox_points = self.bng2lonlat(self.bbox_points)
            self.points = self.bng2lonlat(self.points)
            
        self.bbox_path = mplPath.Path(self.bbox_points)
        self.shape_path = mplPath.Path(self.points)
        self.nsettled = 0
        
    def get_bbox(self):
        return self.bbox
        
    def get_points(self):
        return self.points
        
    def get_sitename(self):
        if self.area_type == 'OFF_SAC':
            return self.record[1]
        if self.area_type == 'MPA':
            return self.record[0]
        if self.area_type == 'MAR_SAC':
            return self.record[0]
        if self.area_type == 'SPA':
            return self.record[0]
        
    def get_settled(self):
        return self.nsettled
        
    def settles(self,larva):
        # tests if an object of class Larva settles in the mpa
        if larva.ready_to_settle():
            pos = larva.get_position()
            x = pos[0]
            y = pos[1]
            if self.bbox_path.contains_point((x,y)):
                if self.shape_path.contains_point((x,y)):
                    self.nsettled = self.nsettled + 1
                    return True
        return False
        
    def bng2lonlat(self,bng):
        lonlat = []
        for i in range(len(bng)):
            x = OSGB36toWGS84(bng[i][0],bng[i][1])
            y = [x[1],x[0]]
            lonlat.append(y)
        return lonlat
        
    def bng2lonlat_bbox(self,bng):
        ll = OSGB36toWGS84(bng[0],bng[1])
        ur = OSGB36toWGS84(bng[2],bng[3])
        bbox = [ll[1],ll[0],ur[1],ur[0]]
        return bbox
        
    def plot_shape(self):
        x = []
        y = []
#        print self.shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
        m.plot(x,y, latlon = True, color = 'red')   
        
class Larva:
    
    # each larva will be an instance of this class

    def __init__(self, pos, vel, source):
        self.rundays = 0.0
        self.at_bed = False
        self.pos = np.array([pos[0], pos[1], pos[2]])
        #if pos[2] is negative, bed release
        
        if self.pos[2] < 0.0:
            self.bed_release(self.pos)
            
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
        
        
    def bed_release(self, position):
        
        # i coordinate nearest larva
        i = round((position[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((position[1] - minlat) / latstep)
          
        # determines k box with larva in
        zlevs = gridk[1]
        nz = len(zlevs)
        
        k = nz

        for layer in range(nz):
            if not mask[i,j,layer]:
                k = layer
                break
        
        self.pos[2] = zlevs[k - 1] - 5.0
        self.at_bed = True
        
        
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
            self.fu2 = interp1d(zlevs, u2, kind = 'linear')
            self.fv2 = interp1d(zlevs, v2, kind = 'linear')
            self.fw2 = interp1d(zlevs, w2, kind = 'linear')

    def advection(self):
        
        # presently no interpolation in the horizontal direction
        # interpolation in the vertical
        # this was to save me coding time, reduce run-time
        
        # update velocity
        # i coordinate nearest larva
        i = round((self.pos[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((self.pos[1] - minlat) / latstep)
        
        if VERTICAL_INTERP:
            # interpolate velocities in the vertical
            self.vertical_interpolation(i,j)
                  
            self.vel[0] = self.fu2(self.pos[2])
            self.vel[1] = self.fv2(self.pos[2])
            self.vel[2] = self.fw2(self.pos[2])
            
        else:    
            # or just go with box larva is in
            # k box with larva in
            zlevs = gridk[1]
            nz = len(zlevs)
            z = self.pos[2]
    
            for layer in range(nz):
                if z <= zlevs[layer]:
                    k = layer
                    break
            
            # interpolate velocities
            # don't bother for now, just use the gridbox the point is in
            
            self.vel[0] = u[i,j,k]
            self.vel[1] = v[i,j,k]
            self.vel[2] = w[i,j,k]            
        

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
        
        k = nz

        for layer in range(nz):
            if z <= zlevs[layer]:
                k = layer
                break
        if k == nz:
            return True
        else:
            return not mask[i,j,k]
        
        
    def ready_to_settle(self):
        return ((self.rundays > MINSETTLEAGE) and self.at_bed)
                        
                
    def update(self, rundays):
        
        self.rundays = rundays
        self.at_bed = False
        
        # updates the larva position, returns a boolean value - True if 
        # the larva has hit the bed or left the area, otherwise False
        
        # update position
        m_to_degree_lon = M_TO_DEGREE / np.cos(np.radians(self.pos[1]))
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
            self.at_bed = True
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
            self.at_bed = True
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
            self.at_bed = True
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

# helper functions
        
def group_settle(mpa_sprite_group, larva_object):
    settled = False
    for mpa in set(mpa_sprite_group):
        if mpa.settles(larva_object):
            settled = True
    return settled

def group_group_settle(larval_sprite_group, mpa_sprite_group):
    for larva in set(larval_sprite_group):
        if group_settle(mpa_sprite_group, larva):
            larval_sprite_group.remove(larva)
            settled_group.add(larva)            

np.random.seed(1)

# set up group of mpas

mpa_group = set([])

# offshore SAC
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
    
# SAC with marine components
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
    
# Nature conservation MPA
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
    
## SPAs
#shapes, records = read_shapefile('C:/Users/af26/Shapefiles/SPA_SCOTLAND_ESRI/SPA_SCOTLAND')
#for i in range(len(shapes)):
#    mpa_group.add(Mpa(shapes[i], records[i],'SPA'))
#    
# read in and calculate the model grid variables

gridlon, gridlat, gridk, minlon, maxlon, minlat, maxlat, lonstep, latstep = setupGrid(inputfilename) 
nx, ny, nz, infile = readHeader(inputfilename)

# read in the opening day's data

u, v, w = readHydrographicData(infile,STARTDAY, nx, ny, nz, first = True)
mask = landmask(u, v, w)

# initialise larvae. 
# Using grids of larvae at the same depth around a central point.

larvae_group = set([])
larvae_outofarea = set([])
settled_group = set([])

# seed larvae randomly in a particular mpa

for mpa in mpa_group:
    
    if mpa.get_sitename() == MPA_SOURCE:
        nlarvae = 0
        bbox = mpa.get_bbox()
        points = mpa.get_points()
        path = mplPath.Path(points)
        while nlarvae < NUM_LARVAE:
            x = -200.0
            y = -200.0
            while not path.contains_point((x,y)):
                x = np.random.uniform(bbox[0],bbox[2])
                y = np.random.uniform(bbox[1],bbox[3])
        # check not on land
            i = round((x - minlon) / lonstep)
            j = round((y - minlat) / latstep)
            if mask[i,j,0]:
                larvae_group.add(Larva([x, y, -1.0], [0.0,0.0,0.0],MPA_SOURCE))
                nlarvae = nlarvae + 1

# the main program loop
# initialise map

m = Basemap(projection='lcc',llcrnrlat=54.,llcrnrlon=-16.,urcrnrlat=62.,\
            urcrnrlon=2.,lat_1=50.,lon_0 = -7.0,resolution='h')
            
            
day = STARTDAY
nsteps = int(SECONDS_IN_DAY * 63.0 / DT)
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
#        print day
        u, v, w = readHydrographicData(infile,day, nx, ny, nz, first = False)

        if ANIMATE:
            plot_animate()
            
    if (SETTLING and (rundays > MINSETTLEAGE)):
        group_group_settle(larvae_group, mpa_group)
            
infile.close()

# plot the results

            
# Draw the larvae tracks on the map. End each track with a dot.
# Ugly if loop code allows different colour tracks for different sources.
           
for larva in larvae_group:
    x, y = larva.get_track()
    x1, y1, z1 = larva.get_position()
    col = '#444444'
    m.plot(x,y, latlon = True, color = col)
    m.scatter(x1,y1,latlon = True, color = col)
    
for larva in larvae_outofarea:
    x, y = larva.get_track()
    x1, y1, z1 = larva.get_position()
    col = '#000000'
    m.plot(x,y, latlon = True, color = col)
    m.scatter(x1,y1,latlon = True, color = col)
        
for larva in settled_group:
    x, y = larva.get_track()
    x1, y1, z1 = larva.get_position()
    col = '#880000'
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
plt.title("Larval dispersal 1 Jan release from bed. Larsson et al behaviour")
# draw the mpas

for mpa in mpa_group:
    mpa.plot_shape()
    print mpa.get_settled(), mpa.get_sitename()

#plt.savefig('foo.pdf')



# this draws a t,z plot of the position of the larvae in the water column

plt.figure()
plt.title("Larval dispersal 1 Jan release from bed. Larsson et al behaviour")
for larva in larvae_group:
    z = larva.get_depth_history()
    tstep = np.array(range(len(z)))
    t = tstep / 24.0
    plt.plot(t,z, color = '#444444')
for larva in larvae_outofarea:
    z = larva.get_depth_history()
    tstep = np.array(range(len(z)))
    t = tstep / 24.0
    plt.plot(t,z, color = '#000000')
for larva in settled_group:
    z = larva.get_depth_history()
    tstep = np.array(range(len(z)))
    t = tstep / 24.0
    plt.plot(t,z, color = '#880000')
    
plt.show()
