# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 09:31:17 2014

@author: af26
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import shapefile

STARTDAY = 182        
SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = 3600.0

KM = np.array([1.0, 1.0, 0.0002]) #constant diffusion coefficient m2/s

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
    
    for i in range(nx):
        gridlon[i] = minlon + (i * lonstep)
      
    for i in range(ny):
        gridlat[i] = minlat + (i * latstep)
        
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
    u = np.zeros((nx,ny,nz),dtype = float)
    v = np.zeros((nx,ny,nz),dtype = float)
    w = np.zeros((nx,ny,nz),dtype = float)
    
# scroll through input file to the required day on first call
    
    if first:
        data_lines_per_day = nx * ny * 5
    
        for iday in range(data_lines_per_day * (dayno - 1)):
            infile.readline() 
        
    
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
    return np.logical_or(u != 0.0,v != 0.0,w != 0.0)
    
def plot_shapefiles(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    for i in range(len(shapes)):
        x = []
        y = []
        points = shapes[i].points
        shape_path = mplPath.Path(points)
        print shape_path.contains_point((-14.0,58.0)), records[i][1]
        for point in points:
            x.append(point[0])
            y.append(point[1])
        m.plot(x,y, latlon = True, color = 'red')    
    
def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
    
class Mpa:
    
    def __init__(self, shape, record):
        
        self.shape = shape
        self.record = record
        self.bbox = shape.bbox
        self.points = shape.points
        
    def get_bbox(self):
        return self.bbox
        
    def get_points(self):
        return self.points
        
    def plot_shape(self):
        x = []
        y = []
        shape_path = mplPath.Path(self.points)
        print shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
        m.plot(x,y, latlon = True, color = 'red')    

class Larva:

    def __init__(self, pos, vel, source):
        self.pos = np.array([pos[0], pos[1], pos[2]])
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
                       
    def update(self):
        
        # update position
        m_to_degree_lon = M_TO_DEGREE / np.cos(np.radians(self.pos[0]))
        m_to_degree = np.array([m_to_degree_lon, M_TO_DEGREE, 1.0])
              
        # advection
        self.advection()       
        self.pos = self.pos + self.vel * m_to_degree * DT
        
        # diffusion
        self.diffusion()
        self.pos = self.pos + self.turb * m_to_degree

        # if leaves the surface, bounce back in!        
        self.pos[2] = self.pos[2] * np.sign(self.pos[2])
        
        # i coordinate nearest larva
        i = round((self.pos[0] - minlon) / lonstep)
        # j coordinate nearest larva
        j = round((self.pos[1] - minlat) / latstep)
        
        outofarea = False
        if i < 0 or i >= nx:
            outofarea = True
        if j < 0 or j >= ny:
            outofarea = True
            
        # k box with larva in
        zlevs = gridk[1]
        nz = len(zlevs)
        z = self.pos[2]

        for layer in range(nz):
            if z <= zlevs[layer]:
                k = layer
                break
        
        if not outofarea:
            landed = not mask[i,j,k]
        else:
            landed = True
                
        
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        
        return landed
        
    def diffusion(self):
        # 2-d horizontal diffusion with constant k. Gaussian randon walk.
        # no vertical diffusion yet
    
        rnorm = np.random.normal(0.0,1.0,3)
        self.turb = rnorm * np.sqrt(2.0 * KM * DT)


# set up group of mpas

mpa_group = set([])

shapes, records = read_shapefile('C:/Users/af26/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/UK_SAC_OFFSHORE_20121029')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i]))

np.random.seed(1)

inputfilename = "C:/Users/af26/MSModelData/clim.txt"
gridlon, gridlat, gridk, minlon, maxlon, minlat, maxlat, lonstep, latstep = setupGrid(inputfilename) 
nx, ny, nz, infile = readHeader(inputfilename)

larvae_group = set([])
larvae_onbed = set([])

# Rockall Bank 1
lon = -15.0
lat = 55.56
for i in range(5):
    for j in range(5):
        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
                                lat - 0.125 + j * 0.025, 5.0], 
                                [0.0,0.0,0.0],'rockall1'))

# Rockall Bank 2
#lon = -14.23
#lat = 57.84
#for i in range(10):
#    for j in range(10):
#        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
#                                lat - 0.125 + j * 0.025, 5.0], 
#                                [0.0,0.0,0.0],'rockall2'))
#
## Hebrides Terrace 
#lon = -10.3
#lat = 56.5
#for i in range(10):
#    for j in range(10):
#        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
#                                lat - 0.125 + j * 0.025, 5.0], 
#                                [0.0,0.0,0.0],'hebrides'))
#
## East Mingulay
#lon = -7.405
#lat = 56.78889
#for i in range(10):
#    for j in range(10):
#        larvae_group.add(Larva([lon - 0.25 + i * 0.05, 
#                                lat - 0.125 + j * 0.025, 5.0], 
#                                [0.0,0.0,0.0],'mingulay'))


nsteps = int(SECONDS_IN_DAY * 90.0 / DT)

m = Basemap(projection='lcc',llcrnrlat=54.,llcrnrlon=-16.,urcrnrlat=62.,\
            urcrnrlon=2.,lat_1=50.,lon_0 = -7.0,resolution='c')

u, v, w = readHydrographicData(infile,STARTDAY, nx, ny, nz, first = True)
mask = landmask(u, v, w)

day = STARTDAY
endday = False
runtime = 0.0
for t in range(nsteps):
    
    for larva in set(larvae_group):
        landed = larva.update()
        if landed:
            larvae_onbed.add(larva)
            larvae_group.remove(larva)
    runtime = (t + 1) * DT
    if ((runtime % SECONDS_IN_DAY) < DT/2.0):
        day = day + 1
        u, v, w = readHydrographicData(infile,day, nx, ny, nz, first = False)
            
infile.close()

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
    
for larva in larvae_onbed:
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
        
m.etopo()
m.drawcoastlines()

#m.fillcontinents(color='coral',lake_color='skyblue')
# draw parallels and meridians.
m.drawparallels(np.arange(50.,66.,1.),labels = [1,1,0,0])
m.drawmeridians(np.arange(-16.,17.,2.),labels = [0,0,0,1])
m.drawmapboundary(fill_color='skyblue')
plt.title("Larval dispersal July release at 5 m")


# draw the mpas

for mpa in mpa_group:
    mpa.plot_shape()

#plot_shapefiles('C:/Users/af26/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/UK_SACs_withMarineComponents_20130821')
#plot_shapefiles('C:/Users/af26/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/UK_SAC_OFFSHORE_20121029')
#plt.savefig('foo.pdf')
#plt.figure()
#
#tstep = np.array(range(nsteps+1))
#t = tstep / 24.0
#for larva in larvae_group:
#    z = larva.get_depth_history()
#    plt.plot(t,z)

plt.show()