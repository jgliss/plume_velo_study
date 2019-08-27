"""
 Copyright 2016 Arve Kylling

 This file is part of the COMTESSA project software.

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this software.  If not, see <http://www.gnu.org/licenses/>.
 """
import copy
import math
import numpy as np
import pickle
import os
import shutil
from glob import glob
import multiprocessing
from subprocess import Popen,PIPE, STDOUT, call
import string, random
from netCDF4 import Dataset
import TomoCamera as TC
import time
import sys

def find_nearest_id(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx

def Average_spc_Files(InputFiles, OutputFile, verbose=False):
    # First check that all files have the same number of lines. If not
    # the files are surely different.
    i = 0
    #    for fn in glob(InputFiles):
    for fn in InputFiles:
        with open(fn) as fp:
            nlin = sum(1 for line in fp)
            if i==0:
                nlin0=nlin
            else:
                if nlin != nlin0:
                    print 'nlin: ' + str(nlin) + ', not equal nlin0: ' + str(nlin0)
                    exit(0)
            i = i + 1

    # All well? Combine all the files
    wvl = np.zeros([len(InputFiles),nlin])
    ix  = np.zeros([len(InputFiles),nlin],dtype=int)
    iy  = np.zeros([len(InputFiles),nlin],dtype=int)
    iz  = np.zeros([len(InputFiles),nlin],dtype=int)
    rad = np.zeros([len(InputFiles),nlin])
    s2  = np.zeros([nlin])
    radavg = np.zeros([nlin])
    i = 0
#    for f in  glob(InputFiles):
    for f in  InputFiles:
        (wvl[i],ix[i],iy[i],iz[i],rad[i]) = read_rad_spc(f, verbose=False)
        radavg[:] = radavg[:] + rad[i,:]
        s2[:]     = s2[:] + rad[i,:]*rad[i,:]
        i = i + 1

    s0 = i
    l = 0
    f = open(OutputFile,'w')
    while l < nlin:
        s1        = radavg[l]
        arg       = s0*s2[l] - s1*s1
        if arg < 0.0:
            print >> sys.stderr, l, arg, s0, s1, s2[l]
            arg = 0.0
        std       = (1.0/s0)*math.sqrt(arg)
        f.write('{0:8.2f} {1:3d} {2:3d} {3:3d} {4:9.4f} {5:9.4f}\n'.format(wvl[0,l], ix[0,l], iy[0,l], iz[0,l], s1/s0, std))
        l = l + 1
    f.close()
    return

def CombineSingleProcessOuput(InputFiles, OutputFile, verbose=False):
    fin = InputFiles[0]
    finstd = fin.replace('.spc','.std.spc')

    rad = np.loadtxt(fin)
    std = np.loadtxt(finstd)

    nwvl, ncol = rad.shape

    f = open(OutputFile,'w')

    iwvl=0
    while  iwvl < nwvl:
        f.write('{0:8.2f} {1:3d} {2:3d} {3:3d} {4:9.4f} {5:9.4f}\n'.format(rad[iwvl,0], int(rad[iwvl,1]),
                                                                           int(rad[iwvl,2]), int(rad[iwvl,3]),
                                                                           rad[iwvl,4], std[iwvl,4]))
        iwvl = iwvl + 1
    f.close()

    return

def read_rad_spc(fn, STD=False, verbose=False):
    # Read MYSTIC mc.rad.spc file
    if verbose:
        print "Reading MYSTIC mc.rad.spc file: ", fn
        sys.stdout.flush()
    if STD:
        wvl,ix,iy,iz,rad, std = np.loadtxt(fn, unpack=True)
        return (wvl,ix,iy,iz,rad,std)
    else:
        wvl,ix,iy,iz,rad = np.loadtxt(fn, unpack=True)
        return (wvl,ix,iy,iz,rad)


def Write2DMYSTICElevationFile(filename, nx, ny, dx, dy, elevation, verbose=False):
    if verbose:
        print "Write2DMYSTICElevationFile filename", filename
        sys.stdout.flush()

    fp=open(filename,'w')
    fp.write('{:d} {:d} {:f} {:f}\n'.format(nx, ny, dx, dy))
    ix=0
    while ix < nx:
        iy=0
        while iy < ny:
            fp.write('{:d} {:d} {:f}\n'.format(ix+1, iy+1, elevation ))
            iy=iy+1
        ix=ix+1

    fp.close()

def zenith(lat, lon, year, month, day,hour,min=0, sec=0, stdlong=0,output=1, uvspecpath=''):

    cmd = uvspecpath+'zenith '+str(day)+' '+str(month)+' '+str(hour)+' '+str(min)+' '+str(sec)+\
          ' '+'-y '+str(year)+' -a '+str(lat)+' -o '+str(lon)+' -s '+str(stdlong)+' -q'
    res = Popen(cmd,shell=True,stdout=PIPE)
    res.wait()
    vals = res.communicate()
    vals = vals[0].split()
    sza = float(vals[1])
    return sza

def Write3DMYSTICFile(filename, type='Generic', nx=0, ny=0, nz=0, dx=0, dy=0, z=None, LWCLimit=0.0, extLimit=0.0,
                      flag=0, verbose=False, **kwargs):
    if verbose:
        print "Write3DMYSTICFile filename", filename
        sys.stdout.flush()

    fp=open(filename,'w')
    if type=='Generic':
        fp.write('{:d} {:d} {:d} {:d}\n'.format(nx, ny, nz,flag))
        fp.write('{:f} {:f} '.format(dx, dy))
        for zz in z:
            fp.write('{:f} '.format(zz))
        fp.write('\n')
        if flag==1:
            tmpext  = kwargs['ext']
            tmpgg   = kwargs['gg']
            tmpssa  = kwargs['ssa']
            ix=0
            il=0
            while ix < nx:
                iy=0
                while iy < ny:
                    iz=0
                    while iz < nz:
                        #                        if tmpext[ix,iy,iz] > 0.0:
                        if tmpext[ix,iy,iz] > extLimit:
                            fp.write('{:d} {:d} {:d} {:f} {:f} {:f}\n'.format(ix+1, iy+1, iz+1, tmpext[ix,iy,iz], tmpgg[ix,iy,iz], tmpssa[ix,iy,iz] ))
                            il=il+1
                        iz=iz+1
                    iy=iy+1
                ix=ix+1
            # If no cloud  still include a dummy line with no optical thickness to fool MYSTIC
            # for background simulations.
            if il==0:
                ix=0
                iy=0
                iz=0
                fp.write('{:d} {:d} {:d} {:f} {:f} {:f}\n'.format(ix+1, iy+1, iz+1, 0.0, 0.0, 0.0 ))
        elif flag==3:
            tmpLWC  = kwargs['LWC']
            tmpreff = kwargs['reff']
            ix=0
            il=0
            while ix < nx:
                iy=0
                while iy < ny:
                    iz=0
                    #                    while iz < nz+1:
                    while iz < nz:
                        if tmpLWC[ix,iy,iz] > LWCLimit:
                            fp.write('{:d} {:d} {:d} {:g} {:f}\n'.format(ix+1, iy+1, iz+1, tmpLWC[ix,iy,iz], tmpreff[ix,iy,iz] ))
                            il=il+1
                        iz=iz+1
                    iy=iy+1
                ix=ix+1
            # If no cloud  still include a dummy line with no optical thickness to fool MYSTIC
            # for background simulations.
            if il==0:
                ix=0
                iy=0
                iz=0
                fp.write('{:d} {:d} {:d} {:f} {:f}\n'.format(ix+1, iy+1, iz+1, 0.0, tmpreff[ix,iy,iz] ))

    else:
        print 'Write3DMYSTICFile: Unknown type'
        exit
    fp.close()
    return

def make3DGaussian(sizex, sizey, sizez, sigmax=1, sigmay=1, sigmaz=1, center=None,verbose=False):
    """ Calculate 3D Gaussian distribution with different standard deviations
    in x, y, and z-directions.

    size? is the size in pixels/indices in ?-direction
    sigma ? is standard deviation in ?-direction

    """
    if verbose==True:
        print "make3DGaussian size", sizex, sizey, sizez
        print "make3DGaussian sigma", sigmax, sigmay, sigmaz
        sys.stdout.flush()
    x = np.arange(0, sizex, 1, float)

    y = x[:,np.newaxis]
    if sizey>sizex:
        y = np.insert(y,np.zeros(sizey-sizex),0,axis=0 )
    elif sizey<sizex:
        y = np.delete(y,np.arange(0,sizex-sizey),axis=0 )
    y[:,0] = np.arange(0, len(y), 1, float)

    z = y[:,np.newaxis]
    if sizez>sizey:
        z = np.insert(z,np.zeros(sizez-sizey),0,axis=0 )
    elif sizez<sizey:
        z = np.delete(z,np.arange(0,sizey-sizez),axis=0 )
    z[:,0,0] = np.arange(0, len(z), 1, float)

    if verbose==True:
        print z
        print "make3DGaussian x.shape", x.shape, y.shape, z.shape
        sys.stdout.flush()

    if center is None:
        x0 = sizex // 2
        y0 = sizey // 2
        z0 = sizez // 2
    else:
        x0 = center[0]
        y0 = center[1]
        z0 = center[2]

#    if verbose==True:
#        print "x", x
#        print "y", y
#        print "z", z
#        print "x0, y0, z0", x0, y0, z0

    tmpG = (1/(2.*sigmax*sigmay*sigmaz))*np.exp(-0.5*(((x-x0)**2)/sigmax**2+((y-y0)**2)/sigmay**2+((z-z0)**2)/sigmaz**2))
    # Transpose to comply with rest
    tmpG = tmpG.T

    if verbose==True:
        print "make3DGaussian dens", tmpG.shape
        sys.stdout.flush()

    return tmpG

def make3DEllipsoid(sizex, sizey, sizez, x_centre, y_centre, z_centre,
                    ellipsoid_a, ellipsoid_b, ellipsoid_c,
                    x, y, z, verbose=False):
    """ Set density to shape of a ellipsoid. Density is constant within.

    size? is the size in pixels/indices in ?-direction

    """
    tmpdens = np.zeros((sizex, sizey, sizez))
    ix=0
    for xx in x:
        iy=0
        for yy in y:
            iz=0
            for zz in z[1:z.shape[0]-1]:
                ell = (((xx-x_centre)**2.)/(ellipsoid_a**2.))+\
                      (((yy-y_centre)**2.)/(ellipsoid_b**2.))+\
                      (((zz-z_centre)**2.)/(ellipsoid_c**2.))
                if ell <=1.0:
                    tmpdens[ix,iy,iz]=1.0
                iz=iz+1
            iy=iy+1
        ix=ix+1

    if verbose==True:
        print "make3DEllipsoid dens", tmpdens.shape
        sys.stdout.flush()

    return tmpdens

def make3DVerticalPlume(sizex, sizey, sizez, iz_start, iz_end,
                        iBottomRadius, iTopRadius, sigmax, sigmay, scale_factor=1,
                        center=None, verbose=False):
    """ Calculate 3D Gaussian distribution with different standard deviations
    in x, y, and z-directions.

    size? is the size in pixels/indices in ?-direction
    sigma ? is standard deviation in ?-direction

    """
    if verbose==True:
        print "size", sizex, sizey, sizey
        print "iz_start", iz_start, iz_end, iBottomRadius, iTopRadius
        sys.stdout.flush()

    tmpdens = np.zeros((sizex, sizey, sizez))

    x = np.arange(0, sizex, 1, float)
    y = x[:,np.newaxis]
    if sizey>sizex:
        y = np.insert(y,np.zeros(sizey-sizex),0,axis=0 )
    elif sizey<sizex:
        y = np.delete(y,np.arange(0,sizex-sizey),axis=0 )
    y[:,0] = np.arange(0, len(y), 1, float)

    if center is None:
        x0 = sizex // 2
        y0 = sizey // 2
    else:
        x0 = center[0]
        y0 = center[1]

    if verbose==True:
        print "x", x
        print "y", y
        print "x0, y0, z0", x0, y0
        sys.stdout.flush()

    a = (iBottomRadius-iTopRadius)/(iz_start-iz_end)
    b = iBottomRadius-a*iz_start
    iz = iz_start
    ii =0
    while iz<iz_end:
        ir = a*iz+b
        this_sigmax = ir*sigmax
        this_sigmay = ir*sigmay
#        print iz, ir, a, b
        tmpG = (1/(2.*this_sigmax*this_sigmay))*np.exp(-0.5*(((x-x0)**2)/this_sigmax**2+((y-y0)**2)/this_sigmay**2))
        tmpdens[:,:,ii] = tmpG*scale_factor/tmpG.max()
#        print tmpG.shape, tmpG.max()
#        exit()
        ii=ii+1
        iz=iz+1


    if verbose==True:
        print x
        print x.shape, y.shape, z.shape
        sys.stdout.flush()

    return tmpdens

def make3DHorizontalPlume(sizex, sizey, sizez, iz_start, iz_end,
                        iBottomRadius, iTopRadius, sigmax, sigmay, scale_factor=1,
                        center=None, verbose=False):
    """ Calculate 3D Gaussian distribution with different standard deviations
    in x, y, and z-directions.

    size? is the size in pixels/indices in ?-direction
    sigma ? is standard deviation in ?-direction

    """
    if verbose==True:
        print "size", sizex, sizey, sizey
        print "iz_start", iz_start, iz_end, iBottomRadius, iTopRadius
        sys.stdout.flush()

    tmpdens = np.zeros((sizex, sizey, sizez))

    x = np.arange(0, sizex, 1, float)
    y = x[:,np.newaxis]
    if sizey>sizex:
        y = np.insert(y,np.zeros(sizey-sizex),0,axis=0 )
    elif sizey<sizex:
        y = np.delete(y,np.arange(0,sizex-sizey),axis=0 )
    y[:,0] = np.arange(0, len(y), 1, float)

    if center is None:
        x0 = sizex // 2
        y0 = sizey // 2
    else:
        x0 = center[0]
        y0 = center[1]

    if verbose==True:
        print "x", x
        print "y", y
        print "x0, y0, z0", x0, y0
        sys.stdout.flush()

    a = (iBottomRadius-iTopRadius)/(iz_start-iz_end)
    b = iBottomRadius-a*iz_start
    iz = iz_start
    ii =0
    while iz<iz_end:
        ir = a*iz+b
        this_sigmax = ir*sigmax
        this_sigmay = ir*sigmay
#        print iz, ir, a, b
        tmpG = (1/(2.*this_sigmax*this_sigmay))*np.exp(-0.5*(((x-x0)**2)/this_sigmax**2+((y-y0)**2)/this_sigmay**2))
        tmpdens[:,:,ii] = tmpG*scale_factor/tmpG.max()
#        print tmpG.shape, tmpG.max()
#        exit()
        ii=ii+1
        iz=iz+1


    if verbose==True:
        print x
        print x.shape, y.shape, z.shape
        sys.stdout.flush()

    return tmpdens

def make3DBox(sizex, sizey, sizez, verbose=False, scale_factor=1.0):
    """ Make homogeneous 3D 3x3x3 cube.
    """
    if verbose==True:
        print "make3DBox size", sizex, sizey, sizez
        sys.stdout.flush()

    tmpdens = np.ones((sizex, sizey, sizez))

    tmpdens = tmpdens*scale_factor

    if verbose==True:
        print "make3DBox", tmpdens.shape
        print "make3DBox", tmpdens
        sys.stdout.flush()

    return tmpdens

def make3DCell(sizex, sizey, sizez, verbose=False, scale_factor=1.0):
    """ Calculate 3D distribution for an SO2 cell. Here a 3x3x3 cube
    with a hole in the middle for the camera.
    """
    if verbose==True:
        print "make3DCell size", sizex, sizey, sizez
        sys.stdout.flush()

    tmpdens = np.ones((sizex, sizey, sizez))

    tmpdens = tmpdens*scale_factor
    x0 = sizex // 2
    y0 = sizey // 2
    z0 = sizez // 2
    tmpdens[x0,y0,z0] = 0.0

    if verbose==True:
        print "make3DCell", tmpdens.shape
        print "make3DCell", tmpdens
        sys.stdout.flush()

    return tmpdens

class ConversionFactors:
    def __init__(self):
        self.nm2cm = 1.e-07   # Convert from nm to cm
        self.nm2um = 1.e-03   # Convert from nm to cm
        self.um2m = 1.e-06   # Convert from um to m
        self.cm2km = 1.e-05   # Factor from cm to km
        self.m2km = 1.e-03   # Factor from m to km
        self.m2cm = 1.e+02   # Factor from m to cm
        self.km2m  = 1000.   # Factor from km to m
        self.km2cm  = 100000.   # Factor from km to m
        self.gm3togcm3=1e-06 # Convert from g/m**3 to g/cm**3
        self.kgtog=1.0e+03   # Convert from kg to g
        self.m3tocm3=1.0e+06 # Convert from m**3 to cm**3

class Camera:
    def __init__(self,RandString=''):
        self.Type='Camera'
        self.verbose=False
        self.name='TestCamera'
        self.savename='Camera'
        # Localisation of camera in m
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        # Viewing direction
        self.umu = 1.0 # Viewing the horizon
        self.phi = 0.0 # Looking north
        # Number of pixels in horizontal and vertical
        self.h_pixels=0
        self.v_pixels=0
        # Field of view of camera: phi1 phi2 theta1 theta2
        # phi = 0 is looking to the south
        # phi = 180 is with the sun in the back if phi0=0
        self.phi1 = 0.0
        self.phi2 = 0.0
        self.theta1 = 0.0
        self.theta2 = 0.0
        self.wavelength_grid_file='../Data/XSections/uvspec_SO2_wavelength_grid_file'
        self.wavelength=-9999
        self.wavelength1=-9999
        self.wavelength2=-9999
        if RandString == '':
            self.RandString= 'Cam'+''.join((random.sample(string.ascii_lowercase, 5)))
        else:
            self.RandString= 'Cam'+RandString
        return

    def info(self,latex=False):
        print
        print self.savename+' name:', self.name
        print '(All dimensions are in units of m)'
        print 'Localisation x={:f}, y={:f}, z={:f}'.format(self.x, self.y, self.z)
        print 'Pixels h_pixels={:d}, v_pixels={:f}'.format(self.h_pixels, self.v_pixels)
        print 'FOV phi1={:f}, phi2={:f}, theta1={:f}, theta2={:f}'.format(self.phi1, self.phi2, self.theta1, self.theta2)
        sys.stdout.flush()

        if latex:
            print '& {:5.1f} & {:5.1f} & {:5.1f} & {:5.1f} & {:6.1f} & {:6.1f} & {:6.1f} & {:6.1f} & {:d}& {:d}\\\\'.format(self.wavelength1, self.x, self.y, self.z, self.phi1, self.phi2, self.theta1, self.theta2, self.h_pixels, self.v_pixels)
            sys.stdout.flush()

        print
        return

    def save(self,folder):
        pkl_file = open(folder+self.savename+self.name+'.pkl', 'wb')
        pickle.dump(self,pkl_file )
        pkl_file.close()
        return

    def SetRTInput(self, UVS):
        print "Cam SetRTInput"
        sys.stdout.flush()

        UVS.inp['mc_sensorposition'] = '{:8.1f} {:8.1f} {:8.1f}'.format(self.x, self.y, self.z)
        UVS.inp['mc_sample_grid'] = '{:d} {:d}'.format(self.h_pixels, self.v_pixels)
        UVS.inp['mc_panorama_view'] = '{:f} {:f} {:f} {:f}'.format(self.phi1, self.phi2, self.theta1, self.theta2)
        UVS.inp['mc_panorama_alignment'] = 'mu'
        UVS.inp['umu'] = '{:f}'.format(self.umu)
        UVS.inp['phi'] = '{:f}'.format(self.phi)

#        UVS.inp['mc_panorama_alignment'] = 'sun'
#       UVS.inp['mc_panorama'] = 'weight_with_cos'
#        UVS.inp['mc_panorama'] = 'with_direct_rad'
#        UVS.inp['umu'] = '{:f}'.format((np.cos(np.deg2rad(0.5*(self.theta1+self.theta2)))))
#        UVS.inp['phi'] = '{:f}'.format(0.5*(self.phi1+self.phi2)-UVS.inp['phi0'])

        if self.wavelength1 != self.wavelength2:
            UVS.inp['wavelength_grid_file'] = self.wavelength_grid_file
        if not 'mol_abs_param' in UVS.inp:
            UVS.inp['mol_abs_param']    = 'crs'
        try:
            self.wavelength
            UVS.inp['wavelength']        = self.wavelength
        except:
            pass
        try:
            self.wavelength1
            self.wavelength2
            UVS.inp['wavelength']       = str(self.wavelength1)+' '+str(self.wavelength2)
        except:
            print "Both wavelength1 and wavelength2 must be given"
            exit()
        # try:
        #     self.filterfunction
        #     UVS.inp['filter_function_file'] = self.filterfunction
        #     UVS.inp['output_process'] = 'integrate'
        # except:
        #     pass

        return

class Spectrometer(Camera):
    def __init__(self, RunName=''):
        self.Type='Spectrometer'
        self.verbose=False
        self.name='Spectrometer'
        self.savename='Spectrometer'
        # Localisation of camera in m
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        # Number of pixels in horizontal and vertical
        self.h_pixels=1
        self.v_pixels=1
        # Field of view of camera: phi1 phi2 theta1 theta2
        # phi = 0 is looking to the south
        # phi = 180 is with the sun in the back if phi0=0
        self.phi1 = 0.0
        self.phi2 = 0.0
        self.theta1 = 0.0
        self.theta2 = 0.0
        self.mol_modify_o3=-9999
        self.crs_o3= '../Data/XSections/O3_Serdyuchenko_2014_223K_213-1100nm2013version.txt'
        self.crs_o4= '../Data/XSections/o4_thalman_volkamer_293K.dat'
        self.slitfunction = ''
        self.wavelength_grid_file='../Data/XSections/uvspec_SO2_wavelength_grid_file'
        self.RandString= 'Spec'+RunName+'_'+''.join(random.sample(string.ascii_lowercase, 5))

        return

    def SetRTInput(self, UVS):
        print self.savename+" SetRTInput"
        sys.stdout.flush()

        UVS.inp['mc_sensorposition'] = '{:8.1f} {:8.1f} {:8.1f}'.format(self.x, self.y, self.z)
        UVS.inp['mc_sample_grid'] = '{:d} {:d}'.format(self.h_pixels, self.v_pixels)
        UVS.inp['mc_panorama_view'] = '{:f} {:f} {:f} {:f}'.format(self.phi1, self.phi2, self.theta1, self.theta2)
        UVS.inp['wavelength']       = str(self.wavelength1)+' '+str(self.wavelength2)
        UVS.inp['mol_abs_param']    = 'crs'
        UVS.inp['crs_file O4']    =  self.crs_o4
        UVS.inp['crs_file O3']    =  self.crs_o3
        if self.mol_modify_o3>0.0:
            UVS.inp['mol_modify O3']    =  self.mol_modify_o3+' DU'

# Do this in a separate call to conv and spline after running uvspec
#        if self.slitfunction != '':
#            UVS.inp['slit_function_file'] = self.slitfunction
        UVS.inp['wavelength_grid_file'] = self.wavelength_grid_file

        return

    def CalculateColumnDensity(self):
        """
        NOTE:  It is assumed that integration is along x-axis for the
        center pixels in the y-direction.
        """
        for Plu in self.PlumeList:
            if 'SO2' in Plu.name:
                fact = 100.  # Convert from m to cm to column in cm-2
            else:
                fact = 1.0
                #            print "CalculateCol", Plu.name, fact
            # Calculate line integral for Spectrometer using Tomography software
            nx=Plu.nx
            dx=Plu.dx*1000 # #  In meters
            x0=0
            nz=Plu.nz
            z0=Plu.z[0]*1000.0
            dz=Plu.dz*1000.0 #100. #  In meters
            # RR=ReconstructionRegion
            RR =  TC.Area(nx=nx,dx=dx, x0=x0, nz=nz, dz=dz, z0=z0 )
            RR.zmin = RR.z0
            RR.zmax = RR.z0 + RR.dz*RR.nz
            RR.Image = np.zeros([nz,nx])
            #            print RR.Image.shape, Plu.dens.shape
            ycenter = int(Plu.dens.shape[1]/2) # Do this for center slice in y-direction
            RR.Image=Plu.dens[:,ycenter,1:] # Plu is one pixel larger in z-direction, inconsistent.
            RR.Image=RR.Image.T  # Sigh, why do I use different conventions for x and y......
            indices=np.argmax(RR.Image)
            #maxind=np.unravel_index(indices, RR.Image.shape)
            #print "RR", indices, maxind
            indices=np.argmax(Plu.dens)
            maxind=np.unravel_index(indices, Plu.dens.shape)
            #            print "Plu dens", indices, maxind
            #            print "RR.Image.max()", RR.Image.min(), RR.Image.max(), Plu.dens[maxind]

            theta1 = self.theta1-90  # And different conventions for the angles.....
            theta2 = self.theta2-90
            Nrays=9
            Camera1=TC.TomoCamera(x=-Plu.x_start*1000-Plu.x_length*1000/2, z=self.z, theta1=theta1, theta2=theta2, Nrays=Nrays,Name='Cam 1')
            Camera1.Rays()
            iRay=0
            sumtmpRq=0
            while iRay<Camera1.Nrays:
                tmpRq, tmpTotalLength, tmpind, tmpN =Camera1.CalculateLineIntegral(RR, iRay)
                Camera1.Sinogram[iRay]=tmpRq
                #                print '{:e}'.format(tmpRq)
                sumtmpRq=sumtmpRq+tmpRq
                iRay=iRay+1
            Plu.ColumnDensity = fact*sumtmpRq/Nrays
            #            fnTestLineIntegral='tmpgabba_'
            #            tmpfn = fnTestLineIntegral+'Cam1.dat'
            #            print 'tmpfn', tmpfn
            #            Camera1.WriteLineIntegralToFile(tmpfn)
        return

class Domain:
    def __init__(self):
        self.verbose=False
        self.name='TestDomain'
        # Domain size, all in km
        self.x_start = 0
        self.x_end   = 0.4
        self.dx      = 0.001
        self.nx=0
        self.x =None
        self.y_start = 0
        self.y_end   = 0.8
        self.dy      = 0.001
        self.ny=0
        self.y =None
        self.z_start = 0.150
        self.z_end   = 0.350
        self.dz      = 0.001
        self.nz=0
        self.z =None
        self.x_centre    = 0.5*(self.x_start+self.x_end)
        self.y_centre    = 0.5*(self.y_start+self.y_end)
        self.z_centre    = 0.5*(self.z_start+self.z_end)
        return

    def finalize(self):
        self.nx  = int(np.rint((self.x_end-self.x_start)/self.dx))
        self.x   = np.linspace(self.x_start,self.x_end,self.nx+1)
        self.x_size = self.x_end-self.x_start
        self.ny  = int(np.rint((self.y_end-self.y_start)/self.dy))
        self.y   = np.linspace(self.y_start,self.y_end,self.ny+1)
        self.y_size = self.y_end-self.y_start
        self.nz  = int(np.rint((self.z_end-self.z_start)/self.dz))
        self.z   = np.linspace(self.z_start,self.z_end,self.nz+1)
        self.z_size = self.z_end-self.z_start
        self.x_centre    = 0.5*(self.x_start+self.x_end)
        self.y_centre    = 0.5*(self.y_start+self.y_end)
        self.z_centre    = 0.5*(self.z_start+self.z_end)

        return

    def info(self):
        print
        print 'Domain name:', self.name
        print '(All dimensions are in units of km)'
        print 'x_start {:f}, x_end {:f}, dx {:f}, nx {:d}'.format(self.x_start, self.x_end, self.dx, self.nx)
        print 'y_start {:f}, y_end {:f}, dy {:f}, ny {:d}'.format(self.y_start, self.y_end, self.dy, self.ny)
        print 'z_start {:f}, z_end {:f}, dz {:f}, nz {:d}'.format(self.z_start, self.z_end, self.dz, self.nz)
        print 'x_centre {:f}, y_centre {:f}, z_centre {:f}'.format(self.x_centre, self.y_centre, self.z_centre)
        print
        sys.stdout.flush()

        return

    def save(self,folder):
        pkl_file = open(folder+'Domain.pkl', 'wb')
        pickle.dump(self,pkl_file )
        pkl_file.close()
        return



class Plume:
    def __init__(self):
        self.verbose=False
        self.name='TestPlume'
        self.randname = 'Plume'+''.join((random.sample(string.ascii_lowercase, 5)))
        self.shape= '' # Either set by input to set_density or inside ReadLESNetCDF
        self.x_start = None
        self.x_end   = None
        self.x_length = None
        self.y_start = None
        self.y_end   = None
        self.y_length = None
        self.z_start = None
        self.z_end   = None
        self.z_length = None
        self.set_density_flag=False
        self.ext = 0.0
        self.gg  = 0.0
        self.ssa = 0.0
        self.cdf = ''
        self.MYSTIC_profile_file_flag=0
        return

    def revise_domain(self, Domain,  verbose=False):
        if verbose:
            print "Inside revise_domain"
            sys.stdout.flush()

        # Only need to revise z-direction as x- and y-direction should be ok

        # For z-direction avoid first value as it has a half step.
        tmpdz = self.z[2]-self.z[1]
        tmpz_start=self.z[1]-tmpdz#*self.z.shape[0]
        # Add 1 as we want altitude at levels.
        tmpz = np.arange(tmpz_start, tmpz_start+tmpdz*float(self.z.shape[0]+1), tmpdz)

        Domain.z  = tmpz
        Domain.dz = tmpdz
        Domain.nz  = Domain.z.shape[0]-1 # -1 as this should by number of layers, not levels
        Domain.z_start=Domain.z[0]
        Domain.z_end=Domain.z[Domain.nz]
        Domain.z_size = Domain.z_end-Domain.z_end

        return Domain

    def finalize_size(self, Domain, verbose=False):
        if verbose:
            print "Inside finalize_size"
            sys.stdout.flush()

        self.x_end    = self.x_start+self.x_length
        self.y_end    = self.y_start+self.y_length
        self.z_end    = self.z_start+self.z_length
        self.x_centre    = self.x_start+0.5*self.x_length
        self.y_centre    = self.y_start+0.5*self.y_length
        self.z_centre    = self.z_start+0.5*self.z_length
        self.dx       = Domain.dx
        self.dy       = Domain.dy
        self.dz       = Domain.dz
        self.nx  = int(np.rint((self.x_end-self.x_start)/self.dx))
        self.x   = np.linspace(self.x_start,self.x_end,self.nx+1)
        self.ny  = int(np.rint((self.y_end-self.y_start)/self.dy))
        self.y   = np.linspace(self.y_start,self.y_end,self.ny+1)
        self.nz  = int(np.ceil((self.z_end-self.z_start)/self.dz))
        self.z   = np.linspace(self.z_start,self.z_end,self.nz+1)

        # Check that plume is entirely within domain
        if self.x_start < Domain.x_start or self.x_start > Domain.x_end:
            print 'Plume finalize: x_start {:f} outside domain {:f} {:f}'.format(self.x_start, Domain.x_start, Domain.x_end)
            exit()
        if self.x_end < Domain.x_start or self.x_end > Domain.x_end:
            print 'Plume finalize: x_end {:f} outside domain {:f} {:f}'.format(self.x_end, Domain.x_start, Domain.x_end)
            exit()
        if self.y_start < Domain.y_start or self.y_start > Domain.y_end:
            print 'Plume finalize: y_start {:f} outside domain {:f} {:f}'.format(self.y_start, Domain.y_start, Domain.y_end)
            exit()
        if self.y_end < Domain.y_start or self.y_end > Domain.y_end:
            print 'Plume finalize: y_end {:f} outside domain {:f} {:f}'.format(self.y_end, Domain.y_start, Domain.y_end)
            exit()
        if self.z_start < Domain.z_start or self.z_start > Domain.z_end:
            print 'Plume finalize: z_start {:f} outside domain {:f} {:f}'.format(self.z_start, Domain.z_start, Domain.z_end)
            exit()
        if self.z_end < Domain.z_start or self.z_end > Domain.z_end:
            print 'Plume finalize: z_end {:f} outside domain {:f} {:f}'.format(self.z_end, Domain.z_start, Domain.z_end)
            print self.z_end-Domain.z_start, self.z_end-Domain.z_end
            exit()

        # Get pixel indixes for plume within domain. These pixel indices correspond
        # to those needed by MYSTIC
        if self.shape=='User':
            #            self.ix_start = self.xadd*self.x.shape[0]
            #            self.ix_end   = (self.xadd+1)*self.x.shape[0]
            #            self.ix_start = 0
            #            self.ix_end   = self.x.shape[0]


            id = find_nearest_id(Domain.x,self.x_start)
            self.ix_start = id
            self.ix_end   = id+self.x.shape[0]
            #            self.iy_start = self.yadd*self.y.shape[0]
            #            self.iy_end   = (self.yadd+1)*self.y.shape[0]
            #            self.iy_start = 0
            #            self.iy_end   = self.y.shape[0]
            id = find_nearest_id(Domain.y,self.y_start)
            self.iy_start = id
            self.iy_end   = id+self.y.shape[0]
            id = find_nearest_id(Domain.z,self.z_start)
            self.iz_start = id
            self.iz_end   = id+self.z.shape[0]
#            self.iz_start = 0
#            self.iz_end   = self.z.shape[0]
#            self.iz_start = self.zadd*self.z.shape[0]
#            self.iz_end   = (self.zadd+1)*self.z.shape[0]
        else:
            self.ix_start = int((self.x_start-Domain.x_start)/Domain.dx)
            self.ix_end   = int((self.x_end-Domain.x_start)/Domain.dx+Domain.dx/2.)
            self.ix_size  = self.ix_end-self.ix_start
            self.iy_start = int((self.y_start-Domain.y_start)/Domain.dy)
            self.iy_end   = int((self.y_end-Domain.y_start)/Domain.dy+Domain.dy/2.)
            self.iy_size  = self.iy_end-self.iy_start
            self.iz_start = int((self.z_start-Domain.z_start)/Domain.dz)
            self.iz_end   = int((self.z_end-Domain.z_start)/Domain.dz+Domain.dz/2.)
            self.iz_size  = self.iz_end-self.iz_start

        return

    def info(self,latex=False, Domain=None):
        CF=ConversionFactors()
        print
        print 'Plume name:', self.name
        print '(All dimensions are in units of km)'
        print 'x_start {:f}, x_end {:f}'.format(self.x_start, self.x_end)
        print 'y_start {:f}, y_end {:f}'.format(self.y_start, self.y_end)
        print 'z_start {:f}, z_end {:f}'.format(self.z_start, self.z_end)
        if self.set_density_flag==True and 'Gauss' in self.shape:
            print 'x_sigma: {:f} y_sigma: {:f} z_sigma: {:f}'.format(self.x_sigma, self.y_sigma, self.z_sigma)
        if latex:
            print '& {:4.1f} & {:4.1f} & {:6.4f} & {:4.1f} & {:4.1f} & {:6.4f} & {:4.1f} & {:4.1f} & {:6.4f}& {:8.2f}\\\\'.format(
                self.x_start*CF.km2m, self.x_end*CF.km2m, self.x_sigma*CF.km2m,
                self.y_start*CF.km2m, self.y_end*CF.km2m, self.y_sigma*CF.km2m,
                self.z_start*CF.km2m, self.z_end*CF.km2m, self.z_sigma*CF.km2m, self.dens.max())

        print 'Indices'
        if Domain != None:
            print 'ix_start {:d}, ix_end {:d}, {:f}, {:f}'.format(self.ix_start, self.ix_end, Domain.x[self.ix_start], Domain.x[self.ix_end-1])
            print 'iy_start {:d}, iy_end {:d}, {:f}, {:f}'.format(self.iy_start, self.iy_end, Domain.y[self.iy_start], Domain.y[self.iy_end-1])
            print 'iz_start {:d}, iz_end {:d}, {:f}, {:f}'.format(self.iz_start, self.iz_end, Domain.z[self.iz_start], Domain.z[self.iz_end-1])
        else:
            print 'ix_start {:d}, ix_end {:d}'.format(self.ix_start, self.ix_end)
            print 'iy_start {:d}, iy_end {:d}'.format(self.iy_start, self.iy_end)
            print 'iz_start {:d}, iz_end {:d}'.format(self.iz_start, self.iz_end)
        print 'Dens'
        print 'dens.shape {:s} {:f} {:e}'.format(self.dens.shape, self.dens.max(), self.dens.max())
 #       if self.set_density_flag==True:
 #           print 'ix_sigma: {:d} iy_sigma: {:d} iz_sigma: {:d}'.format(self.ix_sigma, self.iy_sigma, self.iz_sigma)
        print
        sys.stdout.flush()

        return

    def set_density(self, shape, **kwargs):
        """
        Specify the shape and density of the plume. For each shape different
        input is needed. These are given in the kwargs dictionary. For the
        various shapes this is:

        shape=='3DGaussian':
        Required input:
            x_sigma, y_sigma, z_sigma: standard deviations (km) in each direction
            ext:  Extinction (km-1). The maximum value in the density field will
                  be set to this value, and the rest scaled accordingly.
        Optional input:
            The pklfile and pklfile_z_integrated inputs are useful for checking
            of plume shape and magnitude of density field and are used together
            with separate scripts for plotting.
            pklfile: dump pickle of density field to specified file
            pklfile_z_integrated: dump pickle of density field integrated allong
                 z-axis to specified file
        """
        self.set_density_flag=True
        self.shape = shape
        if kwargs=={}:
            print "Plume set_density: No density information provided"
            exit()
        if 'ext' in kwargs and 'gg' in kwargs and 'ssa' in kwargs:
            self.MYSTIC_profile_file_flag =1
        elif 'ext' in kwargs and 'reff' in kwargs:
            self.MYSTIC_profile_file_flag =2
        elif 'LWC' in kwargs and 'reff' in kwargs:
            self.MYSTIC_profile_file_flag =3

        if self.MYSTIC_profile_file_flag == 1:
            self.ext  = kwargs['ext']
            self.gg   = kwargs['gg']
            self.ssa  = kwargs['ssa']
        elif self.MYSTIC_profile_file_flag == 3:
            self.LWC  = kwargs['LWC']
            self.reff = kwargs['reff']

        if shape=='3DGaussian':
            self.x_sigma = kwargs['x_sigma']
            self.y_sigma = kwargs['y_sigma']
            self.z_sigma = kwargs['z_sigma']
            self.ix_sigma= int(self.x_sigma/self.dx)
            self.iy_sigma= int(self.y_sigma/self.dy)
            self.iz_sigma= int(self.z_sigma/self.dz)
            if self.ix_sigma<=0: self.ix_sigma=1
            if self.iy_sigma<=0: self.iy_sigma=1
            if self.iz_sigma<=0: self.iz_sigma=1
            # Gaussian distribution is made in pixel indices coordinates:
            self.dens=make3DGaussian(self.ix_size, self.iy_size, self.iz_size,
                                     sigmax=self.ix_sigma, sigmay=self.iy_sigma, sigmaz=self.iz_sigma, center=None)

            # Scale plume to wanted extinction/optical depth
            if self.MYSTIC_profile_file_flag == 1:
                self.dens = self.dens*self.ext/self.dens.max()
            elif self.MYSTIC_profile_file_flag == 3:
                if self.LWC<0.0:
                    self.dens = self.dens
                else:
                    self.dens = self.dens*self.LWC/self.dens.max()


        elif shape=='3DVerticalPlume':
            self.TopRadius  = kwargs['TopRadius']
            self.BottomRadius  = kwargs['BottomRadius']
            self.x_sigma = kwargs['x_sigma']
            self.y_sigma = kwargs['y_sigma']
            self.iTopRadius = self.TopRadius/self.dx  # TopRadius in pixels
            self.iBottomRadius = self.BottomRadius/self.dx  # BottomRadius in pixels
            if self.MYSTIC_profile_file_flag == 1:
                dummy = self.ext
            elif self.MYSTIC_profile_file_flag == 3:
                dummy = self.LWC
            self.dens=make3DVerticalPlume(self.ix_size, self.iy_size, self.iz_size,
                                          self.iz_start, self.iz_end, self.iBottomRadius, self.iTopRadius,
                                          self.x_sigma, self.y_sigma, scale_factor=dummy)

        elif shape=='Box':
#            print "self.nx", self.nx, self.ny, self.nz, self.x, self.z
#            print "self.nx", self.ix_size, self.ix_end, self.ix_start, \
#                self.iy_size, self.iy_end, self.iy_start, \
#                self.iz_size, self.iz_end, self.iz_start
            if self.MYSTIC_profile_file_flag == 1:
                dummy = self.ext
            elif self.MYSTIC_profile_file_flag == 3:
                dummy = self.LWC
            self.dens=make3DBox(self.ix_size, self.iy_size, self.iz_size,
                                 scale_factor=dummy)
        elif shape=='Cell':
#            print "self.nx", self.nx, self.ny, self.nz, self.x, self.z
#            print "self.nx", self.ix_size, self.ix_end, self.ix_start, \
#                self.iy_size, self.iy_end, self.iy_start, \
#                self.iz_size, self.iz_end, self.iz_start
            if self.MYSTIC_profile_file_flag == 1:
                dummy = self.ext
            elif self.MYSTIC_profile_file_flag == 3:
                dummy = self.LWC
            self.dens=make3DCell(self.ix_size, self.iy_size, self.iz_size,
                                 scale_factor=dummy)

        elif shape=='SingleLayer':
            print self.ix_size, self.iy_size, self.iz_size
            tmpdens = np.zeros((self.ix_size, self.iy_size, self.iz_size))
            print tmpdens.shape
            ix=0
            while  ix < self.ix_size:
                iy=0
                while  iy < self.iy_size:
                    iz=0
                    while  iz < self.iz_size:
                        tmpdens[ix,iy,iz]=self.LWC
                        iz=iz+1
                    iy=iy+1
                ix=ix+1
            self.dens = tmpdens

        elif shape=='Ellipsoid':
            self.x_centre = kwargs['x_centre']
            self.y_centre = kwargs['y_centre']
            self.z_centre = kwargs['z_centre']
            self.ellipsoid_a = kwargs['a']
            self.ellipsoid_b = kwargs['b']
            self.ellipsoid_c = kwargs['c']

            self.ix_centre= int((self.x_centre-self.x_start)/self.dx)
            self.iy_centre= int((self.y_centre-self.y_start)/self.dy)
            self.iz_centre= int((self.z_centre-self.z_start)/self.dz)

            # Ellipsoid shaped constant concentration is made in pixel indices coordinates:
            self.dens=make3DEllipsoid(self.ix_size, self.iy_size, self.iz_size,
                                      self.x_centre, self.y_centre, self.z_centre,
                                      self.ellipsoid_a, self.ellipsoid_b, self.ellipsoid_c,
                                      self.x, self.y, self.z)

            # Scale plume to wanted extinction/optical depth
            if self.MYSTIC_profile_file_flag == 1:
                self.dens = self.dens*self.ext/self.dens.max()
            elif self.MYSTIC_profile_file_flag == 3:
                self.dens = self.dens*self.LWC/self.dens.max()

        elif shape=='User':

            # Scale plume to wanted extinction/optical depth
            if self.MYSTIC_profile_file_flag == 1:
                self.dens = self.dens*self.ext/self.dens.max()
            elif self.MYSTIC_profile_file_flag == 3:
                print "self.dens.max()", self.dens.max()
                sys.stdout.flush()

                if self.LWC<0.0:
                    self.dens = self.dens
                else:
                    if self.dens.max() > 0.0:
                        print "self.dens.max()", self.dens.max()
                        self.dens = self.dens*self.LWC/self.dens.max()
                    else:
                        self.dens = self.dens*self.LWC

        else:
            print 'Plume set_density: Unknown plume density shape: ', shape
            exit()


        if 'pklfile' in kwargs.keys():
            pkl_file = open(kwargs['pklfile'], 'wb')
            pickle.dump(self.dens,pkl_file )
            pkl_file.close()
        if 'pklfile_z_integrated' in kwargs.keys():
            tmpInt= np.trapz(self.dens[:,:,:], self.z[1:], axis=2)
            pkl_file = open(kwargs['pklfile_z_integrated'], 'wb')
            pickle.dump(tmpInt,pkl_file )
            pkl_file.close()

        return

    def save(self,folder):
        pkl_file = open(folder+'Plume'+self.name+'.pkl', 'wb')
        pickle.dump(self,pkl_file )
        pkl_file.close()
        return

    def SetRTInput(self, Domain, UVS):
        print "Plume SetRTInput", UVS.IOdir+self.name+'.profile'
        sys.stdout.flush()

        self.ProfileInputFile = UVS.IOdir+self.name+'.profile'
        UVS.inp['profile_file '+self.randname+' 3d '] = self.ProfileInputFile
        if self.cdf == '':
            UVS.inp['profile_properties '+self.randname] = 'mie interpolate'
        else:
            UVS.inp['profile_properties '+self.randname] = self.cdf+' interpolate'

        # These large arrays are only local as this is the only place they
        # are needed.
        if self.MYSTIC_profile_file_flag==1:
            tmpext = np.zeros([Domain.nx,Domain.ny, len(Domain.z)])
            tmpgg  = np.zeros([Domain.nx,Domain.ny, len(Domain.z)])
            tmpssa = np.zeros([Domain.nx,Domain.ny, len(Domain.z)])

            tmpext[self.ix_start:self.ix_end,self.iy_start:self.iy_end,self.iz_start:self.iz_end] = self.dens
            tmpgg  = np.where(tmpext>0, self.gg, 0.0)
            tmpssa = np.where(tmpext>0, self.ssa, 0.0)

            args = {
                'ext': tmpext,
                'gg' : tmpgg,
                'ssa': tmpssa
                }
        elif self.MYSTIC_profile_file_flag==3:
            tmpLWC = np.zeros([Domain.nx,Domain.ny, len(Domain.z)])
            tmpreff  = np.zeros([Domain.nx,Domain.ny, len(Domain.z)])

            tmpLWC[self.ix_start:self.ix_end,self.iy_start:self.iy_end,self.iz_start:self.iz_end] = self.dens
            tmpreff  = np.where(tmpLWC>0, self.reff, 1.0)

            args = {
                'LWC': tmpLWC,
                'reff' : tmpreff,
                }

        Write3DMYSTICFile(self.ProfileInputFile, type='Generic',
                          nx = Domain.nx, ny = Domain.ny, nz = Domain.nz,
                          dx = Domain.dx, dy = Domain.dy, z = Domain.z, LWCLimit=self.LWCLimit,
                          flag=self.MYSTIC_profile_file_flag, **args)
        return

    def ReadLESNetCDF(self, fn, timestep=0, verbose=False, scale_axis_factor=1.,
                      ROI=False, xlimits=(0,0), zlimits=(0,0), Qbl=-999):
        """
        Qbl = SO2 emission rate in kg/s, see note by Soon-Young
        """

        CF=ConversionFactors()

        if verbose:
            print "ReadLESNetCDF: fn:", fn
            sys.stdout.flush()

        ncfile = Dataset(fn,'r')
        if verbose:
            print ncfile.variables
            sys.stdout.flush()

        self.shape='User'
        self.time   = ncfile.variables['time'][:]
        self.z  = ncfile.variables['zu_3d'][:]*CF.m2km*scale_axis_factor
        self.x  = ncfile.variables['x'][:]*CF.m2km*scale_axis_factor
        self.y  = ncfile.variables['y'][:]*CF.m2km*scale_axis_factor
        self.dens  = ncfile.variables['s'][:]
        self.dens = self.dens[timestep,:,:,:]
        #        print "self.dens", self.dens.shape
        # LES densities are z,y,x, change to x,y,z
        self.dens  = np.swapaxes(self.dens,0,2)
        #        print "self.dens", self.dens.shape
        #        print "self.z", scale_axis_factor, CF.m2km, self.z
        if ROI:
            print "self.dens", self.dens.shape, self.x.shape, self.z.shape
            sys.stdout.flush()

            if xlimits[0]!=xlimits[1]:
                print "GABBA", xlimits
                indx = np.where((self.x >=xlimits[0] ) & (self.x <= xlimits[1]))
                self.dens = self.dens[indx[0][0]:indx[0][len(indx[0])-1],:,:]
                self.x  = self.x[indx[0][0]:indx[0][len(indx[0])-1]]
                print "self.dens", self.dens.shape, self.x.shape
            if zlimits[0]!=zlimits[1]:
                print "GABBA", zlimits
                indz = np.where((self.z >=zlimits[0] ) & (self.z <= zlimits[1]))
                self.dens = self.dens[:,:, indz[0][0]:indz[0][len(indz[0])-1]]
                self.z  = self.z[indz[0][0]:indz[0][len(indz[0])-1]]
                print "self.dens", self.dens.shape, self.x.shape, self.z.shape


        self.x_start=self.x[0]
        self.x_length=self.x[len(self.x)-1]-self.x_start
        self.y_start=self.y[0]
        self.y_length=self.y[len(self.y)-1]-self.y_start
        self.z_start=self.z[0]
        self.z_length=self.z[len(self.z)-1]-self.z_start

        if Qbl>0.0:
            from scipy.constants import N_A
            # 0.00001024 = (H_wt**2)/(H_bl**2) , see note by Soon-Young
            # SO2: 64.066 g mol-1
            molweight=64.066
            ConcFact=Qbl*0.00001024*CF.kgtog*N_A/(CF.m3tocm3*molweight)
            self.dens=self.dens*ConcFact
            print "Qbl=", Qbl, N_A, self.dens.max()


        ncfile.close()
        return

    def test(self):
        print "#################################################################################"
        print "TEST TEST TEST PLUME DENSITY CHANGED"
        print "#################################################################################"
        #        print self.z
        self.dens = self.dens*0.0
        indz = np.where((self.z > 0.05) & (self.z < 0.300))
        #        indy = np.where((self.y > 1.6) & (self.y < 1.8))
        indy = np.where((self.y > .15) & (self.y < .35))
        indx = np.where((self.x > 0.502) & (self.x < 1.4))
        #        print "indz", indz[0]
        #        print "indy", indy[0]
        #        print "indx", indx[0]
        self.dens[indx[0][0]:indx[0][len(indx[0])-1], indy[0][0]:indy[0][len(indy[0])-1], indz[0][0]:indz[0][len(indz[0])-1]]=0.0000025
        #        indz = np.where((self.z > 0.475) & (self.z < 0.625))
        #        indz = np.where((self.z > 0.1) & (self.z < 0.3))
        #        self.dens[indx[0][0]:indx[0][len(indx[0])-1], indy[0][0]:indy[0][len(indy[0])-1], indz[0][0]:indz[0][len(indz[0])-1]]=0.000000025
#        indx = np.where((self.x > 1.2) & (self.x < 1.4))
#        self.dens[indx[0][0]:indx[0][len(indx[0])-1], indy[0][0]:indy[0][len(indy[0])-1], indz[0][0]:indz[0][len(indz[0])-1]]=0.0
        indx = np.where((self.x > .65) & (self.x < 1.0))
        self.dens[indx[0][0]:indx[0][len(indx[0])-1], indy[0][0]:indy[0][len(indy[0])-1], indz[0][0]:indz[0][len(indz[0])-1]]=0.0
        # indy = np.where((self.y > 1.7) & (self.y < 1.8))
        # self.dens[indx[0][0]:indx[0][len(indx[0])-1], indy[0][0]:indy[0][len(indy[0])-1], indz[0][0]:indz[0][len(indz[0])-1]]=0.00025


        indn = np.where(self.dens>0.0)
        #        print "len(indn[0])", len(indn[0])
        indn = np.where(self.dens>0.5)
        #        print "len(indn[0])", len(indn[0])


    def shiftxyz(self, x_shift=0.0, y_shift=0.0, z_shift=0.0):
        """
        Shift location of plume in 3D. All units in km
        """
        self.x = self.x+x_shift
        self.y = self.y+y_shift
        self.z = self.z+z_shift

        self.x_shift = x_shift
        self.y_shift = y_shift
        self.z_shift = z_shift

        self.x_start=self.x[0]
        self.x_length=self.x[len(self.x)-1]-self.x_start
        self.y_start=self.y[0]
        self.y_length=self.y[len(self.y)-1]-self.y_start
        self.z_start=self.z[0]
        self.z_length=self.z[len(self.z)-1]-self.z_start

        return

class Experiment():
    def __init__(self, home):
        self.home = home
        self.uvspecpath = ''
        self.verbose=False
        self.CF     = ConversionFactors()
        self.Domain = Domain()
        self.UVS    = UVspec()
        self.PlumeList = []
        self.CameraList = []
        self.SpectrometerList = []
        self.n_processes=1
        self.folder = './tmp/'
        self.RandString= ''.join(random.sample(string.ascii_lowercase, 5))
        self.UVS.RandString=self.RandString
        return

    def finalize(self):
        try:
            os.makedirs(self.folder)
        except OSError:
            if os.path.exists(self.folder):
                # We are nearly safe
                pass
            else:
                # There was an error on creation, so make sure we know about it
                print 'Could not create {:s}'.format(self.folder)
                raise

        self.UVS.IOdir = self.folder

        return

    def save(self, FileName=''):
        if FileName=='':
            pkl_file = open(self.folder+'Experiment.pkl', 'wb')
        else:
            pkl_file = open(FileName, 'wb')
        pickle.dump(self,pkl_file )
        pkl_file.close()

        self.Domain.save(self.folder)

        for Plu in self.PlumeList:
            Plu.save(self.folder)

        for Cam in self.CameraList:
            Cam.save(self.folder)

        return

    def Run(self, Instrument, RunUvspec=True, OldOutFiles=[], verbose=True, Wait=False):
        if verbose:
            print 'Experiment.Run InputFile {:s}'.format(self.InputFile)
            sys.stdout.flush()

        Convolve=False
        if Instrument.Type=='Spectrometer':
            try :
                Instrument.slitfunction
                Convolve=True
            except:
                Convolve=False
        if Instrument.Type=='Camera':
            try :
                Instrument.filterfunction
                Convolve=True
            except:
                Convolve=False

        if verbose:
            print "Experiment:Run Convolve", Convolve

        if Convolve and Instrument.Type=='Spectrometer':
            self.OutputFile = self.InputFile.replace('.inp','.out')
        else:
            self.OutputFile = self.InputFile.replace('.inp','.out')
        if verbose:
            print 'Experiment.Run OutputFile {:s}'.format(self.OutputFile)
            sys.stdout.flush()

        OutputFiles=[]
        tmpoutputfiles=[]
        (tmpoutputfiles, OutputFiles)=self.UVS.Run(self.InputFile, self.OutputFile, n_processes=self.n_processes,
                                 uvspecpath=self.UVS.uvspecpath, RandString=Instrument.RandString,
                                 Convolve=Convolve, Instrument=Instrument, OldOutFiles=OldOutFiles,
                                 RunUvspec=RunUvspec, Wait=Wait)

        if RunUvspec and Convolve:
            for fnto, fnfrom in zip(OutputFiles, tmpoutputfiles):
                print "fnto, fnfrom", fnto, fnfrom
                sys.stdout.flush()
                shutil.copy(fnfrom,fnto)

        return OutputFiles

    def WriteRTInput(self, Instrument, InputFile=None,verbose=False):
        # Make new UVS for each call to this function to avoid
        # keeping old input parameters from previous calls.
        UVS = copy.deepcopy(self.UVS)

        if verbose:
            print "Exp WriteRTInput"
            sys.stdout.flush()

        Instrument.SetRTInput(UVS)
        for Plume in Instrument.PlumeList:
            print "Exp WriteRTInput", Instrument.name, Plume.name, Plume.dens.shape
            sys.stdout.flush()
            Plume.SetRTInput(self.Domain, UVS)
        if not 'mc_minphotons' in self.UVS.inp:
            UVS.inp['mc_minphotons'] = UVS.inp['mc_photons']

        #        if self.n_processes==1:
        UVS.inp['mc_std'] = ''

        try:
            self.Domain.elevation
            nx=2
            ny=2
            fn = self.folder+'MYSTIC2DElevation.dat'
            Write2DMYSTICElevationFile(fn, nx, ny, self.Domain.x_size, self.Domain.y_size, self.Domain.elevation)
            UVS.inp['mc_elevation_file'] = fn
        except:
            pass

        if InputFile==None:
            InputFile = self.UVS.IOdir+'uvspec'+Instrument.name+'.inp'

        self.InputFile=InputFile
        print "Writing uvspec input to file", self.InputFile
        sys.stdout.flush()
        UVS.WriteInputFile(self.InputFile)
        return

class UVspec:
    def __init__(self):
        # Set some uvspec input that most likely will stay the same
        self.IOdir = './'
        self.inp = {}
        self.inp["mc_backward"] = ''
        self.inp["mc_vroom"] = 'on'
        self.inp["albedo"] = '0.0'

        return

    def add_mc_basename_to_input_file(self,mc_basename,fn):
        f = open(fn,'a')
        f.write('{0:s}\n'.format('mc_basename '+mc_basename))
        f.close()

    def info(self):
        print
        print 'UVspec parameters:'
        print 'IOdir: {:s}'.format(self.IOdir)
        sys.stdout.flush()

        return

    def Run(self,inp, out, verbose=True, n_processes=2, uvspecpath='', RandString='gabba',
            Convolve=False, Instrument=None, OldOutFiles='',RunUvspec=True, Wait=False):
        debug=False # True #
        if verbose:
            if RunUvspec:
                print "Running uvspec with input file: ", inp
            else:
                print "NOT running uvspec, using old output file: ", inp
            print "Output to file                : ", out
            print "Number of processors          : ", n_processes
            print "Convolve                      : ", Convolve
            sys.stdout.flush()

        tmp_out_base = 'tmp_mystic_'+RandString+'.out_'
        tmp_inp_base = 'tmp_mystic_'+RandString+'.inp_'
        # Remove all old files
        # OR NOT: Keep and remove manually in order to be able
        # inspect problems in *.err files, AK 20160526
        #FIXME
        # if RunUvspec:
        #     #        for filename in glob('gabba'+tmp_out_base+"*"):
        #     for filename in glob(tmp_out_base+"*"):
        #         if not debug:
        #             os.remove(filename)
        #     for filename in glob(tmp_inp_base+"*"):
        #         #        for filename in glob('gabba'+tmp_inp_base+"*"):
        #         if not debug:
        #             os.remove(filename)

        if RunUvspec:
            jobs = []
            tmpinputfiles=[]
            tmpoutputfiles=[]
            for i in range(n_processes):
                # Copy input file to temporary input file to be able to add different
                # mc_basenames to the file without destroying the input file
                tmp_inp = tmp_inp_base+str(i)
                tmpinputfiles.append(tmp_inp)
                cmd = 'cp '+inp+' '+tmp_inp
                Popen([r"cp",inp, tmp_inp]).wait()
                mc_basename = tmp_out_base+'NP_'+str(i)
                self.add_mc_basename_to_input_file(mc_basename,tmp_inp)
                tmp_out = tmp_out_base+str(i)
                print "tmp_out", tmp_out
                ips = '{:d}'.format(i)
                tmpoutputfile = tmp_out.replace('out_'+ips,'out_NP_'+ips)+'.rad.spc'
                print "tmpoutputfile", tmpoutputfile
                tmpoutputfiles.append(tmpoutputfile)
                if verbose:
                    print 'Starting process:',i,' inp:',tmp_inp,' out:',tmp_out
                    sys.stdout.flush()

                if not debug:
                    if RunUvspec:
                        p = multiprocessing.Process(target=self.worker, args=(tmp_inp,tmp_out,uvspecpath))
                        jobs.append(p)
                        p.start()
            for j in jobs:
                j.join()
        else:
            tmpoutputfiles=OldOutFiles

        if verbose:
            print 'All processes done. Read output, convolve, average and calculate std.'
            sys.stdout.flush()


        if Wait:
            print "Waiting ....."
            sys.stdout.flush()
            time.sleep(60*3) # Sleep for 3 minutes to assure that files are put in right place

        if Convolve:
            finalrawoutputfiles=[]
            tmpfilestoaverage=[]
            if Instrument.Type=='Spectrometer':
                # Convolve with slit function if given.
                if verbose:
                    print 'Convolving with slit function:', Instrument.slitfunction
                    sys.stdout.flush()

                ip=0
                for tmpoutputfile in tmpoutputfiles:
                    ips = '{:d}'.format(ip)
                    rawoutputfile = inp.replace('.inp','.out_NP_'+ips+'.rad.spc')
                    print tmpoutputfile, rawoutputfile
                    sys.stdout.flush()
                    finalrawoutputfiles.append(rawoutputfile)
                    tmpoutconv='tmpoutconv_'+Instrument.RandString+'_'+ips
                    cmd = '/usr/bin/time -v '+self.uvspecpath+'conv '+tmpoutputfile+' '+Instrument.slitfunction+' > '+tmpoutconv+' 2> '+tmpoutconv+'.err'
                    if verbose:
                        print cmd
                        sys.stdout.flush()
                    p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
                    tmpoutspline='tmpoutspline_'+Instrument.RandString+'_'+ips
                    cmd = '/usr/bin/time -v '+self.uvspecpath+'spline '+'-q -l -b '+str(Instrument.wavelength1)+' -s '+str(Instrument.wavelengthstep)+' '+tmpoutconv+' > ' + tmpoutspline+' 2> '+tmpoutspline+'.err'
                    if verbose:
                        print cmd
                        sys.stdout.flush()
                    p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
                    tmpfilestoaverage.append(tmpoutspline)
                    # Copy MYSTIC output files to final destination
                    shutil.copy(tmpoutputfile,rawoutputfile)
                    ip=ip+1
            elif Instrument.Type=='Camera':
                nx = Instrument.h_pixels
                ny = Instrument.v_pixels
                tmpSplineXFile='tmpSplineXFile'+Instrument.RandString
                # Any output file should do to get wavelength information, there should be at least one.
                tmpdata = np.loadtxt(tmpoutputfiles[0])
                nwl = int(tmpdata.shape[0]/(nx*ny))
                tmpdata = np.reshape(tmpdata,(nwl,nx, ny, tmpdata.shape[1]))
                # Interpolate filter function to MYSTIC output wavelengths
                fx = open(tmpSplineXFile,'w')
                wvls = tmpdata[:,0,0,0]
                for wvl in wvls:
                    fx.write('{:f}\n'.format(wvl))
                fx.close()
                tmpSplineOutputFile='tmpSplineOutputFile'+Instrument.RandString
                cmd = '/usr/bin/time -v '+self.uvspecpath+'spline '+'-q -l -x '+tmpSplineXFile+' '+Instrument.filterfunction+' > ' + tmpSplineOutputFile+' 2> '+tmpSplineOutputFile+'.err'
                if verbose:
                    print cmd
                    sys.stdout.flush()
                p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
                tmpfilterfunctionwvl, tmpfilterfunction = np.loadtxt(tmpSplineOutputFile,unpack=True)

                ###
                # Include loop over all output files.
                ###
                ip=0
                for tmpoutputfile in tmpoutputfiles:
                    ips = '{:d}'.format(ip)
                    rawoutputfile = inp.replace('.inp','.out_NP_'+ips+'.rad.spc')
                    if verbose:
                        print "tmpoutputfile, rawoutputfile", tmpoutputfile, rawoutputfile
                        sys.stdout.flush()
                    finalrawoutputfiles.append(rawoutputfile)

                    tmpdata = np.loadtxt(tmpoutputfile)
                    tmpdata = np.reshape(tmpdata,(nwl,nx, ny, tmpdata.shape[1]))

                    tmpoutputfilefilter = tmpoutputfile.replace('.out','.out_NP_'+ips+'.rad.spc')
                    if verbose:
                        print "tmpoutputfilefilter", tmpoutputfilefilter
                    tmpfilestoaverage.append(tmpoutputfilefilter)
                    f= open(tmpoutputfilefilter,'w')
                    # For each pixel
                    ix=0
                    iz=0
                    while ix<nx:
                        iy=0
                        while iy<ny:
                            ## Multiply MYSTIC radiances with filter function
                            tmprad = tmpdata[:,ix,iy,4]*tmpfilterfunction
#                            tmpstd = tmpdata[:,ix,iy,5]*tmpfilterfunction
                            ## Integrate over wavelength
                            totrad = np.trapz(tmprad, x=wvls)
#                            totstd = np.trapz(tmpstd, x=wvls)
#                            f.write('{0:8.2f} {1:3d} {2:3d} {3:3d} {4:9.4f} {5:11.6f}\n'.format(wvls[0],ix,iy,iz,totrad,totstd))
                            f.write('{0:8.2f} {1:3d} {2:3d} {3:3d} {4:9.4f}\n'.format(wvls[0],ix,iy,iz,totrad))
                            iy=iy+1
                        ix=ix+1
                    f.flush()            # Do this to make sure all os written
                    os.fsync(f.fileno()) # to file before continuing.
                    f.close()
                    ip=ip+1

        else:
            tmpfilestoaverage=tmpoutputfiles
            finalrawoutputfiles=tmpoutputfiles

        InputFiles = tmpfilestoaverage #tmp_out_base+'NP_'+'*'+'.rad.spc'
        if n_processes==1:
            if verbose:
                print "InputFiles, OutputFileRaw", InputFiles, out, tmpoutputfiles
                sys.stdout.flush()

            CombineSingleProcessOuput(tmpoutputfiles, out, verbose=True)
        else:
            if verbose:
                print "finalrawoutputfiles", finalrawoutputfiles
                print "tmpoutputfiles", tmpoutputfiles
                sys.stdout.flush()

            Average_spc_Files(InputFiles, out, verbose=True)
        return (tmpoutputfiles, finalrawoutputfiles)

    def SingleRun(self,inp, out, verbose=False, uvspecpath=''):
        if verbose:
            print "Running uvspec with input file: ", inp
            print "Output to file                : ", out
            sys.stdout.flush()

        if 'xnilu_wrk' in uvspecpath:
            uvspec='uvspec3D'
        else:
            uvspec='uvspec'
        cmd = '/usr/bin/time -v '+uvspecpath+uvspec+' < '+inp+' > '+out+' 2> '+out+'.err'
        if verbose:
            print cmd
            sys.stdout.flush()

        #(uvspec < uvspec.inp > uvspec.out) >& uvspec.err

        #FIXME
        p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
        return

    def worker(self, input,output, uvspecpath=''):
        """thread worker function"""
        verbose = True
        self.SingleRun(input,output,verbose=verbose, uvspecpath=uvspecpath)
        return

    def WriteInputFile(self, InputFile=None, verbose=False):
        if verbose:
            print "Writing uvspec input file", InputFile
            sys.stdout.flush()

        try:
            f = open(InputFile,'w')
        except:
            print "Experiment.WriteRTFile: Not able to open uvspec input file."
            exit()
        for key in self.inp:
            if verbose:
                sys.stdout.write( key + ' ' + str(self.inp[key]) + '\n')
            f.write( key + ' ' + str(self.inp[key]) + '\n')
        f.flush()     # Do this to make sure all input files are written
        os.fsync(f.fileno())    #  to file.
        f.close()
        return


#######################################################################

if __name__ == "__main__":

    Exp = Experiment()
    print Exp.Domain.name

    Exp.CameraList.append(Camera())

    for Cam in Exp.CameraList:
        Cam.info()
