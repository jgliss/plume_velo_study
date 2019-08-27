"""
Copyright 2014-2015 Arve Kylling

This file is part of IRoPLUM.

IRoPLUM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

IRoPLUM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with IRoPLUM.  If not, see <http://www.gnu.org/licenses/>.

"""
# Where are we
import os
myhost = os.uname()[1]
home = os.environ['HOME']
import numpy as np
import math
import pickle
if 'abel' not in home:
    import matplotlib.pyplot as plt
if myhost=='ak'  in home:
    from skimage.transform import warp
    from shapely.geometry import LineString

def find_nearest(array,value):
    idx = np.searchsorted(np.sort(array), value, side="left")
#FIXME
#    idx = np.searchsorted(array, value, side="right")
#    print "idx", idx
    if value > array[idx]:
        below=0
    else:
        below=1
#    print 'value', value, array[idx], array, idx, below
#    exit()
    if math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
        return idx-1
    else:
        return idx-below

def FindIndexClosestToValue(t, x):
    indsort = np.argsort(t)
    ind = t.searchsorted(x,sorter=indsort)
    if ind == len(t):
        return ind - 1 # x > max(t)
    elif ind == 0:
        return 0
    before = ind-1
    if x-t[before] < t[ind]-x:
        ind -= 1
    return ind

class Area:
    def __init__(self, nx=None, dx=None, x0=None,
                 nz=None, dz=None, z0=None, Empty=False):

        if Empty:
            return
        self.nx=nx
        self.dx=dx
        self.x0=x0
        self.nz=nz
        self.dz=dz
        self.z0=z0
        self.xmin = x0 - dx*nx/2
        self.xmax = x0 + dx*nx/2
        self.zmin = z0
        self.zmax = z0 + dz*nz
        # Grid given by
        self.x = np.arange(self.xmin,self.xmax+self.dx,self.dx)
#        print self.x
#        self.x = self.x[::-1]
        self.z = np.arange(self.zmin,self.zmax+self.dz,self.dz)
        # FIXME
        self.z = self.z[::-1]
#        print self.z
        self.x_mid = np.arange(self.xmin+dx/2.,self.xmax+self.dx/2,self.dx)
        self.z_mid = np.arange(self.zmin+dz/2.,self.zmax+self.dz/2,self.dz)
        self.x_center = 0.5*(self.x[0]+self.x[len(self.x)-1])
        self.z_center = 0.5*(self.z[0]+self.z[len(self.z)-1])
        self.x_center_ind = FindIndexClosestToValue(self.x,self.x_center)
        self.z_center_ind = FindIndexClosestToValue(self.z,self.z_center)
        self.MaxLength = np.sqrt((dx*nx)**2 + (dz*nz)**2)
#        print "self.MaxLength", self.MaxLength

        # Pixels are defined by grid. So there is one less pixels than
        # grid values in each direction
        self.Image_nx=self.nx
        self.Image_nz=self.nz
        self.Image = np.ones((self.Image_nz, self.Image_nx))
        return

    def save(self, fn='tmpReconstructed'):
        output = open(fn, 'wb')
        pickle.dump(self, output)
        output.close()

    def read(self, fn):
        pf = open(fn, 'rb')
        self =  pickle.load(pf)
        pf.close()
        return self

    def WriteReconstructedTiff(self, fn):
        fp=open(fn,'w')
        fp.close()
        return

    def sumImage(self):
        return np.sum(self.Image)

    def plot(self, ax, cmap='gray'):
        # generate 2 2d grids for the x & y bounds
#        print "self.Image"
#        print self.Image
        tmpdx = (self.xmax-self.xmin)/(self.nx)
        tmpdz = (self.zmax-self.zmin)/(self.nx)
        z, x = np.mgrid[slice(self.zmin, self.zmax+tmpdz/2 , tmpdz),
                        slice(self.xmin, self.xmax+tmpdx/2 , tmpdx)]
        z_min, z_max = -np.abs(z).max(), np.abs(z).max()
        z_min, z_max = -np.abs(self.Image).max(), np.abs(self.Image).max()
        im=ax.pcolor(x, z, self.Image, cmap=cmap,alpha=0.5) #, vmin=z_min, vmax=z_max)
        return im

class PlotProp:
    def __init__(self):
        self.color='r'
        self.text_x = 0.0
        self.text_y = 0.0
        return

class TomoCamera:
    def __init__(self, x=0, y=0, z=0,
                 theta1=None,theta2=None, Nrays=None, Name='',
                 CameraFile=None, Empty=False,
                 verbose=False):
        """
        phi1, phi2, theta1, and theta2 are defined as for the
        libRadtran panoramaview option.
        phi is azimuth and theta elevation. theta=90 is horizon.
        theta larger than 90 is above horizon.
        """
#        verbose=True
        if Empty:
            self.Sinogram=np.empty(0)
            self.Rays_x=np.empty((2,0))
            self.Rays_y=np.empty((2,0))
            self.Angle=np.empty(0)
            self.Rq=np.empty(0)
            self.TotalLength=np.empty(0)
            self.a=np.empty(0)
            self.b=np.empty(0)
            self.PlotProp = PlotProp()
            self.ind=[]
            self.Name = 'Empty'
            return

        # nx and ny should be odd numbers to put plume in center
        # dx, dy, and dz are in km
        self.PlotProp = PlotProp()
        self.verbose=False
        if CameraFile==None:
            self.x=float(x)
            self.y=float(y)
            self.z=float(z)
            self.theta1=theta1
            self.theta2=theta2
            self.Nrays = Nrays
            self.Name = Name
            self.Rq=np.zeros(self.Nrays)
            self.Sinogram=np.zeros(self.Nrays)
            self.TotalLength=np.zeros(self.Nrays)
            self.dtheta = np.abs((self.theta2-self.theta1))/float((self.Nrays-1))
            self.Angle =np.zeros(self.Nrays)
            i=0
            while i<self.Nrays:
                angle = self.theta1+self.dtheta*float(i)
                self.Angle[i] = angle
                i=i+1
            self.ind=[]
        else:
            self.CameraFile = CameraFile
            self.ReadCameraFile(self.CameraFile)

        self.Name=Name
        # FIXME
        self.l = 20000   # Just some large number in meters for plotting
                         # and ray construction

        if verbose:
            print "theta1, theta2", self.theta1, self.theta2

        return

    def ReadCameraFile(self,CameraFile, verbose=False):

        fp = open(CameraFile)
        tmp = fp.readline()
        tmp = tmp.replace('\n','')
        self.Name = tmp
        tmp = fp.readline()
        tmp = tmp.replace('\n','')
        tmpx, tmpy, tmpz=tmp.split()
        self.x = float(tmpx);
        self.y = float(tmpy);
        self.z = float(tmpz);
        tmp = fp.readline()
        tmp = tmp.replace('\n','')
        tmpx, tmpy, tmpz=tmp.split()
        self.Nrays = int(tmpx)
        self.theta1 = float(tmpy)
        self.theta2 = float(tmpz)
        fp.close()
        tmpdat=np.loadtxt(CameraFile,skiprows=3)
        self.Angle = tmpdat[:,0]
        self.Nrays = len(self.Angle)
        self.Sinogram = tmpdat[:,1]
        self.ind=[]
        return

    def Rays(self, verbose=False):
        self.Rays_x=np.zeros((2,self.Nrays))
        self.Rays_y=np.zeros((2,self.Nrays))
        self.a=np.zeros((self.Nrays))
        self.b=np.zeros((self.Nrays))
        x0 = self.x
        y0 = self.z
        l=self.l

#        verbose=True

        i=0
        while i<self.Nrays:
            angle = self.Angle[i]
#            angle = self.theta1+self.dtheta*float(i)
            radang = np.deg2rad(angle)
            x1 = x0+l*np.cos(radang)
            y1 = y0+l*np.sin(radang)
            self.Rays_x[0,i] = x0
            self.Rays_x[1,i] = x1
            self.Rays_y[0,i] = y0
            self.Rays_y[1,i] = y1
            self.a[i] = (y1-y0)/(x1-x0)
            self.b[i] = y0-self.a[i]*x0

            if verbose:
                print "Rays", i, angle, x0, y0, x1, y1
            i=i+1

        return

    def PixelsCoveredByRays(self, RR,plot=False):
        #FIXME not doing anything as is.

        tmpimage = RR.Image*0
        iRay=0
        while iRay < self.Nrays:
            tmpimage[self.ind[iRay]] = 1.0
            iRay=iRay+1

        if plot:
            fig, (ax1) = plt.subplots(1, 1, figsize=(4, 4))
            tmpdx = (RR.xmax-RR.xmin)/(RR.nx)
            tmpdz = (RR.zmax-RR.zmin)/(RR.nx)
            z, x = np.mgrid[slice(RR.zmin, RR.zmax+tmpdz/2 , tmpdz),
                            slice(RR.xmin, RR.xmax+tmpdx/2 , tmpdx)]
            z_min, z_max = -np.abs(z).max(), np.abs(z).max()
            z_min, z_max = -np.abs(RR.Image).max(), np.abs(RR.Image).max()
            im=ax1.pcolor(x, z, tmpimage, cmap='RdBu') #, vmin=z_min, vmax=z_max)
            plt.show()

        self.PixelsCoveredByRaysImage = tmpimage
#        self.PixelsCoveredByRaysImage = RR.Image*1.0

        return self.PixelsCoveredByRaysImage

    def CalculateLineIntegral(self, RR, iRay, Plot=False, quiet=False, verbose=False,
                              debug=False, ax1=None, iPlotRay=0):

#        Plot=False  # Does not work
#        if Plot:
#            xsize=ysize=6
#            fig, (ax1) = plt.subplots(1, 1, figsize=(ysize, xsize))
#            fig.subplots_adjust(top=0.90, bottom=0.18, left=0.18, right=0.95)

        RRn = len(RR.x)-1
        ### Find entry side
        Ray = LineString([(self.Rays_x[0,iRay],self.Rays_y[0,iRay]),
                          (self.Rays_x[1,iRay],self.Rays_y[1,iRay])])
        # Left side
        BoxLineA = LineString([(RR.x[0],RR.z[RRn]),(RR.x[0],RR.z[0])])
        # Right side
        BoxLineB = LineString([(RR.x[RRn],RR.z[RRn]),(RR.x[RRn],RR.z[0])])
        # Bottom
        BoxLineC = LineString([(RR.x[0],RR.z[RRn]),(RR.x[RRn],RR.z[RRn])])
        # Top
        BoxLineD = LineString([(RR.x[0],RR.z[0]),(RR.x[RRn],RR.z[0])])

        #            if verbose:
        #                print Ray
        #                print BoxLineA
        #                print BoxLineB
        #                print BoxLineC
        #                print BoxLineD

        intersectA = Ray.intersection(BoxLineA)
        intersectB = Ray.intersection(BoxLineB)
        intersectC = Ray.intersection(BoxLineC)
        intersectD = Ray.intersection(BoxLineD)

        # Find the two intersection points
        intsec1=None
        side1=''
        if intersectA.geom_type=='Point': intsec1 = intersectA; side1='L'
        if intersectB.geom_type=='Point' and intsec1==None:
            intsec1 = intersectB
            side1 = 'R'
        elif intersectB.geom_type=='Point':
            intsec2 = intersectB
            side2 = 'R'
        if intersectC.geom_type=='Point' and intsec1==None:
            intsec1 = intersectC
            side1 = 'B'
        elif intersectC.geom_type=='Point':
            intsec2 = intersectC
            side2 = 'B'
        if intersectD.geom_type=='Point' and intsec1==None:
            intsec1 = intersectD
            side1 = 'T'
        elif intersectD.geom_type=='Point':
            intsec2 = intersectD
            side2 = 'T'

#        print intersectA.geom_type, intsec1, side1, side2
#        print intersectB.geom_type, intsec1, side1, side2
#        print intersectC.geom_type, intsec1, side1, side2
#        print intersectD.geom_type, intsec1, side1, side2
        if intsec1==None:
            if not quiet:
                print "Ray does not hit reconstruction area"
            tmpimage = RR.Image * 0.0
            ind = np.where(tmpimage>0)
            return 0, 0.0, ind, 0

        # Decide which side the ray entered first.
        class Entry:
            def __init__(self):
                return
        entry=Entry()

        if side1=='L' and side2=='T':
            if self.Rays_x[0,iRay] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='L' and side2=='R':
            if self.Rays_x[0,iRay] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='L' and side2=='B':
            if self.Rays_x[0,iRay] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='R' and side2=='T':
            if self.Rays_x[0,iRay] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='R' and side2=='L':
            if self.Rays_x[0,iRay] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='R' and side2=='B':
            if self.Rays_x[0,iRay] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='B' and side2=='T':
            if self.Rays_y[0,iRay] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='B' and side2=='L':
            if self.Rays_y[0,iRay] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2
        elif side1=='B' and side2=='R':
            if self.Rays_y[0,iRay] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
            else: entry.side=side2; entry.intsec=intsec2

        ### Set coordinates of entry point
        if verbose:
            print "entry.side", entry.side
        if entry.side=='L':
            ix = 0
            x0 = RR.x[ix]
            z0 = self.a[iRay]*x0 + self.b[iRay]
#FIXME
#            za = np.where(RR.z > z0)
            za = np.where(RR.z < z0)
            iz = za[0][0]-1
            sign=+1
            offset=0
            zoffset=0
        elif entry.side=='R':
            offset=1
            zoffset=0
            ix = RRn
            x0 = RR.x[ix]
            z0 = self.a[iRay]*x0 + self.b[iRay]
#            za = np.where(RR.z > z0)
#FIXME
            za = np.where(RR.z < z0)
            iz = za[0][0]-1
            sign=-1
        elif entry.side=='B':
            offset=0
            zoffset=1
            iz=RRn #0
            z0 = RR.z[iz]
            x0  = (z0-self.b[iRay])/self.a[iRay]
#FIXME            if self.a[iRay]>0:
            if self.a[iRay]<0:
#                sign=+1
                sign=0
#                print "x0, RR.x",x0, RR.x, self.a[iRay]
                xa = np.where(RR.x > x0)
#                print xa
                ix = xa[0][0]-1
#                ix = xa[0][len(xa[0])-1]
#FIXME            elif self.a[iRay]<0:
            elif self.a[iRay]>0:
#                sign=-1
#                print "x0, RR.x",x0, RR.x, self.a[iRay]
                sign=0
                xa = np.where(RR.x < x0)
                ix = xa[0][len(xa[0])-1]
#                ix = xa[0][0]

#        print entry.side
#        print "iz, za, ix, xa", iz,ix-offset, z0, ix, x0,  entry.side, ix+sign*1
#        print "Z0,X0", z0, x0
        MaxLength   = RR.MaxLength
        TotalLength = 0
        TotalColumn = 0
        ib=0
        tmpimage = RR.Image * 0.0
        tmpimageN = RR.Image * 0.0
        tmpN=0
        if debug:
            print "###################################################"
            print "Angle", self.Angle[iRay], entry.side, iz, ix
            print RR.Image
        while TotalLength < MaxLength:
            if debug:
                print "--------------------------------------------------"
            if entry.side=='L' or entry.side=='R':
                x1 = RR.x[ix+sign*1]
                z1 = self.a[iRay]*x1 + self.b[iRay]
            elif entry.side=='B':
                if debug:
                    print "iz", iz
                z1 = RR.z[iz-1]
                x1 = ( z1 -  self.b[iRay])/self.a[iRay]

            if debug:
                print "Z0,X0", z0, x0
                print "Z1,X1", z1, x1

            # Mark traversed pixels
            ixp=ix
            izp=iz
#            print "iz-zoffset,ix-offset", iz-zoffset,ix-offset, z1, x1, ix+sign*1, RR.x
#            print "iz,ix", iz, ix, izp, ixp
            tmpimage[iz-zoffset,ix-offset]=1
#            print "z0, RR.z", z0, RR.z
#FIXME
            if entry.side=='L':
#                print entry.side
#                za  = np.where(RR.z <= z0)
#                zai = za[0][0]
#                if debug:
                if self.a[iRay] > 0.0 :
                    za  = np.where(RR.z > z0)
                    if z0 <= RR.z[len(RR.z)-1]:
                        zai = len(RR.z)-1
                    else:
                        zai = za[0][len(za[0])-1]
                elif self.a[iRay] < 0.0 :
                    za  = np.where(RR.z <= z0)
                    zai = za[0][0]
                if debug:
                    print "z0, RR.z", z0, RR.z, za, len(za[0])-1,zai, RR.z[len(RR.z)-1]
                    print za[0]
                #            za  = np.where(RR.z > z0)
                #            zai = za[0][0]
#                zai = za[0][len(za[0])-1]
                xa  = (RR.z[zai]-self.b[iRay])/self.a[iRay]
                if debug:
                    print "zai", zai, RR.z[zai]
                    print entry.side, self.a[iRay], xa, x1, zai
                # Three possibilities

                # One, ray leaves bottom of pixel before x1 is reached
                if (self.a[iRay] > 0.0 and xa < x1):
                    if debug:
                        print "ONE"
                        iz=iz-1
                    else:
                        iz=iz-1
                    x1 = xa
                    z1 = RR.z[zai]

#                    x1 = xa
#                    z1 = RR.z[zai]
#                    iz=iz-1
                elif self.a[iRay] < 0.0 and xa < x1:
                    if debug:
                        print "TWO"
                    # Leaving pixel before x1 reached
                    x1 = xa
                    z1 = RR.z[zai]
                    iz=iz+1
                else:
                    if debug:
                        print "THREE"
                    ix=ix+1
                if debug:
                    print "z1,x1", z1, x1

            elif entry.side=='R':
                if self.a[iRay] > 0.0 :
                    za  = np.where(RR.z > z0)
                    if z0 <= RR.z[len(RR.z)-1]:
                        zai = len(RR.z)-1
                    else:
                        zai = za[0][len(za[0])-1]
                elif self.a[iRay] < 0.0 :
                    za  = np.where(RR.z > z0)
                    if debug:
                        print "za", za[0], RR.z, z0, RR.z[0]
                    if z0 >= RR.z[0]:
                        zai = 0
                    else:
                        zai = za[0][len(za[0])-1]
#                    zai = za[0][0]
                if debug:
#                    za  = np.where(RR.z <= z0)
                    print "z0, RR.z", z0, RR.z, za, len(za[0])-1
                    print za[0]
                #            za  = np.where(RR.z > z0)
                #            zai = za[0][0]
#                zai = za[0][len(za[0])-1]
                xa  = (RR.z[zai]-self.b[iRay])/self.a[iRay]
                if debug:
                    print "zai", zai, RR.z[zai]
                    print entry.side, self.a[iRay], xa, x1, zai
                # Three possibilities

                # One, ray leaves bottom of pixel before x1 is reached
                if (self.a[iRay] > 0.0 and xa < x1):
                    x1 = xa
                    z1 = RR.z[zai]
                    iz=iz+1
                elif self.a[iRay] < 0.0 and xa > x1:
                    # Leaving pixel before x1 reached
                    x1 = xa
                    z1 = RR.z[zai]
                    iz=iz-1
                else:
                    ix=ix-1
                if debug:
                    print "z1,x1", z1, x1

            elif entry.side=='B':
                if self.a[iRay] > 0.0 :
                    xa  = np.where(RR.x > x0)
                    if x0 >= RR.x[len(RR.x)-1]:
                        xai = len(RR.x)-1
                    else:
                        xai = xa[0][0]
                elif self.a[iRay] < 0.0 :
                    xa  = np.where(RR.x < x0)
                    if x0 <= RR.x[0]:
                        xai = 0
                    else:
                        xai = xa[0][len(xa[0])-1]
#                    xai = xa[0][0]
                if debug:
                    print "xa", x0, RR.x, xa, xa[0]
                    print "x0, RR.x", x0, RR.x, xa, len(xa[0])-1
                    print xa[0]
                    print "len(xa[0])-1", len(xa[0])-1, xai
                    print "xai", xai, RR.x[xai], xa
                za  = RR.x[xai]*self.a[iRay]+self.b[iRay]
                if debug:
                    print "Hmm", entry.side, self.a[iRay], za, z1, xai, RR.x[xai]
                # Three possibilities
                if (self.a[iRay] > 0.0 and za < z1):
                    # One, ray leaves side of pixel before z1 is reached
                    x1 = RR.x[xai]
                    z1 = za
                    if debug:
                        print "z1,x1", z1, x1
                    ix=ix+1
                elif self.a[iRay] < 0.0 and za < z1:
                    # Two, ray leaves other side of pixel before z1 is reached
                    x1 = RR.x[xai]
                    z1 = za
                    ix=ix-1
                else:
                    # Three, ray leaves top of pixel.
                    iz=iz-1

            Length = np.sqrt((z1-z0)**2+(x1-x0)**2)
            if RR.Image[izp-zoffset,ixp-offset]>0.0:
                TotalLength = TotalLength + Length
                tmpimageN[izp-zoffset,ixp-offset]=tmpimageN[izp-zoffset,ixp-offset]+1
                tmpN=tmpN+1

#FIXME            Column = RR.Image[izp-zoffset,ixp-offset]*Length/RR.dx
            Column = RR.Image[izp-zoffset,ixp-offset]*Length

            TotalColumn = TotalColumn + Column

#            print "TotalColumn", RR.Image[izp-zoffset,ixp-offset], izp-zoffset,ixp-offset, Length, Column, TotalColumn

#            TotalColumn = TotalColumn + Length*RR.Image[izp,ixp-offset]
#            print "TC", izp,ixp+offset, TotalColumn, RR.Image[izp,ixp-offset]
#            print "Length", TotalLength, MaxLength, ixp, ix, izp, iz, ib, RRn, entry.side
#            print "Length", z0, z1, x0, x1, entry.side, Length, TotalLength, MaxLength, ix, RRn, iz
            if debug:
                print "TotalColumn", TotalColumn, iz, ix, izp-zoffset,ixp-offset, RR.Image[izp-zoffset,ixp-offset], Column
                print "Length", z0, z1, x0, x1, entry.side, Length, TotalLength, MaxLength, ix, RRn, iz, Length/RR.dx
                print "entry,side", entry.side, ix, iz, RRn
            if Plot and iRay==iPlotRay and ib<60:
                if debug:
                    print "xg", x0, x1, z0, z1, ib
                color='w'
                ax1.plot([x0,x1],[z0,z1], color=color)

            if debug:
                if entry.side=='L' and (ix >= RRn or iz < 0): break
            else:
                if entry.side=='L' and (ix >= RRn or iz < 0): break
#                if entry.side=='L' and (ix >= RRn or iz >= RRn): break
            if entry.side=='R' and (ix <=  0  or iz < 0): break
            if entry.side=='R' and (ix <=  0  or iz >= RRn): break
            if entry.side=='B' and (ix <   0  or iz < 0): break
            if entry.side=='B' and (ix >= RRn or iz <= 0): break
#FIXME            if entry.side=='B' and (ix <   0  or iz >= RRn): break
#FIXME            if entry.side=='B' and (ix >= RRn or iz >= RRn): break


            x0 = x1
            z0 = z1

            ib=ib+1
        #    if ib >= 13: break
#        print "tmpimage", iRay, self.Angle[iRay], tmpN
#        print tmpimage
#        print tmpimageN
#        if debug:
#            exit()
#FIXME        ind = np.where(tmpimage>0)
        ind = np.where(tmpimageN>0)
        self.Rq[iRay] = TotalColumn
#FIXME        TotalLength=TotalLength/RR.dx
        TotalLength=TotalLength
        self.TotalLength[iRay] = TotalLength
#        print "TotalLength", TotalLength
        self.ind.append(ind)

        return TotalColumn, TotalLength, ind, tmpN


    def LineIntegral(self, image, RR, plot=False):
        verbose= False #False #

        iPlotRay=0
        if plot:
            fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(8, 4))
            z, x = np.mgrid[slice(RR.zmin, RR.zmax+RR.dz , RR.dz),
                            slice(RR.xmin, RR.xmax+RR.dx , RR.dx)]
            z_min, z_max = -np.abs(image).max(), np.abs(image).max()
            im=ax1.pcolor(x,z,image, cmap='RdBu')
            ax1.axis([RR.xmin-RR.dx/2, RR.xmax+RR.dx/2, RR.zmin-RR.dz/2, RR.zmax+RR.dz/2])
            fig.colorbar(im)
            print "Lineintegral, image", image

        # Subtract - 0.5 from center to not shift one pixel too many in both directions
        center = image.shape[0] // 2 - 0.5
#        print "center", image.shape[0], center
        shift0 = np.array([[1, 0, -center],
                           [0, 1, -center],
                           [0, 0, 1]])
        shift1 = np.array([[1, 0, center],
                           [0, 1, center],
                           [0, 0, 1]])
        def build_rotation(theta):
            T = np.deg2rad(theta)
#            print "BD", theta, T, np.cos(T), np.sin(T)
            R = np.array([[np.cos(T), np.sin(T), 0],
                          [-np.sin(T), np.cos(T), 0],
                          [0, 0, 1]])
            return shift1.dot(R).dot(shift0)

        # Find intersection between each ray and the reconstruction region
        RRn = len(RR.x)-1
        self.Sinogram=np.zeros(self.Nrays)
        i=0
        while i<self.Nrays:
            Ray = LineString([(self.Rays_x[0,i],self.Rays_y[0,i]),
                                (self.Rays_x[1,i],self.Rays_y[1,i])])
            # Left side
            BoxLineA = LineString([(RR.x[0],RR.z[0]),(RR.x[0],RR.z[RRn])])
            # Right side
            BoxLineB = LineString([(RR.x[RRn],RR.z[0]),(RR.x[RRn],RR.z[RRn])])
            # Top
            BoxLineC = LineString([(RR.x[0],RR.z[RRn]),(RR.x[RRn],RR.z[RRn])])
            # Bottom
            BoxLineD = LineString([(RR.x[0],RR.z[0]),(RR.x[RRn],RR.z[0])])

#            if verbose:
#            print Ray
#            print BoxLineA
#            print BoxLineB
#            print BoxLineC
#            print BoxLineD

            intersectA = Ray.intersection(BoxLineA)
            intersectB = Ray.intersection(BoxLineB)
            intersectC = Ray.intersection(BoxLineC)
            intersectD = Ray.intersection(BoxLineD)

            # Find the two intersection points
            intsec1=None
            side1=''
            if intersectA.geom_type=='Point': intsec1 = intersectA; side1='L'
            if intersectB.geom_type=='Point' and intsec1==None:
                intsec1 = intersectB
                side1 = 'R'
            elif intersectB.geom_type=='Point':
                intsec2 = intersectB
                side2 = 'R'
            if intersectC.geom_type=='Point' and intsec1==None:
                intsec1 = intersectC
                side1 = 'B'
            elif intersectC.geom_type=='Point':
                intsec2 = intersectC
                side2 = 'B'
            if intersectD.geom_type=='Point' and intsec1==None:
                intsec1 = intersectD
                side1 = 'T'
            elif intersectD.geom_type=='Point':
                intsec2 = intersectD
                side2 = 'T'

            # Decide which side the ray entered first.
            class Entry:
                    def __init__(self):
                        return
            entry=Entry()
            if side1=='L' and side2=='T':
                if self.Rays_x[0,i] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='L' and side2=='R':
                if self.Rays_x[0,i] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='L' and side2=='B':
                if self.Rays_x[0,i] < RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='R' and side2=='T':
                if self.Rays_x[0,i] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='R' and side2=='L':
                if self.Rays_x[0,i] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='R' and side2=='B':
                if self.Rays_x[0,i] > RR.x[0]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='B' and side2=='T':
                if self.Rays_y[0,i] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='B' and side2=='L':
                if self.Rays_y[0,i] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2
            elif side1=='B' and side2=='R':
                if self.Rays_y[0,i] < RR.z[RRn]: entry.side=side1; entry.intsec=intsec1
                else: entry.side=side2; entry.intsec=intsec2

            # Find line parameters for this ray
            x0a = self.Rays_x[0,i]
            x1a = self.Rays_x[1,i]
            y0a = self.Rays_y[0,i]
            y1a = self.Rays_y[1,i]
            a = (y1a-y0a)/(x1a-x0a)
            b = y1a - a*x1a

            # Find coordinates for rotated line
            x0 = 0
            y0 = a*0+b

            # Find pixel id on entry side.
            x0r=RR.x[RRn/2]
            y0r=RR.z[RRn/2]
            # Below equation from
            # http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
            ratio = (x0r+a*y0r-a*b)/(a**2+1)
            d = np.sqrt((ratio -x0r)**2 + (a*ratio+b-y0r)**2)

            if entry.side=='L':
                if y0 < y0r:    d = y0r+d
                elif y0 >= y0r:    d = y0r-d
                yid = find_nearest(RR.z, d)
                if yid > RRn-1: yid=RRn-1
                pid = yid
            elif entry.side=='R':
#                print "d", d
                if y0 < y0r:    d = y0r-d
                elif y0 >= y0r:    d = y0r+d
                yid = find_nearest(RR.z, d)
                pid = yid
                print "d", d, y0, yid, pid, RR.z
            elif entry.side=='B':
                sign = a/np.sqrt(a**2)
                d = y0r-sign*d
                yid = find_nearest(RR.z, d)
                pid = yid

            # Warp image according to ray angle
            # warp by defaul clips to the range [0,1], we do not want that.
            # So normalize and denormalize after warp
            tmpimg = image/image.max()
#            print "image"
#            print image
#            print tmpimg
#            print "Rotation matrix"
#            print build_rotation(-self.Angle[i])
#            rotated = warp(tmpimg, build_rotation(-self.Angle[i]))
            rotated = warp(tmpimg, build_rotation(-self.Angle[i]))
#            print rotated
            rotated = rotated*image.max()
#            print "i",i, entry.side
#            print rotated
            if plot and i==iPlotRay:
                print "iPlotRay", iPlotRay
                print "image min", image.min(), image.max(), rotated.min(), rotated.max(), self.Angle[i]
                print image
                print rotated
                ax1.plot([x0a,x1a],[y0a,y1a])
                ax2.pcolor(x,z,rotated, cmap='RdBu')
                ax2.axis([RR.xmin-RR.dx/2, RR.xmax+RR.dx/2, RR.zmin-RR.dz/2, RR.zmax+RR.dz/2])
                x0 = RR.x.min()
                y0 = RR.z[pid]
                x1 = RR.x.max()
                y1 = RR.z[pid]
#                print "x0", x0, y0, x1, y1
#                print RR.z
#                print 'rotated[pid]', rotated[pid]
#                print 'rotated[pid]', rotated.sum(1)[pid]
                ax2.plot([x0,x1],[y0,y1])

            # Get sum for correct id.
            self.Sinogram[i] = rotated.sum(1)[pid]
            if verbose:
                print "col", self.Name, i, self.Angle[i], entry.side, rotated.sum(1)[pid], self.Angle[i]#, pid, RR.z[pid]
                print 'rotated[pid]', pid, rotated[pid]

            i=i+1

        if plot:
            plt.show()

        return

    def TotalDensity(self):
        """ Calculate total density for this camera.
        Basically use the sum of all lineintegrals"""
        return np.sum(self.Sinogram)

    def WriteLineIntegralToFile(self, fn):
        fp=open(fn,'w')
        # Write location
        fp.write('{0:s}\n'.format(self.Name))
        fp.write('{0:f} {1:f} {2:f}\n'.format(self.x,self.y,self.z))
        fp.write('{0:d} {1:f} {2:f}\n'.format(self.Nrays,self.theta1,self.theta2))
        # Write line integral
        for Ang,SG in zip(self.Angle,self.Sinogram):
            fp.write('{0:f} {1:f}\n'.format(Ang,SG))
        fp.close()
        return


    def plot(self, ax, step=1):
        i=0
        ax.text(self.PlotProp.text_x, self.PlotProp.text_y, self.Name)
        while i<self.Nrays:
          x0 = self.Rays_x[0,i]
          x1 = self.Rays_x[1,i]
          y0 = self.Rays_y[0,i]
          y1 = self.Rays_y[1,i]
          print "TomoCamera x0,x1,y0,y1",x0,x1,y0,y1
          color = self.PlotProp.color
          if i==0:
              color='k'
          ax.plot([x0,x1],[y0,y1], color=color)
          i=i+step
        return
