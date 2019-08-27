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
import numpy as np
import UVspec
import matplotlib.pyplot as plt
from scipy import ndimage

class Image:

    def __init__(self):
        self.Plotxsize=10
        self.Plotysize=7
        self.AspectRatio = 1.0
        self.title='Title'
        self.xlabel='X-axis'
        self.ylabel='Y-axis'
        self.cblabel='CB-label'
        self.fsize=20
        self.vmin=-999
        self.vmax=-999
        self.HideAxis=False
        return

    def ReadMYSTICradspcfile(self, fn, flip_y=False, Filter=None, RemoveGhosts=False,
                             AddRows=False,  timestart='005'):
        """
        Input is uvspec output file name as provided by MYSTIC.
        See libRadtran documentation for details, specifically
        mc_basename option.

        The following variables are set:

        rad: radiance
        std: standard deviation of Rad.

        The latter is only set if mc_std is included in
        uvspec input file.
        """
        RAD = np.zeros((self.ny,self.nx,self.nrgb))
        STD = np.zeros((self.ny,self.nx,self.nrgb))
        self.fn=fn
        f  = open(fn,'r')
        ir = 0
        it = 0
        for line in f:
            ls = line.split()
            ix = int(ls[1])
            iy = int(ls[2])
            RAD[iy,ix,ir] = float(ls[4])
            if  len(ls) > 5:
                STD[iy,ix,ir] = float(ls[5])

            if self.nrgb == 3 and it >= self.ny*self.nx:
                ir = ir + 1
                it = 0
            else:
                it = it + 1
        f.close()
        if self.nrgb==1:
            self.rad = RAD[:,:,0]
            self.std = STD[:,:,0]
            if flip_y:
                self.rad=self.rad[::-1,:]
                self.std=self.std[::-1,:]
        else:
            self.rad = RAD
            self.std = STD

        if Filter != None:
            self.rad = ndimage.gaussian_filter(self.rad, sigma=0.75)

        if RemoveGhosts:
            print "timestart", timestart
            if timestart=='001':
                self.rad[0:3, 197:203] = self.rad[0:3, 107:113]
            elif timestart=='002':
                self.rad[0:3, 204:222] = self.rad[0:3, 104:122]
            elif timestart=='003':
                self.rad[0:3, 210:225] = self.rad[0:3, 110:125]
            elif timestart=='004':
                self.rad[0:3, 220:240] = self.rad[0:3, 120:140]
            elif timestart=='005':
                self.rad[0:3, 229:239] = self.rad[0:3, 210:220]
            elif timestart=='006':
                self.rad[0:3, 230:245] = self.rad[0:3, 130:145]
            #            self.rad[0:5, 219:249] = self.rad[0:5, 190:220]
            else:
                1;
            pass

        if AddRows:
            print "self.ny,self.nx,self.nrgb", self.ny,self.nx,self.nrgb
            tmpRAD = np.zeros((self.nx,self.nx))+80
            tmpSTD = np.zeros((self.nx,self.nx))+0.8
            print self.rad.shape, tmpRAD.shape
            tmpRAD[150:150+self.ny,:] = self.rad[0:self.ny,:]
            tmpSTD[150:150+self.ny,:] = self.std[0:self.ny,:]
            self.rad = tmpRAD
            self.std = tmpSTD

    def MYSTIC_statistics(self, std_limit=0.01, verbose=True):
        indx=np.where(self.std/self.rad< std_limit )
        npix=len(indx[0])
        ntot=self.ny*self.nx
        if verbose:
            print 'Total number of pixels:', ntot
            print 'Number of pixels with std< {:f}: {:d}'.format(std_limit, npix)
            print 'Percentage: {:f}'.format(100*float(npix)/ntot)

    def GetUVSPECImageInputVals(self, fi):
        """
        Input is uvspec input file name.

        The following variables are set

        nx,ny,xpos,ypos,altitude,nrgb

        Their definitions are given in the libRadtran user's guide
        for the options mc_sample_grid and mc_sensorposition.

        Note that this function does not work properly if comments
        included on the same line as the above options, in the
        uvspec input file given as input.
         """

        vals = UVspec.get_vals(fi,'mc_sample_grid')
        if len(vals) == 2:
            self.nx = int(vals[0])
            self.ny = int(vals[1])
            self.AspectRatio = float(self.ny)/float(self.nx)
        else:
            print "Input file ",fi, " did not contain mc_sample_grid"
            print "This is not a panorama input file. Exiting."
            exit(0)
        vals = UVspec.get_vals(fi,'mc_panorama_view')
        if len(vals) == 4:
            self.phi1 = float(vals[0])
            self.phi2 = float(vals[1])
            self.theta1 = float(vals[2])
            self.theta2 = float(vals[3])
            self.theta1_org = float(vals[2])
            self.theta2_org = float(vals[3])
        else:
            print "Input file ",fi, " did not contain mc_panorama_view"
            print "This is not a panorama input file. Exiting."
            exit(0)

        vals = UVspec.get_vals(fi,'umu')
        if len(vals) == 1:
            self.umu = float(vals[0])

            # Shift theta to align with viewing direction.
            self.thetacenter = 0.5*(self.theta1+self.theta2)
            self.umudeg = np.rad2deg(np.arccos(self.umu))
            self.theta1 = self.theta1-self.thetacenter + self.umudeg
            self.theta2 = self.theta2-self.thetacenter + self.umudeg

        vals = UVspec.get_vals(fi,'mc_panorama_alignment')
        if len(vals)== 1:
            self.mc_panorama_alignment=True
        else:
            self.mc_panorama_alignment=False
            print "Input file ",fi, " with camera not aligned"
            print "Sure this is correct?"


        vals = UVspec.get_vals(fi,'mc_sensorposition')
        if len(vals) == 3:
            self.xpos     = float(vals[0])
            self.ypos     = float(vals[1])
            self.altitude = float(vals[2])


        vals = UVspec.get_vals(fi,'source')
        if vals[0] == 'thermal':
            self.thermal=True
        else:
            self.thermal=False

        vals = UVspec.get_vals(fi,'output')
        self.nrgb = 1

        if len(vals) == 1:
            if vals[0] == 'rgb':
                self.nrgb = 3

    def PlotRaw(self, Annotate=False, Extent=None,std=False):
        import matplotlib.pyplot as plt

        # Make pixel size the same in x- and y-directions.
        self.Plotysize=self.Plotxsize*self.AspectRatio
        fig = plt.figure(figsize=(self.Plotxsize,self.Plotysize))
        if std==True:
            tmpdat = self.STD
            print "Mean STD", np.mean(tmpdat), np.mean(self.STD/self.rad)

        else:
            tmpdat = self.rad

        if Extent==None:
            if self.vmin==-999 and self.vmax==-999:
                im = plt.imshow(tmpdat,origin='lower',aspect='auto')
            else:
                im = plt.imshow(tmpdat,origin='lower',aspect='auto',vmin=self.vmin, vmax=self.vmax)
        else:
            if self.vmin==-999 and self.vmax==-999:
                im = plt.imshow(tmpdat,origin='lower',aspect='auto',extent=Extent)
            else:
                im = plt.imshow(tmpdat,origin='lower',aspect='auto',extent=Extent,vmin=self.vmin, vmax=self.vmax)
        if Annotate:
            from matplotlib.legend import rcParams
            rcParams.update({'font.size':self.fsize})
            fig.subplots_adjust(top=0.94, right=0.96,left=0.11,bottom=0.11)
            cb=plt.colorbar()
            cb.set_label(self.cblabel)
            plt.xlabel(self.xlabel, fontsize = self.fsize)
            plt.ylabel(self.ylabel, fontsize = self.fsize)
            plt.title(self.title, fontsize = self.fsize)
        if self.HideAxis:
            frame1 = plt.gca()
            frame1.axes.get_xaxis().set_visible(False)
            frame1.axes.get_yaxis().set_visible(False)

        plt.show()

    def Plot310330tau(self, Annotate=False, Extent=None,std=False):

        # Make pixel size the same in x- and y-directions.
        self.Plotysize=self.Plotxsize*self.AspectRatio
        fig = plt.figure(figsize=(4*self.Plotxsize,self.Plotysize))
        if std==True:
            tmpdat = self.STD
            print "Mean STD", np.mean(tmpdat), np.mean(self.STD/self.rad)

        else:
            tmpdat = self.rad

        print "Extent", Extent
        ax = fig.add_subplot(131)
        self.cblabel='Radiance at 310 nm (mW m$^{-2}$ nm$^{-1}$)'
        self.PlotSingle(ax, fig, self.rad, Annotate=Annotate, Extent=Extent)
        ax = fig.add_subplot(132)
        self.cblabel='Radiance at 330 nm (mW m$^{-2}$ nm$^{-1}$)'
        self.PlotSingle(ax, fig, self.radBG, Annotate=Annotate, Extent=Extent)
        ax = fig.add_subplot(133)
        self.cblabel='Apparant absorbance'
        self.PlotSingle(ax, fig, self.tau, Annotate=Annotate, Extent=Extent)

        plt.show()

    def PlotSingle(self, ax, fig, dat, Extent=None, Annotate=False):
        if Extent==None:
            if self.vmin==-999 and self.vmax==-999:
                im = ax.imshow(dat,origin='lower',aspect='auto')
            else:
                im = ax.imshow(dat,origin='lower',aspect='auto',vmin=self.vmin, vmax=self.vmax)
        else:
            if self.vmin==-999 and self.vmax==-999:
                im = ax.imshow(dat,origin='lower',aspect='auto',extent=Extent)
            else:
                im = ax.imshow(dat,origin='lower',aspect='auto',extent=Extent,vmin=self.vmin, vmax=self.vmax)
        if Annotate:
            from matplotlib.legend import rcParams
            rcParams.update({'font.size':self.fsize})
            fig.subplots_adjust(top=0.94, right=0.96,left=0.05,bottom=0.11)
            cb=plt.colorbar(im)
            cb.set_label(self.cblabel)
            plt.xlabel(self.xlabel, fontsize = self.fsize)
            plt.ylabel(self.ylabel, fontsize = self.fsize)
            plt.title(self.title, fontsize = self.fsize)
        if self.HideAxis:
            frame1 = plt.gca()
            frame1.axes.get_xaxis().set_visible(False)
            frame1.axes.get_yaxis().set_visible(False)
        return

    def LineIntegrals(self, Camfile, Case, Slice=0):
        print self.fn
        ixc=Slice
        fp = open(Camfile,'w')
        fp.write("{:s}\n".format( Case))
        dtheta=np.abs(self.theta2-self.theta1)/self.ny
        if 'A' in Case:
            offset=90
        else:
            offset=270
        fp.write("{:f} {:f} {:f}\n".format( self.ypos, self.xpos, self.altitude))
#        fp.write("{:f} {:f} {:f}\n".format( self.xpos, self.ypos, self.altitude))
        if 'A' in Case:
            fp.write("{:d} {:f} {:f}\n".format( self.ny, self.theta1-offset, self.theta2-offset-dtheta)) #, self.nx,
        else:
            fp.write("{:d} {:f} {:f}\n".format( self.ny, offset-self.theta2-dtheta, offset-self.theta1)) #, self.nx,
        print dtheta
        if 'A' in Case:
            iy=0
            while iy<self.ny:
                fp.write("{:f} {:e}\n".format(self.theta1-offset+iy*dtheta, self.tau[iy,ixc]))
                #            fp.write("{:f} {:e}\n".format(offset-self.theta1+iy*dtheta, self.tau[iy,ixc]))
                iy=iy+1
        else:
            iy=self.ny-1
            while iy>=0:
                fp.write("{:f} {:e}\n".format(offset-self.theta1-iy*dtheta, self.tau[iy,ixc]))
                iy=iy-1
        fp.close()
        return


###########################################################

if __name__ == "__main__":
    verbose = True

    ART=False #True#
    Annotate=True#False #
    CalcTau=True #False#
    LineIntegrals=True
    Case='' #'_B_5.0_gabba'#'_D'#'_A'#
    wavelengthA=310
    wavelengthB=330
    plottype='310330tau' # ''
    sza = 40
    std=False #True #
    Case ='LocC'#'LocAAerosol_A0_AP0_R01' #'LocAAerosol_A1_KSTGeo01' #  'Aerosol_None'#'LocAAerosol_A1_KSTGeo02'#'AerosolDefault'# 'Aerosol_None'#aerosol_default'


    # Value of extraterrestrial spectrum from
    # libRadtran/data/solar_flux/atlas_plus_modtran
    # used to normalize radiances
    if wavelengthA==310:
        nf = 346.590*np.cos(np.pi*sza/180)
    elif wavelengthA==550:
        nf = 1885.9515*np.cos(np.pi*sza/180)
    nf=1

#    filebase='/home/aky/NILU/xnilu_wrk/aky/Projects/COMTESSA/Experiments/testJG/' #'Guallatiri/'#TEST/'#AerosolDefault/'#Ellipsoid/'#VerticalPlume/'#ART/' # tmp
#    filebase='/home/aky/Projects/COMTESSA/Experiments/GuallatiriBlueSky/' #Cell/' #ART/' #Ellipsoid/'# tmp/'# #TEST/'#
    filebase='/home/aky/Projects/COMTESSA/Experiments/Ellipsoid_3Cams/' #
#    filebase='/home/aky/NILU/xnilu_wrk/aky/Projects/COMTESSA/bin/Examples/HorizontalPlume/' #Sphere/'
#    uvspecinputfile = filebase+'Wavelength_{:5.1f}{:2s}.inp'.format(wavelength,Case)
#    uvspecoutputfile = filebase+'Wavelength_{:5.1f}{:2s}_mc.rad.spc'.format(wavelength,Case)
#    uvinp = 'Wavelength_{:5.1f}-Sza_{:4.1f}-{:s}.inp'.format(wavelengthA,sza,Case)
#    uvout = 'Wavelength_{:5.1f}-Sza_{:4.1f}-{:s}_mc.rad.spc'.format(wavelengthA,sza,Case)
#    uvinp = 'CamW{:3.0f}{:s}uvspec.inp'.format(wavelengthA,Case)
#    uvout = 'CamW{:3.0f}{:s}uvspec.out'.format(wavelengthA,Case)
    uvinp = 'uvspecCamW{:3.0f}{:s}.inp'.format(wavelengthA,Case)
    uvout = 'uvspecCamW{:3.0f}{:s}.out'.format(wavelengthA,Case)
#    uvinp = 'uvspecW{:3.0f}{:s}.inp'.format(wavelengthA,Case)
#    uvout = 'uvspecW{:3.0f}{:s}.out'.format(wavelengthA,Case)
#    uvinp = 'Cam1uvspec.inp'
#    uvout = 'Cam1uvspec.out'
#    uvinp = 'uvspecW310LocA.inp'
#    uvout = 'uvspecW310LocA.out'

    uvspecinputfile = filebase+uvinp
    uvspecoutputfile = filebase+uvout

#    sza = UVspec.get_vals(uvspecinputfile,'sza')
#    sza = float(sza[0])
#    wavelength = UVspec.get_vals(uvspecinputfile,'wavelength')
#    wavelength = float(wavelength[0])

#uvspecoutputfile = './test_cloud/mc.rad.spc'
#    print uvspecoutputfile
    ImgAM = Image()
    if ART:
        ImgAM.Plotxsize=20
        ImgAM.Plotysize=14
        ImgAM.vmin=3.6
        ImgAM.vmax=10.0
        ImgAM.HideAxis=True
        Annotate=False

#    ImgAM.title='Wavelength {:5.1f} nm, solar zenith angle {:4.1f}'.format(wavelength,sza)
    ImgAM.title='Solar zenith angle {:4.1f}, {:s}'.format(sza,Case)
    ImgAM.xlabel='Azimuth angle ($^{\circ}$)'
    ImgAM.ylabel='Polar angle ($^{\circ}$)'
#    Img.cblabel='Normalized radiance'
    ImgAM.cblabel='Radiance (mW m$^{-2}$ nm$^{-1}$)'

    ImgAM.GetUVSPECImageInputVals(uvspecinputfile)
    extent=[ImgAM.phi1,ImgAM.phi2,ImgAM.theta1,ImgAM.theta2]

    ImgAM.ReadMYSTICradspcfile(uvspecoutputfile)
#    Img.rad = Img.rad/nf

    if CalcTau:

        # The naming of the various Image objects corresponds to the indices used
        # by Lubcke et al (2013) in Eq. (3)
        # AM: Wavelength A influenced by SO2, viewing plume
        # A0: Wavelength A influenced by SO2, background image without SO2 plume
        # BM: Wavelength B not influenced by SO2, viewing plume
        # B0: Wavelength B notinfluenced by SO2, background image without SO2 plume


        ImgA0= Image()
        uvspecinputfileA0=uvspecinputfile.replace(Case,Case+'BG')
        uvspecoutputfileA0=uvspecinputfileA0.replace('.inp','.out')
#        uvspecoutputfileA0=uvspecoutputfile.replace('_mc.rad.spc','-BG_mc.rad.spc')
        if verbose:
            print uvspecoutputfileA0
        ImgA0.GetUVSPECImageInputVals(uvspecinputfileA0)
        ImgA0.ReadMYSTICradspcfile(uvspecoutputfileA0)

        ImgBM= Image()
        uvspecinputfileBM=uvspecinputfile.replace(str(wavelengthA),str(wavelengthB))
        uvspecoutputfileBM=uvspecoutputfile.replace(str(wavelengthA),str(wavelengthB))
        if verbose:
            print uvspecoutputfileBM
        ImgBM.GetUVSPECImageInputVals(uvspecinputfileBM)
        ImgBM.ReadMYSTICradspcfile(uvspecoutputfileBM)

        ImgB0= Image()
        uvspecinputfileB0=uvspecinputfileBM.replace(Case,Case+'BG')
        uvspecoutputfileB0=uvspecinputfileB0.replace('.inp','.out')
#        uvspecoutputfileB0=uvspecoutputfileBM.replace('_mc.rad.spc','-BG_mc.rad.spc')
        if verbose:
            print uvspecoutputfileB0
        ImgB0.GetUVSPECImageInputVals(uvspecinputfileB0)
        ImgB0.ReadMYSTICradspcfile(uvspecoutputfileB0)

        ImgAM.radBG = ImgBM.rad
#        ImgAM.radorg = ImgAM.rad[:]
#        ImgAM.rad = np.log((ImgBM.rad/ImgAM.rad)*(ImgA0.rad/ImgB0.rad))
        ImgAM.tau = np.log((ImgBM.rad/ImgAM.rad)*(ImgA0.rad/ImgB0.rad))
#       ImgAM.rad = np.log((1./ImgAM.rad)*(ImgA0.rad/1.))
        ImgAM.cblabel='Apparant absorbance'
        ImgAM.title='Solar zenith angle {:4.1f}, {:s}'.format(sza,Case)


    if LineIntegrals:
        ImgAM.LineIntegrals(Case)

    if plottype=='310330tau':
        ImgAM.Plot310330tau(Annotate=Annotate,Extent=extent,std=std)
    else:
        ImgAM.PlotRaw(Annotate=Annotate,Extent=extent,std=std)
