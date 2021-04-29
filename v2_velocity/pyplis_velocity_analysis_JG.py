import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import argparse
from matplotlib import rcParams
import shutil
CLEAR = True

# -*- coding: utf-8 -*-
"""
Modified heavily from pyplis example script no. 3 - Plume background analysis

Arve Kylling, NILU.

This script .
"""
import numpy as np
import pyplis
import matplotlib.pyplot as plt

import datetime
import Image

from glob import glob

def Centerline(xs, ys, col2d):
    # Equation from Dosio and de Arellano (2006).
    Centerline = np.zeros(xs.shape[0])
    for i in (list(range(xs.shape[0]))):
        colsum= np.sum(col2d[:,i])
        if colsum > 0.0:
            Centerline[i] = np.sum(ys*col2d[:,i])/colsum

    # This is shorter, but gives problems for colsum==0
    #    Centerline = np.array([np.sum(ys*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))])

    return Centerline, np.nanmean(Centerline)

def AbsoluteDispersion(xs, ys, col2d, mean):
    # i=0
    Dispersion = np.zeros(xs.shape[0])
    # for x in xs:
    #     j=0
    #     sum1=0
    #     sum2=0
    #     for y in ys:
    #         zp=y-mean
    #         sum1=sum1+zp*zp*col2d[j,i]
    #         sum2=sum2+col2d[j,i]
    #         j=j+1
    #     Dispersion[i]=np.sqrt(sum1/sum2)
    #     i=i+1

    for i in (list(range(xs.shape[0]))):
        colsum= np.sum(col2d[:,i])
        if colsum > 0.0:
            Dispersion[i] = np.sqrt(np.sum((ys-mean)*(ys-mean)*col2d[:,i])/colsum)

    # This is shorter, but gives problems for colsum==0
    #Dispersion = np.sqrt(np.array([np.sum((ys-mean)*(ys-mean)*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))]))

    return Dispersion

def RelativeDispersion(xs, ys, col2d, yc):
    # i=0
    Dispersion = np.zeros(xs.shape[0])
    # for x in xs:
    #     j=0
    #     sum1=0
    #     sum2=0
    #     for y in ys:
    #         zp=y-yc[i]
    #         sum1=sum1+zp*zp*col2d[j,i]
    #         sum2=sum2+col2d[j,i]
    #         j=j+1
    #     Dispersion[i]=np.sqrt(sum1/sum2)
    #     i=i+1
    # print ys.shape, yc.shape

    for i in (list(range(xs.shape[0]))):
        colsum= np.sum(col2d[:,i])
        if colsum > 0.0:
            Dispersion[i] = np.sqrt(np.sum((ys-yc[i])*(ys-yc[i])*col2d[:,i])/colsum)

    # This is shorter, but gives problems for colsum==0
    #    Dispersion = np.sqrt(np.array([np.sum((ys-yc[i])*(ys-yc[i])*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))]))

    return Dispersion

def Skewness(xs, ys, col2d, mean):
    #     i=0
    Skewness  = np.zeros(xs.shape[0])
    #     for x in xs:
    #         j=0
    #         sum1=0
    #         sum2=0
    #         for y in ys:
    #             zp=y-mean
    #             sum1=sum1+zp*zp*zp*col2d[j,i]
    #             sum2=sum2+col2d[j,i]
    #             j=j+1
    #         Skewness[i]=sum1/sum2
    # #        Skewness[i]=cubic_root(sum1/sum2)
    #         i=i+1

    for i in (list(range(xs.shape[0]))):
        colsum= np.sum(col2d[:,i])
        if colsum > 0.0:
            Dispersion = np.sqrt(np.sum((ys-mean)*(ys-mean)*col2d[:,i])/colsum)
            Skewness[i] = (np.sum((ys-mean)*(ys-mean)*(ys-mean)*col2d[:,i])/colsum)/np.power(Dispersion, 1.5)

    # This is shorter, but gives problems for colsum==0
    #Dispersion = np.sqrt(np.array([np.sum((ys-mean)*(ys-mean)*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))]))
    #Skewness = np.array([np.sum((ys-mean)*(ys-mean)*(ys-mean)*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))])/np.power(Dispersion, 1.5)

    return Skewness

def AbsoluteSkewness(xs, ys, col2d, mean):
    #    i=0
    Skewness  = np.zeros(xs.shape[0])
    #     for x in xs:
    #         j=0
    #         sum1=0
    #         sum2=0
    #         for y in ys:
    #             zp=y-mean
    #             sum1=sum1+zp*zp*zp*col2d[j,i]
    #             sum2=sum2+col2d[j,i]
    #             j=j+1
    #         Skewness[i]=sum1/sum2
    # #        Skewness[i]=cubic_root(sum1/sum2)
    #         i=i+1

    for i in (list(range(xs.shape[0]))):
        colsum= np.sum(col2d[:,i])
        if colsum > 0.0:
            Skewness[i] =np.sum((ys-mean)*(ys-mean)*(ys-mean)*col2d[:,i])/colsum

    # This is shorter, but gives problems for colsum==0
    #    Skewness = np.array([np.sum((ys-mean)*(ys-mean)*(ys-mean)*col2d[:,i])/np.sum(col2d[:,i]) for i in (list(range(xs.shape[0])))])

    return Skewness


class Img(pyplis.Img):
    def __myinit__(self):
        return

    def statistics(self, axis=0, verbose=False):
        xp = np.arange(0,self.meta["pix_width"])
        yp = np.arange(0,self.meta["pix_heigth"])
        yp=yp[::-1] # Reverse since image pixels go the other way around
        if verbose:
            print("statistics", xp, yp, self.img.shape)

        self.Centerline, self.mean = Centerline(xp, yp, self.img)
        self.Dispersion = AbsoluteDispersion(xp, yp, self.img, self.mean)
        self.RelativeDispersion = RelativeDispersion(xp, yp, self.img, self.Centerline)
        self.AbsoluteSkewness = self.Centerline*0.0
        self.Skewness = self.Centerline*0.0
        i=0
        for AbsSkew, Skew, Disp3 in zip(AbsoluteSkewness(xp, yp, self.img, self.mean), Skewness(xp, yp, self.img, self.mean), np.power(self.Dispersion,3)):
            if Disp3>0.0:
                self.AbsoluteSkewness[i] = AbsSkew/Disp3
                self.Skewness[i] = Skew/Disp3
            i=i+1

        # Simpler, but risk dividing by zero.
        #        self.AbsoluteSkewness = AbsoluteSkewness(xp, yp, self.img, self.mean)/np.power(self.Dispersion,3)
        #        self.Skewness = Skewness(xp, yp, self.img, self.mean)/np.power(self.Dispersion,3)
        return

    def diff(self, imgB):
        diff=self
        diff.img = self.img-imgB.img
        diff.Centerline = self.Centerline - imgB.Centerline
        diff.Dispersion = self.Dispersion - imgB.Dispersion
        diff.RelativeDispersion = self.RelativeDispersion - imgB.RelativeDispersion
        diff.AbsoluteSkewness = self.AbsoluteSkewness - imgB.AbsoluteSkewness
        diff.Skewness = self.Skewness - imgB.Skewness

        return diff


def load_MYSTIC_image(filebase, wavelengthA,loc,timestartstr,experimenttype,
                      wavelengthB, filebaseBG, timestampBG,experimenttypeBG,
                      start_acq=datetime.datetime(2016, 10, 10, 13, 15, 12),
                      print_MYSTIC_statistics=False, flip_y=False, Filter=None,
                      AddRows=False, verbose=True):
    """Load images

    """

    uvinp = 'uvspecCamW{:s}{:s}{:s}{:s}.inp'.format(wavelengthA,loc,timestartstr,experimenttype)
    uvout = 'uvspecCamW{:s}{:s}{:s}{:s}.out'.format(wavelengthA,loc,timestartstr,experimenttype)
    uvspecinputfile = filebase+uvinp
    uvspecoutputfile = filebase+uvout
    if verbose:
        #        print("uvspecinputfile", uvspecinputfile)
        #        print("uvspecoutputfile", uvspecoutputfile)
        print('uvspecinputfile,{:s}'.format(uvspecinputfile))
        print('uvspecoutputfile,{:s}'.format(uvspecoutputfile))
    ImgAM = Image.Image()
    ImgAM.GetUVSPECImageInputVals(uvspecinputfile)
    ImgAM.ReadMYSTICradspcfile(uvspecoutputfile, flip_y=flip_y, Filter=Filter,
                               RemoveGhosts=True, AddRows=AddRows, timestart=timestartstr)
    if print_MYSTIC_statistics:
        ImgAM.MYSTIC_statistics(verbose=True)

    MYSTIC_AM = Img(start_acq = start_acq)
    MYSTIC_AM.img=ImgAM.rad
    MYSTIC_AM.meta["bit_depth"]=12  # Just to get show_histogram working.
    MYSTIC_AM.meta["pix_heigth"]=ImgAM.rad.shape[0]
    MYSTIC_AM.meta["pix_width"]=ImgAM.rad.shape[1]
    if verbose:
        print("pix_heigth, pix_width", MYSTIC_AM.meta["pix_heigth"],
              MYSTIC_AM.meta["pix_width"])

    ImgBM= Image.Image()
    uvspecinputfileBM=uvspecinputfile.replace(str(wavelengthA),str(wavelengthB))
    uvspecoutputfileBM=uvspecoutputfile.replace(str(wavelengthA),str(wavelengthB))
    ImgBM.GetUVSPECImageInputVals(uvspecinputfileBM)
    ImgBM.ReadMYSTICradspcfile(uvspecoutputfileBM, flip_y=flip_y, Filter=Filter,
                               RemoveGhosts=True, AddRows=AddRows, timestart=timestartstr)
    if print_MYSTIC_statistics:
        ImgBM.MYSTIC_statistics(verbose=True)
    MYSTIC_BM = Img(start_acq = start_acq)
    MYSTIC_BM.img=ImgBM.rad
    MYSTIC_BM.meta["bit_depth"]=12  # Just to get show_histogram working.
    MYSTIC_BM.meta["pix_heigth"]=ImgBM.rad.shape[0]
    MYSTIC_BM.meta["pix_width"]=ImgBM.rad.shape[1]

    # Background files
    uvspecinputfileA0 = filebaseBG+uvinp
    uvspecoutputfileA0 = filebaseBG+uvout
    uvspecinputfileA0=uvspecinputfileA0.replace(loc+timestartstr+experimenttype,loc+timestampBG+experimenttypeBG)
    uvspecoutputfileA0=uvspecinputfileA0.replace('.inp','.out')
    if verbose:
        #        print("uvspecinputfileA0", uvspecinputfileA0)
        #        print("uvspecoutputfileA0", uvspecoutputfileA0)
        print('uvspecinputfileA0,{:s}'.format(uvspecinputfileA0))
        print('uvspecoutputfileA0,{:s}'.format(uvspecoutputfileA0))
    ImgA0= Image.Image()
    ImgA0.GetUVSPECImageInputVals(uvspecinputfileA0)
    ImgA0.ReadMYSTICradspcfile(uvspecoutputfileA0, flip_y=flip_y, Filter=Filter,
                               RemoveGhosts=True, AddRows=AddRows, timestart=timestartstr)
    if print_MYSTIC_statistics:
        ImgA0.MYSTIC_statistics(verbose=True)
    MYSTIC_A0 = Img(start_acq = start_acq)
    MYSTIC_A0.img=ImgA0.rad
    MYSTIC_A0.meta["bit_depth"]=12  # Just to get show_histogram working.
    MYSTIC_A0.meta["pix_heigth"]=ImgA0.rad.shape[0]
    MYSTIC_A0.meta["phi1"]=ImgA0.phi1
    MYSTIC_A0.meta["phi2"]=ImgA0.phi2

    ImgB0= Image.Image()
    uvspecinputfileB0 = filebaseBG+uvinp
    uvspecoutputfileB0 = filebaseBG+uvout
    if verbose:
        #        print("uvspecinputfileB0", uvspecinputfileB0)
        #        print("uvspecoutputfileB0", uvspecoutputfileB0)
        print('uvspecinputfileB0,{:s}'.format(uvspecinputfileB0))
        print('uvspecoutputfileB0,{:s}'.format(uvspecoutputfileB0))
    uvspecinputfileB0=uvspecinputfileB0.replace(loc+timestartstr+experimenttype,loc+timestampBG+experimenttypeBG)
    uvspecinputfileB0=uvspecinputfileB0.replace(str(wavelengthA),str(wavelengthB))
    uvspecoutputfileB0=uvspecinputfileB0.replace('.inp','.out')
    uvspecoutputfileB0=uvspecoutputfileB0.replace(str(wavelengthA),str(wavelengthB))
    ImgB0.GetUVSPECImageInputVals(uvspecinputfileB0)
    ImgB0.ReadMYSTICradspcfile(uvspecoutputfileB0, flip_y=flip_y, Filter=Filter,
                               RemoveGhosts=True, AddRows=AddRows, timestart=timestartstr)
    if print_MYSTIC_statistics:
        ImgB0.MYSTIC_statistics(verbose=True)
    MYSTIC_B0 = Img(start_acq = start_acq)
    MYSTIC_B0.img=ImgB0.rad
    MYSTIC_B0.meta["bit_depth"]=12  # Just to get show_histogram working.
    MYSTIC_A0.meta["pix_heigth"]=ImgA0.rad.shape[0]
    MYSTIC_A0.meta["pix_width"]=ImgA0.rad.shape[1]

    return MYSTIC_AM, MYSTIC_A0, MYSTIC_BM, MYSTIC_B0


def load_MYSTIC_images(filebase, wavelengthA,
                       loc, timestart, timeend, timestep, experimenttype,
                       wavelengthB, filebaseBG,
                       timestampBG, experimenttypeBG, flip_y=False, Filter=None,
                       AddRows=False):

    ImagesAA=[]

    if int(timestep)>0:
        timestamps = list(range(int(timestart), int(timeend)+1, int(timestep)))
    else:
        timestamps = list(range(int(timestart), int(timeend)-1, int(timestep)))

    k=0
    #    print "timestamps", timestamps
    for timestamp in timestamps:
        timestampstr   = '{:0003d}'.format(int(timestamp))

        MYSTIC_AM, MYSTIC_A0, MYSTIC_BM, MYSTIC_B0 = load_MYSTIC_image(filebase, wavelengthA,
                                                                       loc, timestampstr, experimenttype,
                                                                       wavelengthB, filebaseBG,
                                                                       timestampBG, experimenttypeBG,
                                                                       flip_y=flip_y, Filter=Filter,
                                                                       AddRows=AddRows)
#                                                                        timestampBG, experimenttypeBG, start_acq=times[k])
        tau_A = MYSTIC_AM.to_tau(MYSTIC_A0) # Note that Jonas definition of tau differs
                                            # with a minus sign from Lubcke et al. 2013 definition.
        tau_B = MYSTIC_BM.to_tau(MYSTIC_B0)
        tau_A.img = tau_A.img-tau_B.img
        # Zero all values below threshold
        tau_A.threshold=0.03 #0.00000001 #
        tau_A.set_val_below_thresh(val=0, threshold=tau_A.threshold)
        tau_A.statistics()
        ImagesAA.append(tau_A)
        k=k+1

    return ImagesAA

def downscale_images_and_pcs_lines(images, pcs_lines, pyrlevels_down=3):
    imgs_downscaled = []
    lines_downscaled = []
    for img in images:
        imgs_downscaled.append(img.duplicate().pyr_down(pyrlevels_down))
    for line in pcs_lines:
        lines_downscaled.append(line.convert(to_pyrlevel=3))
    return imgs_downscaled, lines_downscaled

### SCRIPT MAIN FUNCTION
if __name__=="__main__":
    import os

    FONTSIZE=20

    parser = argparse.ArgumentParser(description='Velocity analysis of simulated images.')

    parser.add_argument('--COL_NUM1', type=int, default=200,
                        help='Number of column 1')
    parser.add_argument('--COL_NUM2', type=int, default=202,
                        help='Number of column 2')
    parser.add_argument('--COL_NUM3', type=int, default=204,
                        help='Number of column 3')
    parser.add_argument('--Camera', type=str,  default='LocA',
                        help='Camera location to be plotted')
    parser.add_argument("--timestart", type=str, help="Start time stamp of LES simulation", default='001')
    #    parser.add_argument("--timeend", type=str, help="End time stamp of LES simulation", default='010')
    parser.add_argument("--timestep", type=str, help="Time step in stamp units of LES simulation", default='1')

    args = parser.parse_args()
    COL_NUM1 = args.COL_NUM1
    COL_NUM2 = args.COL_NUM2
    COL_NUM3 = args.COL_NUM3
    loc             = args.Camera
    timestart   = args.timestart
    timestep    = args.timestep

    CLEAR_OUTPUT = True
    plt.close("all")

    PrintInput=False #True
    if PrintInput:
        print('Camera', loc)
        print('COL_NUM1', COL_NUM1)
        print('COL_NUM2', COL_NUM2)
        print('COL_NUM3', COL_NUM3)
        print('timestart', timestart)
        print('timestep', timestep)
#        exit()

    outdir = './output_new'+'/'+loc+'_'+str(COL_NUM2)+'/'
    txtfile = outdir+'velocities.txt'
    # output for results withoug multigauss analysis
    txtfile_noMG = outdir+'velocities_noMG.txt'


    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    ftxt = open(txtfile,'w')
    ftxt_noMG = open(txtfile_noMG,'w')

    if CLEAR_OUTPUT:

        for f in glob('{}*.png'.format(outdir)):
            os.remove(f)

    #    loc           = 'LocA'
    if 'LocV' in loc:
        FONTSIZE=14
        #    filebase      = './'
        #    filebaseBG    = './'
        #        filebase='/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena_HT_Velo/'
        #        filebaseBG='/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena_HT_Velo/'
        filebase='../images/'
        filebaseBG='../images/'
        experimenttype  = 'SUNSZA40WQBLC'
        experimenttypeBG  = 'SUNSZA40WQBLC_BG'
        timeend       = '010'
        # For LocV?? columns run from 0 to 99.
        #        COL_NUM1 = 20#5
        #        COL_NUM2 = 23 #60#0
        #        COL_NUM3 = 26 #60#0
        ROW_BOTTOM = 0+60# 150
        ROW_TOP    = 60 + 150 #110+150

    elif 'LocA' in loc:
        FONTSIZE=14
#        filebase='/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena_high_time_resolution/'
#        filebaseBG='/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena/'
        filebase='../images/'
        filebaseBG='../images/'
        experimenttype  = 'SUNSZA40WQBLC'
        experimenttypeBG  = 'SUNSZA40W_BG'
        timeend       = '008'
        # For LocA columns run from 0 to 399. First about 10 have  no signal.
        #        COL_NUM1 = 200#5
        #        COL_NUM2 = 203 #10 #60#0
        #        COL_NUM3 = 206 #10 #60#0
        ROW_BOTTOM = 0+20# 150
        ROW_TOP    = 20 + 64 #110+150

    DEPTH = COL_NUM3 - COL_NUM1
    wavelengthA   = '310'
    wavelengthB   = '330'

    timestampBG    = '005'
    experimenttypesensitivity  = ''
    experimenttypesensitivityBG  = ''
    Filter='Gaussian'#None#

    # Noise added to images (im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape))
    noiseamp = 0
    # time step between images
    DT=6.25/4 # 6.25 #s

    T0 = datetime.datetime(2017,1,1,0,0,0) # define arbitrary start time

    # Create an example plume intersection
    # Camera geometry
    #    DY=(3.828125+3.028125)/2.
    #    DZ=(3.903125+3.100000)/2.
    if 'LocV' in loc:
        DY=0.250 # Horizontal distance from camera to plume, from Fig 2., in km
        DZ=0.150 # Vertical distance from camera to plume, from Fig 2., in km
        dxhalfang = 11.5/2 # Size of half image in the horizontal, in degrees
        dzhalfang = 30/2    # Size of half image in the vertical, in degrees
        # 100/2 pixels in the horizontal gives horizontal pixel size
        npixelshor = 100
        npixelsver = 264
        npixelshorhalf = npixelshor/2
        if 'LocV01' in loc:   vmax=0.35
        elif 'LocV02' in loc:   vmax=0.15
        elif 'LocV03' in loc:   vmax=0.15

        # Size of half image in the vertical, 5. = (96-86)/2.
        #        ZS=np.sin(np.deg2rad(5))*DS
        # 88/2 pixels in the vertical gives vertical pixel size
        #        dz = ZS / 44.0 #=0.0063961804401032104 #km =6.396180440m
        # From uvspec input file: mc_panorama_view 157.000000 203.000000 86.000000 96.000000

    elif 'LocA' in loc:
        DY=0.250 # Horizontal distance from camera to plume, from Fig 2., in km
        DZ=0.025 # Vertical distance from camera to plume, from Fig 2., in km
        dxhalfang = 46/2  # Size of half image in the horizontal, in degrees
        dzhalfang = 10/2   # Size of half image in the vertical, in degrees
        npixelshor = 400
        npixelsver = 88
        npixelshorhalf = npixelshor/2
        vmax=1.2

    DS = np.sqrt(DZ*DZ+DY*DY)   # in km
    XS = np.sin(np.deg2rad(dxhalfang))*DS  # in km
    DX = (XS / npixelshorhalf) * 1000   # in m
    print("DX, ", DX, DY, DZ, DS, XS, npixelshorhalf)


    # Pixel size in x-direction without correcting for vertical  distance
    dx = DY * np.tan(np.deg2rad(dxhalfang)) *1000/npixelshorhalf
    # Correct for vertical distance
    dx = dx/np.cos(np.arctan(DZ/DY))
    DX=dx
    print("dx", dx,  1/np.cos(np.arctan(DZ/DY)))
#    exit()


    # Clear sky ROI
    BG_ROI = [0, 10, 0, 40]
    # Setup optical flow class (you may play with the parameters)

    ### DEFINE PCS LINES
    pcs1 = pyplis.LineOnImage(COL_NUM1, ROW_BOTTOM, COL_NUM1, ROW_TOP,
                             line_id="L", normal_orientation="right", #"left",#
                              color="#ff0000", linestyle=':')
    pcs2 = pyplis.LineOnImage(COL_NUM2, ROW_BOTTOM, COL_NUM2,  ROW_TOP,
                             line_id="M", normal_orientation="right", #"left",
                             color="#ffa500")
    pcs3 = pyplis.LineOnImage(COL_NUM3, ROW_BOTTOM, COL_NUM3,  ROW_TOP,
                             line_id="R", normal_orientation="right", #"left",
                             color="#ff0000", linestyle=':')

    ### LOAD IMAGES
    ImagesAA = load_MYSTIC_images(filebase, wavelengthA,
                                  loc, timestart, timeend, timestep, experimenttype,
                                  wavelengthB, filebaseBG,
                                  timestampBG, experimenttypeBG, flip_y=True,
                                  Filter=Filter,
                                  AddRows=False)
#                                  AddRows=True)

    times = [T0 + datetime.timedelta(x * DT / 86400.0) for x in range(len(ImagesAA))]

    ### LOAD SMALL IMAGES (NOT USED CURRENTLY)
    #    ImagesAA_small, (pcs1_small, pcs2_small) = downscale_images_and_pcs_lines(
    #            ImagesAA, (pcs1, pcs2))

    ### COMPUTE ICA TIMESERIES FOR BOTH LINES BASED ON N INPUT IMAGES
    # (includes plotting of AA images with lines in it)
    ts_pcs1 = []
    ts_pcs2 = []
    ts_pcs3 = []

    images = []

    for i, im in enumerate(ImagesAA):
        im.img = im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape)
        # noisy = f + 0.4 * f.std() * np.random.random(f.shape)
        #            print im
        rcParams.update({'font.size':FONTSIZE})
        kwargs={'vmin':0.0,
                'vmax':vmax}
        ax = im.show_img(**kwargs)

        F = plt.gcf()
        DefaultSize = F.get_size_inches()
        print( "Default size in inches", DefaultSize, im.img.shape)
        x1 = DefaultSize[0]; y1 = DefaultSize[1]
        x2 = im.img.shape[1]; y2 = im.img.shape[0]
        if y2 > x2:
            xn = 1.7*y1*(x2/y2)
            yn = y1
            left=0.0
            right=0.95
        else:
            xn = x1
            yn = x1*(y2/x2)
            left=0.05
            right=0.98
        F.set_size_inches( (xn, yn ))
        F.subplots_adjust(left=left, right=right,bottom=0.05,top=0.95)
        print( "New size in Inches", xn, yn, y2/x2)

        im.meta['start_acq'] = times[i]
        im.edit_log['is_aa'] = True
        #        ax.set_title("AA image {} (pyrlevel {})".format(i+1, im.pyrlevel))
        ax.set_title("AA image {}".format(i+1))

        pcs1.plot_line_on_grid(ax=ax)
        pcs2.plot_line_on_grid(ax=ax)
        pcs3.plot_line_on_grid(ax=ax)
        ax.legend(loc='upper right')
        images.append(im)
        rcParams.update({'font.size':FONTSIZE})
        ax.figure.savefig('{}aa_img{}_N{}.png'.format(outdir, i,
                          noiseamp))

        ts_pcs1.append(sum(pcs1.get_line_profile(im.img)))
        ts_pcs2.append(sum(pcs2.get_line_profile(im.img)))
        ts_pcs3.append(sum(pcs3.get_line_profile(im.img)))


    ### INTERPOLATE THE 2 PCS timeseries to better register the lag
    from scipy.interpolate import interp1d
    flow = pyplis.OptflowFarneback(disp_skip=10)

    line = pcs2
    line.set_rect_roi_rot(depth=DEPTH)

    props = pyplis.plumespeed.LocalPlumeProperties(line.line_id,
                                                   color=line.color)


    veff = {0 : [], # without multigauss fit
            1: []} # without multigauss fit
    veff_err = {0 : [], 1: []}
    mus= {0 : [], 1: []}



    for i in range(len(images)-1):
        flow.calc_flow(images[i], images[i+1])
        flow.roi_abs = line._roi_from_rot_rect()
        # JGLISS: the following function will do the histogram analysis of
        # the retrieved flow field and based on that will calculate the
        # average velocity in the rectangle around the line. Use
        # book "dir_multi_gauss" to activate / deactivate multigauss fit
        # retrieval in histogram analysis.
        for withMG in (0, 1):
            props.get_and_append_from_farneback(flow, line=line,
                                                dir_multi_gauss=withMG)

            veff1, veff_err1 = props.get_velocity(pix_dist_m=DX,
                                                  normal_vec=pcs2.normal_vector,
                                                  sigma_tol=flow.settings.hist_sigma_tol)
            print("veff1,  veff_err1", i, veff1,  veff_err1)
            fig, axd = plt.subplots(1,1,figsize=(16,10))
            (_,  mu,  sigma) = flow.plot_orientation_histo(pix_mask=None,
                                                           apply_fit=True,
                                                           ax=axd, color='red')
            plt.close()
            print('mu, sigma', mu, sigma);
            if i>0:
                veff[withMG].append(veff1)
                veff_err[withMG].append(veff_err1)
                mus[withMG].append(mu)

    print("veff_avg (NO MULTIGAUSS FIT)", np.nanmean(veff[0]), np.std(veff[0]))
    print("veff_avg (WITH MULTIGAUSS FIT)", np.nanmean(veff[1]), np.std(veff[1]))

    fig, ax = plt.subplots(1,1,figsize=(16,10))
    flow.plot(ax=ax)
    fig.savefig(outdir+'Farneback_example_output.png')


    fig = flow.plot_flow_histograms()
    fig.savefig(outdir+'Farneback_example_histograms.png')

    fig = props.plot()
    fig.savefig(outdir+'Farneback_plume_props_timeseries.png')

    ax = props.plot_velocities(pix_dist_m=DX,
                               normal_vec=line.normal_vector)

    ax.figure.savefig(outdir+'Farneback_velocities.png')
    # original index
    x = range(len(ts_pcs1))

    res_fac = 100
    last_idx = x[-1]
    num_itp = res_fac * last_idx
    xnew = np.linspace(0, last_idx, num=num_itp, endpoint=True)

    # the two timeseries (each 5 timestamps)
    ts_pcs1 = np.asarray(ts_pcs1)
    ts_pcs3 = np.asarray(ts_pcs3)


    ts_pcs1_highres = interp1d(x, ts_pcs1, kind='linear')(xnew)
    ts_pcs3_highres = interp1d(x, ts_pcs3, kind='linear')(xnew)

    ### PLOT ORIGINAL AND INTERPOLATED ICA TIMESERIES
    fig, ax = plt.subplots(1,1)
    ax.plot(x, ts_pcs1, 'o ', label='data PCS 1', color=pcs1.color)
    ax.plot(xnew, ts_pcs1_highres, ' x', label='interpolated', color=pcs1.color)

    ax.plot(x, ts_pcs3, 'o ', label='data PCS 2', color=pcs3.color)
    ax.plot(xnew, ts_pcs3_highres, ' x', label='interpolated', color=pcs3.color)
    ax.legend(loc='best')

    fig.savefig('{}PCS_timeseries_overlay.png'.format(outdir))

    ### RUN PYPLIS CROSSCORRELATION ANALYSIS
    (lag,
     coeffs,
     s1_ana,
     s2_ana,
     max_coeff_signal,
     (ax1, ax2, ax3)) = pyplis.plumespeed.find_signal_correlation(
                                         first_data_vec=ts_pcs1_highres,
                                         next_data_vec=ts_pcs3_highres,
                                         time_stamps=None, reg_grid_tres=None,
                                         freq_unit="S", itp_method="linear",
                                         max_shift_percent=50,
                                         sigma_smooth=0, plot=True)


    ax3.set_xlabel('Retrieved lag in units of timestep indices')
    ax3.figure.savefig('{}result_crosscorr_analysis.png'.format(outdir))

    ### COMPUTE VELOCITY
    dist_pcs_lines = pcs1.dist_other(pcs3)

    lag = float(lag)/res_fac
    print('dist_pcs_lines {}\n'
          'lag {}\n'
          'res_fac {}\n'
          .format( dist_pcs_lines, lag, res_fac ))


    velocity = dist_pcs_lines *DX / (lag*DT)
    print('Results:\n'
          'Lag between 2 integrated ICA time-series: {} [time indices]\n'
          'DT: {} s\n'
          'DX: {} m\n'
          'Retrieved velocity: {:.3f} m/s'
          .format(lag, DT, DX, velocity))

    print('Summary {} {} {} {} {}\n'
          .format(lag, COL_NUM1, DT, DX, velocity))


    print('VELOCITIES (with MULTIGAUSS) {:s} {:003d} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}'
          .format(loc,  COL_NUM2, np.nanmean(veff[1]), np.nanstd(veff[1]), velocity, np.nanmean(mus[1]), np.nanstd(mus[1])))
    print('txtfile', txtfile)
    ftxt.write('{:003d} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}\n'.format(COL_NUM2, np.nanmean(veff[1]), np.nanstd(veff[1]), velocity, np.nanmean(mus[1]), np.nanstd(mus[1])))
    ftxt.close()

    ftxt_noMG.write('{:003d} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}\n'.format(COL_NUM2, np.nanmean(veff[0]), np.nanstd(veff[0]), velocity, np.nanmean(mus[0]), np.nanstd(mus[0])))
    ftxt_noMG.close()
