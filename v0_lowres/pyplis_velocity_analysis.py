import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# -*- coding: utf-8 -*-
"""
Modified heavily from pyplis example script no. 3 - Plume background analysis

Arve Kylling, NILU.

This script .
"""
import numpy as np
import pyplis
#from matplotlib.pyplot import show, subplots, close
import matplotlib.pyplot as plt

import datetime
import Image

from glob import glob


FONTSIZE=20



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
                      print_MYSTIC_statistics=True, flip_y=False, Filter=None,
                      AddRows=False, verbose=True):
    """Load images

    """

    uvinp = 'uvspecCamW{:s}{:s}{:s}{:s}.inp'.format(wavelengthA,loc,timestartstr,experimenttype)
    uvout = 'uvspecCamW{:s}{:s}{:s}{:s}.out'.format(wavelengthA,loc,timestartstr,experimenttype)
    uvspecinputfile = filebase+uvinp
    uvspecoutputfile = filebase+uvout
    if verbose:
        print("uvspecinputfile", uvspecinputfile)
        print "uvspecoutputfile", uvspecoutputfile
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
        print "pix_heigth, pix_width", MYSTIC_AM.meta["pix_heigth"], MYSTIC_AM.meta["pix_width"]

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
        print "uvspecinputfileA0", uvspecinputfileA0
        print "uvspecoutputfileA0", uvspecoutputfileA0
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
    print "timestamps", timestamps
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
    
    CLEAR_OUTPUT = True
    plt.close("all")
    
    outdir = './output_new/'

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    if CLEAR_OUTPUT:
        
        for f in glob('{}*.png'.format(outdir)):
            os.remove(f)
    
    filebase      = './'
    filebaseBG    = './'
    loc           = 'LocA'
    wavelengthA   = '310'
    wavelengthB   = '330'
    timestart     = '001'
    timestep      = '1'
    timeend       = '005'
    
    timestampBG    = '005'
    experimenttype  = 'SUNSZA40WQBLC'
    experimenttypeBG  = 'SUNSZA40W_BG'
    experimenttypesensitivity  = ''
    experimenttypesensitivityBG  = ''
    Filter='Gaussian'#None#

    # Noise added to images (im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape))
    noiseamp = 0
    # time step between images
    DT=6.25 #s

    # Create an example plume intersection
    # Camera geometry
    DY=(3.828125+3.028125)/2.
    DZ=(3.903125+3.100000)/2.
    DY=0.250 # From Fig 2.
    DZ=0.125
    DS = np.sqrt(DZ*DZ+DY*DY)

    # Size of half image in the vertical, 5. = (96-86)/2.
    ZS=np.sin(np.deg2rad(5))*DS
    # 88/2 pixels in the vertical gives vertical pixel size
    dz = ZS / 44.0 #=0.0063961804401032104 #km =6.396180440m

    # From uvspec input file: mc_panorama_view 157.000000 203.000000 86.000000 96.000000
    # Size of half image in the horizontal 46/2
    XS=np.sin(np.deg2rad(23))*DS

    # 400/2 pixels in the horizontal gives horizontal pixel size
    DX = XS / 200. * 1000


    COL_NUM1 = 250
    COL_NUM2 = 290
    ROW_BOTTOM = 0+150
    ROW_TOP    = 110+150

    # Clear sky ROI
    BG_ROI = [0, 10, 0, 40]
    # Setup optical flow class (you may play with the parameters)
    
    ### DEFINE PCS LINES
    pcs1 = pyplis.LineOnImage(COL_NUM1, ROW_BOTTOM, COL_NUM1, ROW_TOP,
                             line_id="A", normal_orientation="right", #"left",#
                             color="#e67300")
    pcs2 = pyplis.LineOnImage(COL_NUM2, ROW_BOTTOM, COL_NUM2,  ROW_TOP,
                             line_id="B", normal_orientation="right", #"left",
                             color="#1a1aff")

    ### LOAD IMAGES
    ImagesAA = load_MYSTIC_images(filebase, wavelengthA,
                                  loc, timestart, timeend, timestep, experimenttype,
                                  wavelengthB, filebaseBG,
                                  timestampBG, experimenttypeBG, flip_y=True, 
                                  Filter=Filter,
                                  AddRows=True)
    
    ### LOAD SMALL IMAGES (NOT USED CURRENTLY)
    ImagesAA_small, (pcs1_small, pcs2_small) = downscale_images_and_pcs_lines(
            ImagesAA, (pcs1, pcs2))
    
    ### COMPUTE ICA TIMESERIES FOR BOTH LINES BASED ON 5 INPUT IMAGES
    # (includes plotting of AA images with lines in it)
    ts_pcs1 = []
    ts_pcs2 = []
    for i, im in enumerate(ImagesAA):
        im.img = im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape)
        # noisy = f + 0.4 * f.std() * np.random.random(f.shape)
        #            print im
        ax = im.show_img()
        ax.set_title("AA image {} (pyrlevel {})".format(i+1, im.pyrlevel))
        
        pcs1.plot_line_on_grid(ax=ax)
        pcs2.plot_line_on_grid(ax=ax)
        ax.legend()
        ax.figure.savefig('{}aa_img{}_N{}.png'.format(outdir, i,
                          noiseamp))
        
        ts_pcs1.append(sum(pcs1.get_line_profile(im.img)))
        ts_pcs2.append(sum(pcs2.get_line_profile(im.img)))
    
    
    ### INTERPOLATE THE 2 PCS timeseries to better register the lag 
    from scipy.interpolate import interp1d
    
    # original index
    x = range(len(ts_pcs1))
    
    res_fac = 100
    last_idx = x[-1]
    num_itp = res_fac * last_idx
    xnew = np.linspace(0, last_idx, num=num_itp, endpoint=True)
    
    # the two timeseries (each 5 timestamps)
    ts_pcs1 = np.asarray(ts_pcs1)
    ts_pcs2 = np.asarray(ts_pcs2)
    
    
    ts_pcs1_highres = interp1d(x, ts_pcs1, kind='linear')(xnew)
    ts_pcs2_highres = interp1d(x, ts_pcs2, kind='linear')(xnew)

    ### PLOT ORIGINAL AND INTERPOLATED ICA TIMESERIES
    fig, ax = plt.subplots(1,1)
    ax.plot(x, ts_pcs1, 'o ', label='data PCS 1', color=pcs1.color)
    ax.plot(xnew, ts_pcs1_highres, ' x', label='interpolated', color=pcs1.color)
    
    ax.plot(x, ts_pcs2, 'o ', label='data PCS 2', color=pcs2.color)
    ax.plot(xnew, ts_pcs2_highres, ' x', label='interpolated', color=pcs2.color)
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
                                         next_data_vec=ts_pcs2_highres, 
                                         time_stamps=None, reg_grid_tres=None, 
                                         freq_unit="S", itp_method="linear",
                                         max_shift_percent=50,
                                         sigma_smooth=0, plot=True)
    
    
    ax3.set_xlabel('Retrieved lag in units of timestep indices')
    ax3.figure.savefig('{}result_crosscorr_analysis.png'.format(outdir))
    
    ### COMPUTE VELOCITY
    dist_pcs_lines = pcs1.dist_other(pcs2)
    
    lag = float(lag)/res_fac
    velocity = dist_pcs_lines *DX / (lag*DT)
    print('Results:\n'
          'Lag between 2 integrated ICA time-series: {} [time indices]\n'
          'DT: {} s\n'
          'DX: {} m\n'
          'Retrieved velocity: {:.3f} m/s'
          .format(lag, DT, DX, velocity))
    