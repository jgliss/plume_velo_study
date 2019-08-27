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
from os.path import join, basename
import pyplis
#from matplotlib.pyplot import show, subplots, close
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
from matplotlib.legend import rcParams
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.pyplot import close, sca, figure, yticks, xticks, subplots
from scipy import stats
import argparse
import datetime
import Image
from SETTINGS import IMG_DIR
from scipy import ndimage
import sys
import Experiment as Exp
from glob import glob
import copy

FONTSIZE=20

def LES2Img(path, fns, axis=1, verbose=False):
    zoom=True
    Plume=Exp.Plume()
    for fn in fns:
        # Read as x,y,z
        Plume.ReadLESNetCDF(fn,verbose=False, scale_axis_factor=312.5)

        if zoom:
            indx = np.where((Plume.x > 0.055) & (Plume.x < 0.275))
            Plume.dens = Plume.dens[indx[0][0]:indx[0][len(indx[0])-1],:,:]
            Plume.x  = Plume.x[indx[0][0]:indx[0][len(indx[0])-1]]
            indz = np.where((Plume.z > 0.0) & (Plume.z < 0.05))
            Plume.dens = Plume.dens[:,:,indz[0][0]:indz[0][len(indz[0])-1]]
            Plume.z  = Plume.z[indz[0][0]:indz[0][len(indz[0])-1]]

        img=np.trapz(Plume.dens,axis=axis)

        IMG = Img()
        IMG.img = img.T
        IMG.x  = Plume.x
        IMG.z  = Plume.z
        IMG.meta["pix_width"]  = IMG.img.shape[1]
        IMG.meta["pix_heigth"] = IMG.img.shape[0]
        IMG.statistics()
        IMG.threshold = 0.00

    return IMG

def fractal_dimension(Z, threshold=0.9, ax=None):
    # Only for 2d image
    assert(len(Z.shape) == 2)

    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)

        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])

    # Transform Z into a binary array
    Z = (Z < threshold)

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    print "fractal_dimension p,n", p, n, Z.shape

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)
    # FIXME!!!!!!!!!!!!!!!!!!!!!!
    sizes = np.array([64, 32, 24, 16, 14, 12, 10, 8, 6, 4, 2])
    sizes = np.array([64, 32, 24, 16, 14, 12, 10, 8, 6, 4, 2])
    sizes = np.array([64, 32, 24, 14, 12,  6, ])
    sizes = np.array([16,  10, 8,  4])
    print "fractal_dimension p,n", n, sizes

    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        count=boxcount(Z, size)
        counts.append(count)
        print size, count

    #    print("sizes", sizes, np.log(sizes), np.log(counts))
    #    print("counts", counts)

    # Fit the successive log(sizes) with log (counts)
    #    ret = np.polyfit(np.log(sizes), np.log(counts), 1, full=True)
    #    print ret
    #    exit()
    fit_sizes = sizes#[:3]
    fit_counts = counts#[:3]

    coeffs, residuals, rank, singular_values, rcond = np.polyfit(np.log(fit_sizes), np.log(fit_counts), 1, full=True)
#    print sizes, counts
#    print np.log(sizes), np.log(counts)
    print "coeffs, stats", coeffs, residuals, threshold

    plot=False#True #
    if ax !=None:
        plot=True
    if plot:
        p=ax.plot(np.log(sizes), np.log(counts), linestyle='-', marker='o')
        pf=ax.plot(np.log(fit_sizes), np.log(fit_counts), linestyle='-', marker='x', color='r')
        a=coeffs[0]
        b=coeffs[1]
        c=p[0].get_color()
        ax.plot(np.log(sizes), a*np.log(sizes)+b, linestyle=':', color=c)

    return -coeffs[0]

def plot_histo_result(flow, lines):
    """Plot result of histogram analysis

    Parameters
    ----------
    flow : OptflowFarneback
        optical flow calculation class containing pre-computed flow field
    lines : list
        list containing LineOnImage objects

    Returns
    -------
    figure
        matplotlib figure instance
    """
    if isinstance(lines, pyplis.LineOnImage):
        lines = [lines]
    fig = figure(figsize=(14,8))
    ax0 = fig.add_axes([0.01, 0.15, 0.59, 0.8])
    #ax0.set_axis_flowf()
    ax1 = fig.add_axes([0.61, 0.15, 0.16, 0.8])
    ax2 = fig.add_axes([0.78, 0.15, 0.16, 0.8])
    kwargs = {"linewidth": 1}
    flow.plot(ax=ax0, **kwargs)#, in_roi=True)
    for line in lines:
        m = line.get_rotated_roi_mask(flow.flow.shape[:2])
        line.plot_line_on_grid(ax=ax0, include_normal=1,
                               include_roi_rot=1)

        try:
            _, mu, sigma = flow.plot_orientation_histo(pix_mask=m,
                                                       apply_fit=True,
                                                       ax=ax1,
                                                       color=line.color)
            ax1.legend_.remove()
            low, high = mu-sigma, mu+sigma
            flow.plot_length_histo(pix_mask=m, ax=ax2, dir_low=low, label=line.line_id,
                                 dir_high=high, color=line.color)
        except Exception as e:
            warnings.warn(repr(e))
    #pyplis.helpers.set_ax_lim_roi(roi_disp, ax0)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    ax0.set_title("")
    #ax1.set_title("Directions", fontsize=12)
    ax1.set_title("")
    ax1.set_xlabel(r"$\varphi\,[^\circ]$", fontsize=20)
    #ax1.set_ylim([0, 1600]) #hard coded
    ax1.get_yaxis().set_ticks([])
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_xlabel(r"$|\mathbf{f}|$ [pix]", fontsize=20)
    ax2.set_ylabel("Counts / bin", fontsize=20)

    #ax1.yaxis.set_visible(False)

    #fig.tight_layout()
    #ax2.set_title("Magnitude (pyr=%d)" %fl.pyrlevel, fontsize=12)
    ax2.set_title("")
    #ax2.legend_.remove()
    #subplots_adjust()

    sca(ax1)
    xticks(rotation=40, ha="right")
    sca(ax2)
    yticks(rotation=90, va="center")

    #fig.tight_layout()
    return fig

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
            print "statistics", xp, yp, self.img.shape

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
        print "uvspecinputfile", uvspecinputfile
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
    T0 = datetime.datetime(2017,1,1,0,0,0) # define arbitrary start time
    times = [T0 + datetime.timedelta(x * DT / 86400.0) for x in range(len(timestamps))]

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



def plot_img(img, **kwargs):
    if 'figsize' in kwargs:
        figsize=kwargs['figsize']
    else:
        figsize=(18,7.2)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    rcParams.update({'font.size':FONTSIZE})
#    fig.subplots_adjust(top=0.99, right=0.98,left=0.09,bottom=0.02)
    fig.subplots_adjust(top=0.99, right=0.98,left=0.09,bottom=0.02)

    if 'extent' in kwargs:
        extent=kwargs['extent']
    else:
        extent=None
    if 'flip_x' in kwargs:
        if kwargs['flip_x']:
            img=img[:,::-1]
    if kwargs['log_img']:
        data = np.log10(img)
        data = np.where(data>np.log10(kwargs['threshold']), data, np.log10(kwargs['threshold']))
    else:
        data = img
        data = np.where(data>kwargs['threshold'], data, kwargs['threshold'])

    print("extent", extent)
    if 'vmin' in kwargs:
        im = ax.imshow(data, origin='lower', cmap=kwargs['cmap'], extent=extent,
                       vmin=kwargs['vmin'], vmax=kwargs['vmax'])
    else:
        im = ax.imshow(data, origin='lower', cmap=kwargs['cmap'], extent=extent)
    ax.set_title(kwargs['title'], fontsize=FONTSIZE)
    ax.set_xlabel(kwargs['xlabel'], fontsize = FONTSIZE)
    ax.set_ylabel(kwargs['ylabel'], fontsize = FONTSIZE)
    if 'cb_shrink' in kwargs:
        shrink=kwargs['cb_shrink']
    else:
        shrink=0.64

    if 'cb_ticks' in kwargs:
        plt.colorbar(im, shrink=shrink, pad=0.02, ticks=kwargs['cb_ticks'])
    else:
        plt.colorbar(im, shrink=shrink, pad=0.02)

def plot_images_column(**kwargsall):

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    nplots=len(kwargsall)

    # nplots subplots, the axes array is 1-d
    fig, axarr = plt.subplots(nplots, sharex=True, figsize=(14,15))

    ip=0
    for items in sorted(kwargsall.items()):
        kwargs = items[1]
        img=kwargs['data']
        if kwargs['log_img']:
            data = np.log10(img)
            data = np.where(data>np.log10(kwargs['threshold']), data, np.log10(kwargs['threshold']))
        else:
            data = img
            data = np.where(data>kwargs['threshold'], data, kwargs['threshold'])

        yp = np.arange(0, kwargs['data'].shape[1])

        if 'vmin' in kwargs:
            im=axarr[ip].imshow(data, origin='lower', aspect='auto', extent=kwargs['extent'],
                                cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'])
        else:
            im=axarr[ip].imshow(data, origin='lower', aspect='auto', extent=kwargs['extent'], cmap=kwargs['cmap'])
        if ip==nplots-1:
            axarr[ip].set_xlabel(kwargs['xlabel'], fontsize=FONTSIZE)
        axarr[ip].set_ylabel(kwargs['ylabel'], fontsize=FONTSIZE)
        axarr[ip].set_title(kwargs['title'], fontsize=FONTSIZE)

        # Colorbars are a bit of a mess, from:
        # https://matplotlib.org/gallery/axes_grid1/demo_axes_divider.html#sphx-glr-gallery-axes-grid1-demo-axes-divider-py
        divider = make_axes_locatable(axarr[ip])
        ax_cb = divider.new_horizontal(size="2%", pad=0.05)
        fig1 = axarr[ip].get_figure()
        fig1.add_axes(ax_cb)
        plt.colorbar(im, cax=ax_cb)
        ax_cb.yaxis.tick_right()
        ax_cb.yaxis.set_tick_params(labelright=True)

        ip=ip+1

    rcParams.update({'font.size':FONTSIZE})
    fig.subplots_adjust(top=0.96, right=0.93,left=0.08,bottom=0.05)


def plot_img_statistics(Img, **kwargs):
    fig = plt.figure(1, figsize=(14,15))
    rcParams.update({'font.size':FONTSIZE})
    fig.subplots_adjust(top=0.96, right=0.93,left=0.08,bottom=0.1)
    fig.subplots_adjust(hspace=0.4)
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    yp = np.arange(0, Img.shape[1])

    includeLES=False
    if kwargs['LESdata'] != None:
        includeLES=True
        LESImg = kwargs['LESdata']

    colorline1='b'
    colorline2='r'
    colorline3='g'
    colorline4='k'
    if kwargs['log_img']:
        data = np.log10(Img.img)
        data = np.where(data>np.log10(Img.threshold), data, np.log10(Img.threshold))
    else:
        data = Img.img
        data = np.where(data>Img.threshold, data, Img.threshold)

    print("data min/max", data.min(), data.max())
    ax1.imshow(data, aspect='auto', origin='lower', extent=kwargs['extent'], cmap=cmaps.Greys)
    ax1.set_ylabel('Vertical pixel number', fontsize=FONTSIZE)
    ax1.set_xlabel('Horizontal pixel number', fontsize=FONTSIZE)
    if kwargs['title']!='':
        ax1.set_title(kwargs['title'], fontsize=FONTSIZE)

    box = ax1._position.bounds
    height = box[3]
    tmpax = fig.add_axes([box[0], box[1] , box[2], height])
    tmpax.set_axis_off()
    linewidth=3
    Img.phis = np.arange(Img.meta["phi1"], Img.meta["phi2"], (Img.meta["phi2"]-Img.meta["phi1"])/Img.img.shape[1])-180
    PlotAng=True #False #
    indx = np.where(Img.Centerline>0.0)
    Img.Centerline =  Img.Centerline[indx[0]]
    if PlotAng:
        tmpy = Img.phis[indx[0]]
        AngleShift=1.0
        tmpax.plot(tmpy, Img.Centerline,c=colorline2, lw=linewidth, linestyle=':')
        tmpax.set_xlim(Img.phis[0],Img.phis[len(Img.phis)-1])
    else:
        tmpyp = yp[indx[0]]
        tmpax.plot(tmpyp, Img.Centerline,c=colorline2, lw=linewidth, linestyle=':')
        tmpax.set_xlim(0,len(yp)-1)
    tmpax.set_ylim(Img.shape[0],0) #Upside down since origin=lower
    #    tmpax.set_xlim(Img.meta["phi1"]-180, Img.meta["phi2"]-180)
    if includeLES:
        D = 0.26  # Distance from camera to plume in km. Roughly estimated from Fig. 2 and tuned
        xcentre =  LESImg.x[int(LESImg.x.shape[0]/2)]
        phis = np.rad2deg(np.arctan((LESImg.x-xcentre)/D) )
        #        print "LESImg.x", LESImg.x-xcentre
        #        print "phis", phis[0], phis[len(phis)-1]
        #        print "Img.thetas", Img.thetas.shape, Img.thetas
        tmpCL = np.interp(Img.phis, phis, LESImg.Centerline)
        #        print "tmpCL", tmpCL.shape, tmpCL.min(), tmpCL.max(), Img.Centerline.min(), Img.Centerline.max()
        xpixoffset=10
        xoffset= 12 #9#*0.115  # 46/400=0.115
        yoffset=9 #+88-78
        yscale = Img.shape[0]/float(LESImg.shape[0])
        #        tmpax.plot(Img.phis+xoffset, tmpCL*yscale+yoffset,c=colorline2, linestyle=':', lw=linewidth)
        indxdiff1=20 #0
        indxdiff2=len(tmpCL)-20
        tmpCL=tmpCL*yscale+yoffset
        if PlotAng:
            tmpax.plot(Img.phis+AngleShift, tmpCL,c=colorline2, linestyle='-', lw=linewidth)
        else:
            tmpax.plot(yp[indxdiff1-xoffset:indxdiff2-xoffset]+xpixoffset, tmpCL[indxdiff1-xoffset:indxdiff2-xoffset],c=colorline2, linestyle='-', lw=linewidth)
        #        CL_diff = (Img.Centerline[indxdiff1:indxdiff2]-tmpCL[indxdiff1-xoffset:indxdiff2-xoffset])/tmpCL[indxdiff1:indxdiff2]
        # CL_diff = (Img.Centerline[indxdiff1:indxdiff2]-tmpCL[indxdiff1-xoffset:indxdiff2-xoffset])
        # fn = 'Stats_CL_V01.txt'
        # fp=open(fn,'w')
        # for  CLimg, CLLES, diff in zip(Img.Centerline[indxdiff1:indxdiff2], tmpCL[indxdiff1-xoffset:indxdiff2-xoffset], CL_diff):
        #    fp.write('{:f} {:f} {:f}\n'.format(CLimg, CLLES, diff))
        # close(fn)

        ximgind1 = 20
        ximgind2 = len(Img.Centerline)#-ximgind1
        ximg = np.arange(ximgind1, ximgind2)
        xlesoffset=0

        xlesind1 = 11
        xlesind2 = len(LESImg.Centerline)-10#-xlesind1

        fn = 'Stats_CL_V03.txt'
        fp=open(fn,'a')
        tmpimg = Img.Centerline[ximgind1:ximgind2]
        tmples = LESImg.Centerline[xlesind1:xlesind2]
        xles =  np.linspace(ximgind1, ximgind2, len(tmples))
        print xles.shape, ximg.shape, tmples.shape, LESImg.Centerline.shape, xles.min(), xles.max(), ximg.min(), ximg.max()
        tmpCL = yoffset+yscale*np.interp(ximg, xles, tmples)
        CL_diff =tmpimg - tmpCL
        for  CLimg, CLLES, diff in zip(tmpimg, tmpCL, CL_diff):
           fp.write('{:f} {:f} {:f}\n'.format(CLimg, CLLES, diff))
        close(fn)

    #    print("Dispersion", np.nanmin(Img.Dispersion), np.nanmax(Img.Dispersion))
    #    print "Img.Dispersion", Img.Dispersion.shape, Img.Dispersion
    indx = np.where(Img.Dispersion>0.0)
    Img.Dispersion =  Img.Dispersion[indx[0]]
    Img.RelativeDispersion =  Img.RelativeDispersion[indx[0]]

    pls = []
    if PlotAng:
        xx = np.arange(0,len(Img.phis))
        tmpxx = np.arange(0,len(Img.phis[indx[0]]))
        AngleOffSet=180
        #        pl, = ax2.plot(Img.phis+AngleOffSet,Img.Dispersion,c=colorline1)
        pl, = ax2.plot(tmpxx, Img.Dispersion,c=colorline1, linestyle=':')
        #pls.append(pl)
        #        pl, = ax2.plot(Img.phis+AngleOffSet,Img.RelativeDispersion,c=colorline2)
        pl, = ax2.plot(tmpxx, Img.RelativeDispersion,c=colorline2, linestyle=':')
        #        pls.append(pl)
        #        print xx[len(xx)-1], xx
        #        ax2.set_xlim(0,xx[len(xx)-1])
        ax2.set_xlim(8,400)
    else:
        tmpyp = yp[indx[0]]
        pl, = ax2.plot(tmpyp,Img.Dispersion,c=colorline1, linestyle=':')
        #pls.append(pl)
        pl, = ax2.plot(tmpyp,Img.RelativeDispersion,c=colorline2, linestyle=':')
        #        pls.append(pl)
        ax2.set_xlim(0,len(yp)-1)
        #        ax2.set_xlim(1,400)

    #    ax2.set_ylim(0,20)
    ax2.set_ylim(0.1,20)
    # tmpax.set_xlim(Img.meta["phi1"]-180, Img.meta["phi2"]-180)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('Dispersion', fontsize=FONTSIZE)
    ax2.set_xlabel('Log$_{10}$ of horizontal pixel number', fontsize=FONTSIZE)
#    ax2.legend((pls[0], pls[1]), ('Absolute dispersion','Relative dispersion',),
#                        prop={'size':FONTSIZE}, loc='lower right')

    if includeLES:
        tmpDIS = np.interp(Img.phis, phis, LESImg.Dispersion)
        # tmpDIS = np.interp(yp, phis, LESImg.Dispersion)
        #        yscale = Img.shape[0]/float(LESImg.shape[0])
        yoffset=0
        #
        tmpDIS=tmpDIS*yscale+yoffset
        tmpDIS2 = np.interp(Img.phis, phis, LESImg.RelativeDispersion)
        tmpDIS2 =tmpDIS2*yscale+yoffset
        if PlotAng:
            PixelShift=1
            pl, = ax2.plot(xx+PixelShift, tmpDIS,c=colorline1, linestyle='-', lw=linewidth)
            #            ax2.plot(Img.phis+AngleShift+AngleOffSet, tmpDIS,c=colorline1, linestyle=':', lw=linewidth)
            pls.append(pl)
            pl, = ax2.plot(xx+PixelShift, tmpDIS2,c=colorline2, linestyle='-', lw=linewidth)
            pls.append(pl)
            #            ax2.plot(Img.phis+AngleShift+AngleOffSet, tmpDIS2*yscale+yoffset,c=colorline2, linestyle=':', lw=linewidth)
        else:
            pl, = ax2.plot(yp+xoffset, tmpDIS,c=colorline1, linestyle='-', lw=linewidth)
            pls.append(pl)
            pl, = ax2.plot(yp+xoffset, tmpDIS2, c=colorline2, linestyle='-', lw=linewidth)
            pls.append(pl)
        indx1=10
        indx2=30 # 60#
        indx = np.where((yp>indx1) & (yp<indx2))
        x = yp[indx]
        y = Img.RelativeDispersion[indx]
        m,b = np.polyfit(x, y, 1)
        print "m,b img", m,b
        indx = np.where((yp>indx1+xoffset) & (yp<indx2+xoffset))
        x = yp[indx]
        y = LESImg.RelativeDispersion[indx]
        m,b = np.polyfit(x, y, 1)
        print "m,b LES", m,b, pls

        ax2.legend((pls[0], pls[1]), ('Absolute dispersion','Relative dispersion',),
                   prop={'size':FONTSIZE}, loc='lower right')
        # AD_diff = (Img.Dispersion[indxdiff1:indxdiff2]-tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset])
        # fn = 'Stats_AD_V01.txt'
        # fp=open(fn,'a')
        # for  ADimg, ADLES, diff in zip(Img.Dispersion[indxdiff1:indxdiff2], tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset], AD_diff):
        #    fp.write('{:f} {:f} {:f}\n'.format(ADimg, ADLES, diff))
        # close(fn)

        # RD_diff = (Img.RelativeDispersion[indxdiff1:indxdiff2]-tmpDIS2[indxdiff1-xoffset:indxdiff2-xoffset])
        # fn = 'Stats_RD_V01.txt'
        # fp=open(fn,'a')
        # for  RDimg, RDLES, diff in zip(Img.RelativeDispersion[indxdiff1:indxdiff2], tmpDIS2[indxdiff1-xoffset:indxdiff2-xoffset], RD_diff):
        #    fp.write('{:f} {:f} {:f}\n'.format(RDimg, RDLES, diff))
        # close(fn)

        fn = 'Stats_AD_V03.txt'
        fp=open(fn,'a')
        tmpimg = Img.Dispersion[ximgind1:ximgind2]
        tmples = LESImg.Dispersion[xlesind1:xlesind2]
        tmpCL =  yscale*np.interp(ximg, xles, tmples)
        CL_diff =tmpimg - tmpCL
        for  CLimg, CLLES, diff in zip(tmpimg, tmpCL, CL_diff):
           fp.write('{:f} {:f} {:f}\n'.format(CLimg, CLLES, diff))
        close(fn)

        fn = 'Stats_RD_V03.txt'
        fp=open(fn,'a')
        tmpimg = Img.RelativeDispersion[ximgind1:ximgind2]
        tmples = LESImg.RelativeDispersion[xlesind1:xlesind2]
        tmpCL =  yscale*np.interp(ximg, xles, tmples)
        CL_diff =tmpimg - tmpCL
        for  CLimg, CLLES, diff in zip(tmpimg, tmpCL, CL_diff):
           fp.write('{:f} {:f} {:f}\n'.format(CLimg, CLLES, diff))
        close(fn)

    print("AbsoluteSkewness", np.nanmin(Img.AbsoluteSkewness), np.nanmax(Img.AbsoluteSkewness))
    print "Img.AbsoluteSkewness", Img.Skewness.shape, Img.Skewness
    indx = np.where(Img.Skewness!=0.0)
    Img.Skewness =  Img.Skewness[indx[0]]
    print "Img.AbsoluteSkewness", Img.Skewness.shape, Img.Skewness

    if PlotAng:
        #        ax3.plot(Img.phis, Img.AbsoluteSkewness,c=colorline1, lw=linewidth, linestyle=':')
        ax3.plot(Img.phis[indx[0]], Img.Skewness,c=colorline1, lw=linewidth, linestyle=':')
        ax3.set_xlim(Img.phis[0],Img.phis[len(Img.phis)-1])
    else:
        #        ax3.plot(yp,Img.AbsoluteSkewness,c=colorline1, linestyle=':')
        ax3.plot(yp[indx[0]],Img.Skewness,c=colorline1, linestyle=':')
        ax3.set_xlim(0,len(yp))

    ax3.set_ylim(-2,2)
#    ax3.set_xlabel('Horizontal pixel number', fontsize=FONTSIZE)
    ax3.set_xlabel('Horizontal viewing angle', fontsize=FONTSIZE)
    ax3.set_ylabel('Skewness', fontsize=FONTSIZE)

    if includeLES:
        tmpDIS = np.interp(Img.phis, phis, LESImg.Skewness)
        #        tmpDIS = np.interp(Img.phis, phis, LESImg.AbsoluteSkewness)
        # tmpDIS = np.interp(yp, phis, LESImg.Dispersion)
        yscale = 1 #Img.shape[0]/float(LESImg.shape[0])
        yoffset=0
        #        ax2.plot(Img.phis+xoffset, tmpDIS*yscale+yoffset,c=colorline1, linestyle=':', lw=linewidth)
        if PlotAng:
            ax3.plot(Img.phis+AngleShift, tmpDIS*yscale+yoffset,c=colorline1, linestyle='-', lw=linewidth)
        else:
            ax3.plot(yp[indxdiff1-xoffset:indxdiff2-xoffset]+xpixoffset, tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset]*yscale+yoffset,c=colorline1, linestyle='-', lw=linewidth)

        #        SK_diff = (Img.AbsoluteSkewness[indxdiff1:indxdiff2]-tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset])/tmpDIS[indxdiff1:indxdiff2]
        # SK_diff = (Img.AbsoluteSkewness[indxdiff1:indxdiff2]-tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset])
        # fn = 'Stats_Skewness_V01.txt'
        # fp=open(fn,'a')
        # for  SKimg, SKLES, diff in zip(Img.AbsoluteSkewness[indxdiff1:indxdiff2], tmpDIS[indxdiff1-xoffset:indxdiff2-xoffset], SK_diff):
        #    fp.write('{:f} {:f} {:f}\n'.format(SKimg, SKLES, diff))
        # close(fn)

        fn = 'Stats_Skewness_V03.txt'
        fp=open(fn,'a')
        #        tmpimg = Img.AbsoluteSkewness[ximgind1:ximgind2]
        tmpimg = Img.Skewness[ximgind1:ximgind2]
        #        tmples = LESImg.AbsoluteSkewness[xlesind1:xlesind2]
        tmples = LESImg.Skewness[xlesind1:xlesind2]
        tmpCL = yscale*np.interp(ximg, xles, tmples)
        CL_diff =tmpimg - tmpCL
        for  CLimg, CLLES, diff in zip(tmpimg, tmpCL, CL_diff):
           fp.write('{:f} {:f} {:f}\n'.format(CLimg, CLLES, diff))
        close(fn)


def plot_img_statistics_5(Img, **kwargs):
    fig = plt.figure(1, figsize=(14,18))
    rcParams.update({'font.size':FONTSIZE})
    fig.subplots_adjust(top=0.96, right=0.93,left=0.08,bottom=0.05)
    ax1 = plt.subplot(511)
    ax2 = plt.subplot(512)
    ax3 = plt.subplot(513)
    ax4 = plt.subplot(514)
    ax5 = plt.subplot(515)
    yp = np.arange(0, Img.shape[1])

    colorline1='b'
    colorline2='r'
    if kwargs['log_img']:
        data = np.log10(Img.img)
        data = np.where(data>np.log10(Img.threshold), data, np.log10(Img.threshold))
    else:
        data = Img.img
        data = np.where(data>Img.threshold, data, Img.threshold)

    print("data min/max", data.min(), data.max())
    #    print ("Centerline", Img.Centerline)
    ax1.imshow(data, aspect='auto', origin='lower', extent=kwargs['extent'], cmap=kwargs['cmap'])#, cmap=cmaps.gist_ncar_r)#
    ax1.set_ylabel('Vertical pixel number', fontsize=FONTSIZE)
    ax1.set_title(kwargs['title'], fontsize=FONTSIZE)

    box = ax1._position.bounds
    height = box[3]
    tmpax = fig.add_axes([box[0], box[1] , box[2], height])
    tmpax.set_axis_off()
    tmpax.plot(yp,Img.Centerline,c=colorline2)
    #    tmpax.plot(np.sin(np.linspace(0,np.random.randint(20,1000),1000))*0.4)
    tmpax.set_ylim(Img.shape[0],0) #Upside down since origin=lower
    tmpax.set_xlim(0,len(yp)-1)

    #    yval=Img.shape[0]-Img.Centerline # Turn upside down
    yval=Img.Centerline
    #    ax2.plot(yp,Img.Centerline,c=colorline1)
    ax2.plot(yp, yval, c=colorline1)
    #    ax2.set_ylim(ax2.get_ylim()[::-1])
    #    ax2.set_ylim(20,40)
    ax2.set_ylabel('Centerline', fontsize=FONTSIZE)
    data = Img.Centerline
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    print('{:20s} {:7.3f} {:7.2f} {:7d} {:7d} {:7.3f}'.format('Centerline', np.nanmin(Img.Centerline), np.nanmax(Img.Centerline), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0])))
    #    ax2.set_ylim(Img.shape[0],0) #Upside down since origin=lower
#    ax2.set_ylim(0,Img.shape[0])
    ax2.set_ylim(-2,2)
    ax2.plot([0,len(yp)-1],[0,0],c='k')
    ax2.plot([0,len(yp)-1],[1,1],c='k',linestyle=':')
    ax2.plot([0,len(yp)-1],[-1,-1],c='k',linestyle=':')
    ax2.set_xlim(0,len(yp)-1)

    data = Img.Dispersion
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    print('{:20s} {:7.3f} {:7.2f} {:7d} {:7d} {:7.3f}'.format('Dispersion',np.nanmin(Img.Dispersion), np.nanmax(Img.Dispersion), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0])))
    ax3.plot(yp,Img.Dispersion,c=colorline1)
    #    ax3.set_ylim(0,20)
    ax3.set_ylim(-2,2)
    ax3.plot([0,len(yp)-1],[0,0],c='k')
    ax3.plot([0,len(yp)-1],[1,1],c='k',linestyle=':')
    ax3.plot([0,len(yp)-1],[-1,-1],c='k',linestyle=':')
    ax3.set_ylabel('Abs. dispersion', fontsize=FONTSIZE)
    ax3.set_xlim(0,len(yp)-1)

    data = Img.RelativeDispersion
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    print('{:20s} {:7.3f} {:7.2f} {:7d} {:7d} {:7.3f}'.format('RelativeDispersion', np.nanmin(Img.RelativeDispersion), np.nanmax(Img.RelativeDispersion), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0])))
    ax4.plot(yp,Img.RelativeDispersion,c=colorline1)
    #    ax4.set_ylim(0,15)
    ax4.set_ylim(-2,2)
    ax4.plot([0,len(yp)-1],[0,0],c='k')
    ax4.plot([0,len(yp)-1],[1,1],c='k',linestyle=':')
    ax4.plot([0,len(yp)-1],[-1,-1],c='k',linestyle=':')
    ax4.set_ylabel('Rel. dispersion', fontsize=FONTSIZE)
    ax4.set_xlim(0,len(yp)-1)

    data = Img.AbsoluteSkewness
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    print('{:20s} {:7.3f} {:7.2f} {:7d} {:7d} {:7.3f}'.format('AbsoluteSkewness', np.nanmin(Img.AbsoluteSkewness), np.nanmax(Img.AbsoluteSkewness), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0])))
    ax5.plot(yp,Img.AbsoluteSkewness,c=colorline1)
    ax5.set_ylim(-4,4)
    ax5.plot([0,len(yp)-1],[0,0],c='k')
    ax5.plot([0,len(yp)-1],[1,1],c='k',linestyle=':')
    ax5.plot([0,len(yp)-1],[-1,-1],c='k',linestyle=':')
    ax5.set_xlabel('Horizontal pixel number', fontsize=FONTSIZE)
    ax5.set_ylabel('Abs. skewness', fontsize=FONTSIZE)
    ax5.set_xlim(0,len(yp)-1)


def print_img_statistics(Img, **kwargs):

    mintest=-0.0
    maxtest= 0.0

    print('{:20s} {:7s} {:7s} {:7s} {:7s} {:7s} {:7s} {:7s} {:7s} {:7s} {:7s}'.format('', 'Img min', 'Img max', 'Img avg','Img std', 'Img med', '  N>1', 'Img sha', 'N>1 (%)', 'Img<min', 'Img>max'))


    if 'tau_A' in kwargs.keys():
        tau = kwargs['tau_A']
        data = tau.img
        data = data[~np.isnan(data)]
        indx= np.where(np.absolute(data)>1.0)
        indn = np.where(data<mintest)
        indp = np.where(data>maxtest)
        print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('tau', np.nanmin(tau.img), np.nanmax(tau.img), np.nanmean(tau.img), np.nanstd(tau.img), np.nanmedian(tau.img), len(indx[0]),\
                                                                              data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

        indy = np.where(np.abs(data) > 0.03)
        data = data[indy[0]]
        indx= np.where(np.absolute(data)>1.0)
        indn = np.where(data<mintest)
        indp = np.where(data>maxtest)
        print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('tau plume only', np.nanmin(data), np.nanmax(data), np.nanmean(data), np.nanstd(data), np.nanmedian(data), len(indx[0]),\
                                                                              data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))


    data = Img.img
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('Data', np.nanmin(Img.img), np.nanmax(Img.img), np.nanmean(Img.img), np.nanstd(Img.img), np.nanmedian(Img.img), len(indx[0]),\
                                                                              data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

    data = Img.img
    indy = np.where(np.abs(data) > 0.0)
    #    print "Img",  Img.img.shape, Img.img.min(), Img.img.max(), Img.img
    #   print data.shape, data.min(), data.max()
    data = data[indy[0]]
    #   print data.shape, data.min(), data.max()
    #data = data[~np.isnan(data)]
    #    print data.shape, data.min(), data.max()
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('Data plume only', np.nanmin(data), np.nanmax(data), np.nanmean(data), np.nanstd(data), np.nanmedian(data), len(indx[0]),\
                                                                              data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))


    data = Img.Centerline
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('Centerline', np.nanmin(Img.Centerline), np.nanmax(Img.Centerline), np.nanmean(Img.Centerline), np.nanstd(Img.Centerline), np.nanmedian(Img.Centerline), len(indx[0]),\
                                                                              data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

    data = Img.Dispersion
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('Dispersion',np.nanmin(Img.Dispersion), np.nanmax(Img.Dispersion), np.nanmean(Img.Dispersion), np.nanstd(Img.Dispersion), np.nanmedian(Img.Dispersion), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

    data = Img.RelativeDispersion
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('RelativeDispersion', np.nanmin(Img.RelativeDispersion), np.nanmax(Img.RelativeDispersion), np.nanmean(Img.RelativeDispersion), np.nanstd(Img.RelativeDispersion), np.nanmedian(Img.RelativeDispersion), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

    data = Img.AbsoluteSkewness
    data = data[~np.isnan(data)]
    indx= np.where(np.absolute(data)>1.0)
    indn = np.where(data<mintest)
    indp = np.where(data>maxtest)
    print('{:20s} {:7.3f} {:7.2f} {:7.4f} {:7.4f} {:7.4f} {:7d} {:7d} {:7.4f} {:7d} {:7d}'.format('AbsoluteSkewness', np.nanmin(Img.AbsoluteSkewness), np.nanmax(Img.AbsoluteSkewness) , np.nanmean(Img.AbsoluteSkewness), np.nanstd(Img.AbsoluteSkewness), np.nanmedian(Img.AbsoluteSkewness), len(indx[0]),\
        data.shape[0], 100*len(indx[0])/float(data.shape[0]), len(indn[0]), len(indp[0])))

    return

### SCRIPT MAIN FUNCTION
if __name__=="__main__":

    plt.close("all")
    
    OUT_DIR = './output/'
    import os, glob

    parser = argparse.ArgumentParser(description='Plot camera images and apparant absorbance. '+
                                 'Script requires that both foreground and background simulations'+
                                 'are available in folder specified by InputFolder.')


    parser.add_argument('InputFolder', type=str,
                        help='Full folder name holding experiment')
    parser.add_argument('--PlotType', type=str, default='Images',
                        help='Type of plot/analysis to be made')
    parser.add_argument('--InputFolderBG', type=str, default='',
                        help='Full folder name holding experiment')
    parser.add_argument('Camera', type=str,
                        help='Camera location to be plotted')
    parser.add_argument("--timestart", help="Start time stamp of LES simulation", default='')
    parser.add_argument("--timeend", help="End time stamp of LES simulation", default='')
    parser.add_argument("--timestep", help="Time step in stamp units of LES simulation", default='1')
    parser.add_argument("--timestampBG", help="Time stamp of BG simulation", default='')
    parser.add_argument("--wavelengthA", help="wavelength name of camera A", default='310')
    parser.add_argument("--wavelengthB", help="wavelength name of camera B", default='330')
    parser.add_argument('--experimenttype', type=str, default='',
                        help='String defining type of experiment. See code for more info.')
    parser.add_argument('--experimenttypeBG', type=str, default='-999',
                        help='String defining type of BG experiment. See code for more info.')
    parser.add_argument('--experimenttypesensitivity', type=str, default='',
                        help='String defining type of sensitivity experiment. See code for more info.')
    parser.add_argument('--experimenttypesensitivityBG', type=str, default='',
                        help='String defining type of BG sensitivity experiment. See code for more info.')
    parser.add_argument("--pngfolder", help="Name of png file folder", default='')

    args = parser.parse_args()
    filebase      = args.InputFolder
    filebaseBG    = args.InputFolderBG
    PlotType      = args.PlotType
    loc           = args.Camera
    pngfolder     = args.pngfolder
    wavelengthA   = args.wavelengthA
    wavelengthB   = args.wavelengthB
    timestart     = args.timestart
    timestep      = args.timestep
    timeend      = args.timeend
    Filter='Gaussian'#None#

    # time step between images
    DT=6.25 #s


    if timestart != '':
        timestartstr   = '{:0003d}'.format(int(timestart))
    else:
        timestartstr = ''
    if timeend != '':
        timeendstr   = '{:0003d}'.format(int(timeend))
    else:
        timeendstr = ''
    timestampBG    = args.timestampBG
    experimenttype  = args.experimenttype
    experimenttypeBG  = args.experimenttypeBG
    experimenttypesensitivity  = args.experimenttypesensitivity
    experimenttypesensitivityBG  = args.experimenttypesensitivityBG

    if PlotType=='Velocity':

        # Create an example plume intersection
        # Camera geometry
        DY=(3.828125+3.028125)/2.
        DZ=(3.903125+3.100000)/2.
        DY=0.250 # From Fig 2.
        DZ=0.125
        DS = np.sqrt(DZ*DZ+DY*DY)
        #        DS=1.4
        print "DS", DY, DZ, DS

        # Size of half image in the vertical, 5. = (96-86)/2.
        ZS=np.sin(np.deg2rad(5))*DS
        # 88/2 pixels in the vertical gives vertical pixel size
        dz = ZS / 44.0 #=0.0063961804401032104 #km =6.396180440m

        # From uvspec input file: mc_panorama_view 157.000000 203.000000 86.000000 96.000000
        # Size of half image in the horizontal 46/2
        XS=np.sin(np.deg2rad(23))*DS

        # 400/2 pixels in the horizontal gives horizontal pixel size
        DX = XS / 200. * 1000
        print "DX", DX

        case_here='CASEC'
        if case_here=='CASEA':
            COL_NUM1L = 200 #200 #180 #75 #150 #50 #
            COL_NUM1R = 200 #200 #180 #75 #150 #50 #
            COL_NUM2 =  175 #150 #180 #75 #150 #50 #
            COL_NUM3 =  225 #125 # 0 #180 #75 #150 #50 #
            ROW_BOTTOM = 0+150
            ROW_TOP    = 100+150
        elif case_here=='CASEB':
            COL_NUM1L = 100 #200 #180 #75 #150 #50 #
            COL_NUM1R = 100 #200 #180 #75 #150 #50 #
            COL_NUM2 =  125 #150 #180 #75 #150 #50 #
            COL_NUM3 =  150 #125 # 0 #180 #75 #150 #50 #
            ROW_BOTTOM = 0+150
            ROW_TOP    = 90+150
        elif case_here=='CASEC':
            COL_NUM1L = 330
            COL_NUM1R = 330
            COL_NUM2 = 300
            COL_NUM3 = 250
            ROW_BOTTOM = 0+150
            ROW_TOP    = 110+150

        # Clear sky ROI
        BG_ROI = [0, 10, 0, 40]
        # Setup optical flow class (you may play with the parameters)
        FLOW = pyplis.OptflowFarneback()
        
        PCS1 = pyplis.LineOnImage(COL_NUM1L, ROW_BOTTOM, COL_NUM1R, ROW_TOP,
                                 line_id="A", normal_orientation="right", #"left",#
                                 color="#e67300")
        PCS2 = pyplis.LineOnImage(COL_NUM2, ROW_BOTTOM, COL_NUM2,  ROW_TOP,
                                 line_id="B", normal_orientation="right", #"left",
                                 color="#1a1aff")
        PCS3 = pyplis.LineOnImage(COL_NUM3, ROW_BOTTOM, COL_NUM3,  ROW_TOP,
                                 line_id="C", normal_orientation="right", #"left",
                                 color="#1a1a00")

        pcs1 = PCS1
        pcs2 = PCS2
        pcs3 = PCS3
        # in this object the results from the optical flow histogram analysis
        # are stored (i.e. time series of local average displacement vector)
        props1 = pyplis.LocalPlumeProperties(roi_id=pcs1.line_id,
                                            color=pcs1.color)
        props2 = pyplis.LocalPlumeProperties(roi_id=pcs2.line_id,
                                            color=pcs2.color)
        props3 = pyplis.LocalPlumeProperties(roi_id=pcs3.line_id,
                                            color=pcs3.color)
        
        close("all")
        flow = FLOW

        ImagesAA = load_MYSTIC_images(filebase, wavelengthA,
                                      loc, timestart, timeend, timestep, experimenttype,
                                      wavelengthB, filebaseBG,
                                      timestampBG, experimenttypeBG, flip_y=True, Filter=Filter,
                                      AddRows=True)
        print('No. of images: {}'.format(len(ImagesAA)))
        
        NOISEAMPS =  [0, .4, .8, 2, 4]
        FLOW_WINS = np.arange(10, 111, 25)
        FLOW_ITER = [5, 10]
        FLOW_LEVS = [4, 6, 10]
        
        for noiseamp in NOISEAMPS:
            outdir = '{}NOISE{}/'.format(OUT_DIR, noiseamp)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            
            files = glob.glob('{}*.png'.format(outdir))
            for f in files:
                print('Deleting file {}'.format(f))
                os.remove(f)
            for i, im in enumerate(ImagesAA):
                im.img = im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape)
                # noisy = f + 0.4 * f.std() * np.random.random(f.shape)
                #            print im
                ax = im.show_img()
                ax.set_title("AA image {}".format(i+1))
                ax.figure.savefig('{}aa_img{}_N{}.png'.format(outdir, i,
                                  noiseamp))
                
                
                print(noiseamp)
            do_ana = True
            if do_ana:
                for flowwin in FLOW_WINS:
                    for flowiter in FLOW_ITER:
                        for flowlev in FLOW_LEVS:
                            figs = []
                            S = FLOW.settings
                            S._display["disp_skip"] = 4
                            S.winsize = flowwin
                            S.iterations = flowiter
                            S.levels = flowlev
                            S.hist_sigma_tol = 3
                            S.min_length = 1.5
                            S.hist_dir_gnum_max = 5
                            
                            T0 = datetime.datetime(2017,1,1,0,0,0) # define arbitrary start time
                            times = [T0 + datetime.timedelta(x * DT / 86400.0) for x in range(len(ImagesAA))]
                    
                            veff_arr1, veff_arr_err1 = [], []
                            veff_arr2, veff_arr_err2 = [], []
                            veff_arr3, veff_arr_err3 = [], []
                            lens = []
                            for k in range(len(ImagesAA)-1):
                                #            print("k", k)
                                this = ImagesAA[k]
                                this.edit_log["is_aa"]=True
                                this.meta["start_acq"]=times[k]
                                next_img = ImagesAA[k+1]
                                next_img.edit_log["is_aa"]=True
                                next_img.meta["start_acq"]=times[k+1]
                                flow.set_images(this, next_img)
                                #            print "times", times[k], times[k+1]
                                #            help(flow)
                                flow.calc_flow(this_img=this, next_img=next_img)
                                #            flow.calc_flow()
                                # fig = figure(figsize=(14,8))
                                # ax0 = fig.add_axes([0.01, 0.15, 0.59, 0.8])
                                # flow.plot(ax=ax0)
                                # pcs.plot_line_on_grid(ax=ax0, include_normal=1,
                                #                    include_roi_rot=1)
                                # plt.show()
                    
                                #            print("Optical flow", k)
                                sys.stdout.flush()
                    
                                lens.append(pyplis.Img(flow.get_flow_vector_length_img()).\
                                            crop(BG_ROI).mean())
                                # perform optical flow histogram analysis and append results
                                # to LocalPlumeProperties object
                                kwargs={'dir_multi_gauss': False}
                                props1.get_and_append_from_farneback(flow, line=pcs1, **kwargs)
                                props2.get_and_append_from_farneback(flow, line=pcs2, **kwargs)
                                props3.get_and_append_from_farneback(flow, line=pcs3, **kwargs)
                    
                                # calculate effective velocity through PCS for current image and
                                # based on results from histogram analysis
                                veff1, err1 = props1.get_velocity(pix_dist_m=DX,
                                                                  normal_vec=pcs1.normal_vector,
                                                                  sigma_tol=flow.settings.hist_sigma_tol)
                                veff2, err2 = props2.get_velocity(pix_dist_m=DX,
                                                                  normal_vec=pcs2.normal_vector,
                                                                  sigma_tol=flow.settings.hist_sigma_tol)
                                veff3, err3 = props3.get_velocity(pix_dist_m=DX,
                                                                  normal_vec=pcs3.normal_vector,
                                                                  sigma_tol=flow.settings.hist_sigma_tol)
                                print "k", k
                                print "veff1, err1", veff1, err1
                                print "veff2, err2", veff2, err2
                                print "veff3, err3", veff3, err3
                                sys.stdout.flush()
                    
                                veff_arr1.append(veff1)
                                veff_arr_err1.append(err1)
                                veff_arr2.append(veff2)
                                veff_arr_err2.append(err2)
                                veff_arr3.append(veff3)
                                veff_arr_err3.append(err3)
                            veff_arr1 = np.asarray(veff_arr1)
                            veff_arr_err1 = np.asarray(veff_arr_err1)
                            print "veff_arr1", veff_arr1
                            veff_arr2 = np.asarray(veff_arr2)
                            veff_arr_err2 = np.asarray(veff_arr_err2)
                            print "veff_arr2", veff_arr2
                            veff_arr3 = np.asarray(veff_arr3)
                            veff_arr_err3 = np.asarray(veff_arr_err3)
                            print "veff_arr3", veff_arr3
                            ax0 = this.show(zlabel="AA", zlabel_size=18)
    
                            fig0 = plot_histo_result(flow, [pcs1,pcs2,pcs3])
                            #        fig0 = plot_histo_result(flow, pcs2)
                            fig, axes = subplots(2,1, figsize=(8,10))#, sharex=True)
                            props1.plot_velocities(pix_dist_m=DX, pix_dist_m_err=1,
                                                  normal_vec=pcs1.normal_vector, color=pcs1.color,
                                                  date_fmt="%H:%M:%S",
                                                   ax=axes[1])
                            props1.plot_directions(ax=axes[0])
                            props2.plot_velocities(pix_dist_m=DX, pix_dist_m_err=1,
                                                  normal_vec=pcs2.normal_vector, color=pcs2.color,
                                                  date_fmt="%H:%M:%S",
                                                   ax=axes[1])
                            props2.plot_directions(ax=axes[0])
                            props3.plot_velocities(pix_dist_m=DX, pix_dist_m_err=1,
                                                  normal_vec=pcs3.normal_vector, color=pcs3.color,
                                                  date_fmt="%H:%M:%S",
                                                   ax=axes[1])
                            props3.plot_directions(ax=axes[0])
                            axes[1].set_ylabel(r"$v_{eff}$ [m/s]")
                            axes[1].set_ylabel(r"$v_{eff}$ [m/s]")
                            axes[0].set_xticks([])
                            axes[0].set_title("Time series of displacement orientation angle")
                            img=flow.get_flow_vector_length_img()
                            figs.extend([fig0, fig])
                    
                            #        veff_str = "veff (avg) = %.2f +/- %.2f" %(veff_arr1.mean(),
                            #                      veff_arr_err1.mean())
                            veff_str = "veff$_{young}$(avg) = %.2f +/- %.2f, veff$_{old}$(avg) = %.2f +/- %.2f" %(veff_arr1.mean(),
                                                                                                                  veff_arr_err1.mean(),
                                                                                                                  veff_arr2.mean(),
                                                                                                                  veff_arr_err2.mean())
                            print "Max (mean) displ. length in BG ROI: %.1f" %max(lens)
                            print veff_str
                            axes[1].set_title(veff_str)
                            SAVE_DIR = outdir
                            SAVE_PLOT_TYPE = "png"
                            SAVE_PLOT_DPI = 100
                    
                            SCRIPT_ID = basename(__file__).split("_")[0]
                            for k, fig in enumerate(figs):
                                s = ('{}{}_win{}_iter{}_levs{}_{}.{}'
                                     .format(outdir, case_here, flowwin, flowiter, 
                                             flowlev, (k+1), SAVE_PLOT_TYPE))
                                fig.savefig(s, dpi=SAVE_PLOT_DPI)
                
        exit()

    elif PlotType=='Fractaldim':

        ImagesAA = load_MYSTIC_images(filebase, wavelengthA,
                                      loc, timestart, timeend, timestep, experimenttype,
                                      wavelengthB, filebaseBG,
                                      timestampBG, experimenttypeBG, flip_y=True, Filter=Filter,
                                      AddRows=True)

        #        ths =  np.arange(5,100)/255. #245,5)
        ths = np.array([0.001])
        indx = np.where((ths>35) & (ths<50) )
        fdim=[]
        fdimstd=[]

        for k in range(len(ImagesAA)):
            fds = []
            this = ImagesAA[k]
            this.edit_log["is_aa"]=True
            img = np.where(this.img<=0.0, 1e-10, this.img)
            print "A: min,max ", np.min(img), np.max(img), img.shape
#            img = np.exp(img)
            print "B: min,max ", np.min(img), np.max(img), img.shape
            a = 1/(np.max(img)-np.min(img))
            #            a = 255/(np.max(img)-np.min(img))
            #            a = 255/(200-np.min(img))
            img = img*a-a*np.min(img)
            print "C: min,max ", np.min(img), np.max(img), img.shape

            plothist=False #True#
            if plothist:
                fig, axs = plt.subplots(nrows=2, ncols=1)
                ax = axs[0]
                nxs, xbins, ptchs = ax.hist(img.flatten(), 50, color='g', histtype='bar')
                ax.set_yscale('log', nonposy='clip')
                ax = axs[1]
                im = ax.imshow(img, origin='lower')
                plt.show()

            plotdd=True#False#
            if plotdd:
                fig, axfd = plt.subplots(1, 1, figsize=(6,6))
            else:
                axfd=None

            for th in ths:
                plotimg=False#True
                if plotimg:
                    figimg, aximg = plt.subplots(1, 1, figsize=(5,10))
                fd=fractal_dimension(img, threshold=th, ax=axfd)
                fds.append(fd)
                if plotimg:
                    sizes = np.array([64, 32, 24])#, 16, 14, 12, 10, 8, 6, 4, 2])
                    Z=img
                    xmin, xmax = 0, Z.shape[1]
                    ymin, ymax = 0, Z.shape[0]
                    for i, size in enumerate(sizes):
                        ax = plt.subplot(len(sizes), 1, i+1, frameon=False)
                        ax.imshow(1-Z, plt.cm.gray, interpolation="bicubic", vmin=0, vmax=1,
                                  extent=[xmin, xmax, ymin, ymax], origin="upper")
                        ax.set_xticks([])
                        ax.set_yticks([])
                        for y in range(Z.shape[0]//size+1):
                            for x in range(Z.shape[1]//size+1):
                                s = (Z[y*size:(y+1)*size, x*size:(x+1)*size] > th).sum()
                                if s > 0 and s < size*size:
                                    rect = patches.Rectangle(
                                        (x*size, Z.shape[0]-1-(y+1)*size),
                                        width=size, height=size,
                                        linewidth=.5, edgecolor='.25',
                                        facecolor='.75', alpha=.5)
                                    ax.add_patch(rect)

                    plt.tight_layout()
                    plt.show()

            fds=np.array(fds)
            if plotdd:
                axfd.set_xlabel(r"$\log(sizes)$", fontsize=20)
                axfd.set_ylabel(r"$\log(counts)$", fontsize=20)
                axfd.text(3, 7,'$D=${:5.3f}'.format(fds[0]), fontsize =FONTSIZE, color='black')
                xmin=0.
                xmax=6
                axfd.set_xlim(xmin,xmax)
                ymin=3
                ymax=8
                axfd.set_ylim(ymin,ymax)
                plt.show()

            print "fds", fds
            indx = np.where(fds>1.0 )

            fdsmean = np.mean(fds[indx])
            fdsstd = np.std(fds[indx])
            print "fds mean std", np.mean(fds[indx]), fdsstd
            fdim.append(fdsmean) #fds[indx])
            fdimstd.append(fdsstd) #fds[indx])

            plotfd=False#True#
            if plotfd:
                fig, ax = plt.subplots(1, 1, figsize=(6,6))
                p=ax.plot(ths, fds, linestyle='-')
                p=ax.plot([ths[0],ths[-1]], [fdsmean,fdsmean] , linestyle=':')
                plt.show()
#                exit()

        fdim=np.array(fdim)
        fdimstd=np.array(fdimstd)
        xs = np.arange(0,len(ImagesAA))
        for x, fd, fdstd in zip(xs, fdim, fdimstd):
            print '{:3d} {:6.4f} {:6.4f}\n'.format(x, fd, fdstd)
            #            fp.write( '{:3d} {:6.4f} {:6.4f}\n'.format(x, fd, fdstd))


        LESfdx, LESfd, LESfdstd = np.loadtxt('./LES_Fractal_Dimension_V02.txt', unpack=True)

        fig, axfdim = plt.subplots(1, 1, figsize=(6,6))
        p=axfdim.errorbar(xs, fdim, yerr=fdimstd, fmt='o')
        p=axfdim.errorbar(LESfdx, LESfd, yerr=LESfdstd, fmt='v')
        #        p=ax.plot([ths[0],ths[-1]], [fdsmean,fdsmean] , linestyle=':')
        plt.show()

        exit()

    elif PlotType=='OnBandRadiance':
        uvinp = 'uvspecCamW{:s}{:s}{:s}{:s}.inp'.format(wavelengthA,loc,timestartstr,experimenttype)
        uvout = 'uvspecCamW{:s}{:s}{:s}{:s}.out'.format(wavelengthA,loc,timestartstr,experimenttype)
        uvspecinputfile = filebase+uvinp
        uvspecoutputfile = filebase+uvout
        ImgAM = Image.Image()
        ImgAM.GetUVSPECImageInputVals(uvspecinputfile)
        ImgAM.ReadMYSTICradspcfile(uvspecoutputfile, AddRows=False)
        start_acq=datetime.datetime(2016, 10, 10, 13, 15, 12)
        MYSTIC_AM = Img(start_acq = start_acq)
        MYSTIC_AM.img=ImgAM.rad
        MYSTIC_AM.meta["bit_depth"]=12  # Just to get show_histogram working.
        MYSTIC_AM.meta["pix_heigth"]=ImgAM.rad.shape[0]
        MYSTIC_AM.meta["pix_width"]=ImgAM.rad.shape[1]

    else:
        MYSTIC_AM, MYSTIC_A0, MYSTIC_BM, MYSTIC_B0 = load_MYSTIC_image(filebase, wavelengthA,
                                                                        loc, timestartstr, experimenttype,
                                                                        wavelengthB, filebaseBG,
                                                                        timestampBG, experimenttypeBG, Filter=Filter)

        tmpImg = copy.deepcopy(MYSTIC_A0)  # Just to angles etc. included

        tau_A = MYSTIC_AM.to_tau(MYSTIC_A0) # Note that Jonas definition of tau differs
                                            # with a minus sign from Lubcke et al. 2013 definition.
        tau_B = MYSTIC_BM.to_tau(MYSTIC_B0)
        tau_A.img = tau_A.img-tau_B.img
        # Zero all values below threshold
        tau_A.threshold=0.03
        tau_A.set_val_below_thresh(val=0, threshold=tau_A.threshold)
        tau_A.statistics()
        tau_A.meta["phi1"] = tmpImg.meta["phi1"]
        tau_A.meta["phi2"] = tmpImg.meta["phi2"]

        if experimenttypesensitivity != '':
            MYSTIC_AMS, MYSTIC_A0S, MYSTIC_BMS, MYSTIC_B0S = load_MYSTIC_image(filebase, wavelengthA,
                                                                            loc, timestartstr, experimenttypesensitivity,
                                                                            wavelengthB, filebaseBG,
                                                                            timestampBG, experimenttypesensitivityBG, Filter=Filter)

            tau_AS = MYSTIC_AMS.to_tau(MYSTIC_A0S)
            tau_BS = MYSTIC_BMS.to_tau(MYSTIC_B0S)
            tau_AS.img = tau_AS.img-tau_BS.img
            # Zero all values below threshold
            tau_AS.threshold=0.03
            tau_AS.set_val_below_thresh(val=0, threshold=tau_A.threshold)
            tau_AS.statistics()


    if PlotType=='Statistics':
        includeLES=True# False#
        if includeLES:
            path = '../Experiments/palm_tomo/'
            fns =sorted(glob(path+'LES_25_'+str(timestep)+'.nc'))
            LESdata=LES2Img(path, fns)
        else:
            LESdata=None

        kwargs = {'title': '',
                  'xlabel': 'Horizontal pixel number',
                  'ylabel': 'Vertical pixel number',
                  'cmap': cmaps.Greys,
                  'log_img': True,
                  'threshold':  tau_A.threshold,
                  'extent': [0,MYSTIC_AM.meta["pix_width"],0,MYSTIC_AM.meta["pix_heigth"]],
                  'data': tau_A.img,
                  'LESdata': LESdata
              }

        plot_img_statistics(tau_A, **kwargs)

        if pngfolder != '':
            pngfile=pngfolder+experimenttype+timestartstr+'_stats_V1.png'
            print("pngfile", pngfile)
            plt.savefig(pngfile)
        else:
            plt.show(block=True)

        exit()

    elif PlotType=='StatisticsDiff':
        tau_diff = tau_A.diff(tau_AS)
        kwargs = {'title': 'a) On-band radiance (mW m$^{-2}$ nm$^{-1}$)',
                  'xlabel': 'Horizontal pixel number',
                  'ylabel': 'Vertical pixel number',
                  'cmap': cmaps.Greys,
                  'log_img': False,
                  'threshold': 0.0,
                  'extent': [0,MYSTIC_AM.meta["pix_width"],0,MYSTIC_AM.meta["pix_heigth"]]
              }

        kwargs['title']='a) Difference in apparent absorbance'
        plot_img_statistics_5(tau_diff, **kwargs)

        plt.show()
        exit()

    elif PlotType=='ImageDiff':
        kwargs = {'tau_A': tau_AS }
        tau_diff = tau_A.diff(tau_AS)
        # Minus to put standard at end.
        tau_diff.img = -tau_diff.img
        print_img_statistics(tau_diff, **kwargs)

        print "GABBA", experimenttypesensitivity, experimenttypeBG
        if 'AEROATMA' in experimenttypesensitivity:
            #            cb_ticks= [-0.5, -0.25, 0.0, 0.25, 0.5]
            cb_ticks= [-0.04, -0.02, 0.0, 0.02, 0.04]
            title=r'a) Difference in apparent absorbance, "background aerosol" $-$ "no aerosol"'
            png_name = 'SUNSZA40W_AEROATMA_SUNSZA40W_005_AADiff.png'
            vmin=-0.04#-0.5
            vmax=0.04#0.5
        elif 'AEROPLU' in experimenttypesensitivity:
            #            cb_ticks= [-0.5, -0.25, 0.0, 0.25, 0.5]
            if 'AEROPLUA' in experimenttypesensitivity:
                title=r'b) Difference in apparent absorbance, $\tau(\sim 0.5)-\tau(0.0)$'
                png_name = 'SUNSZA40W_AEROPLUA_SUNSZA40W_005_AADiff.png'
            elif 'AEROPLUC' in experimenttypesensitivity:
                title=r'c) Difference in apparent absorbance, $\tau(\sim 5.0)-\tau(0.0)$'
                png_name = 'SUNSZA40W_AEROPLUC_SUNSZA40W_005_AADiff.png'
            elif 'AEROPLUB' in experimenttypesensitivity:
                title=r'c) Difference in apparent absorbance, $\tau(\sim 50.0)-\tau(0.0)$'
                png_name = 'SUNSZA40W_AEROPLUB_SUNSZA40W_005_AADiff.png'

            fact=1.
            cb_ticks= np.array([-0.04, -0.02, 0.0, 0.02, 0.04])*fact
            vmin=-0.04*fact#-0.5
            vmax=0.04*fact#0.5
        elif 'ALBEDO' in experimenttypesensitivity or 'ALBEDO' in experimenttypeBG:
            #            cb_ticks= [-2.0, -1.5, -1.0, -0.5, 0.0]
            cb_ticks= [-2.0,  -1.0,  0.0, 1, 2]
            if 'SUNSZA40WQBLC' in experimenttypesensitivity:
                title=r'a) Difference in apparent absorbance, $A(0.0)-A(0.05)$'
                png_name = 'SUNSZA40W_ALBEDOC_SUNSZA40W_005_AADiff.png'
            elif 'SUNSZA40WALBEDOB' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $A(0.1)-A(0.05)$'
                png_name = 'SUNSZA40W_ALBEDOB_SUNSZA40W_005_AADiff.png'
            elif 'SUNSZA40WALBEDOA' in experimenttypesensitivity:
                title=r'c) Difference in apparent absorbance, $A(1.0)-A(0.05)$'
                png_name = 'SUNSZA40W_ALBEDOA_SUNSZA40W_005_AADiff.png'
            vmin=-2.0#-0.5
            vmax=2.0#0.5
            vmin=-0.04#-0.5
            vmax=0.04#0.5
        elif 'SUNSZA40' in experimenttypesensitivity:
            cb_ticks= [-0.04, -0.02, 0.0, 0.02, 0.04]
            if 'SUNSZA40NWQBLC' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $\tau(\phi_0=90)-\tau(\phi_0=45)$'
                png_name = 'SUNSZA40WQBLC_SUNSZA40NWQBLC_005_AADiff.png'
            elif 'SUNSZA40NQBLC' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $\tau(\phi_0=90)-\tau(\phi_0=0)$'
                png_name = 'SUNSZA40WQBLC_SUNSZA40NQBLC_005_AADiff.png'
            elif 'SUNSZA40SWQBLC' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $\tau(\phi_0=90)-\tau(\phi_0=135)$'
                png_name = 'SUNSZA40WQBLC_SUNSZA40SWQBLC_005_AADiff.png'
            elif 'SUNSZA40SQBLC' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $\tau(\phi_0=90)-\tau(\phi_0=180)$'
                png_name = 'SUNSZA40WQBLC_SUNSZA40SQBLC_005_AADiff.png'
            elif 'SUNSZA40EQBLC' in experimenttypesensitivity:
                title=r'Difference in apparent absorbance, $\tau(\phi_0=90)-\tau(\phi_0=270)$'
                png_name = 'SUNSZA40WQBLC_SUNSZA40EQBLC_005_AADiff.png'
            vmin=-0.04#-0.5
            vmax=0.04#0.5
        elif 'SUNSZA60W' in experimenttypesensitivity:
            cb_ticks= [-0.04, -0.02, 0.0, 0.02, 0.04]
            title=r'Difference in apparent absorbance, $\tau(\theta_0=60)-\tau(\theta_0=40)$'
            png_name = 'SUNSZA40W_SUNSZA60W_005_AADiff.png'
            vmin=-0.04#-0.5
            vmax=0.04#0.5
        else:
            cb_ticks= [-0.15,-0.1,-0.05,0.0, 0.05, 0.1, 0.15]
            title='a) Difference in apparent absorbance'
            png_name = 'tmp.png'#
            vmin=-0.15
            vmax=0.15
        kwargs = {'xlabel': 'Horizontal pixel number',
                  'ylabel': 'Vertical pixel number',
                  'cmap': cmaps.seismic, #cmaps.spring, #cmaps.Greys,#
                  'log_img': False,
                  'cb_shrink':0.5,
                  'cb_ticks': cb_ticks,
                  'flip_x': False,
                  'figsize':(18,5.2),
                  'threshold':  -9999,
                  'extent': [0,MYSTIC_AM.meta["pix_width"],0,MYSTIC_AM.meta["pix_heigth"]]
              }

        print "experimenttypesensitivity", experimenttypesensitivity

        kwargs['title']=title
        kwargs['vmax'] = vmax
        kwargs['vmin'] = vmin

        print("diff min/max", tau_diff.img.min(), tau_diff.img.max())
        plot_img(tau_diff.img, **kwargs)
        if pngfolder != '':
            pngfile=pngfolder+png_name
            print("pngfile", pngfile)
            plt.savefig(pngfile)
        else:
            plt.show()

        exit()

    elif PlotType=='OnBandRadiance':
        # Do plotting outside pyplis to get more control of plot.
        kwargs = {'title': 'a) On-band radiance (mW m$^{-2}$ nm$^{-1}$)',
                  'xlabel': 'Horizontal pixel number',
                  'ylabel': 'Vertical pixel number',
                  'cmap': cmaps.Greys_r,
                  'log_img': False,
                  'flip_x': False,
                  'threshold':  0.0,
#                  'extent': [0,MYSTIC_AM.meta["pix_width"],MYSTIC_AM.meta["pix_heigth"],0]
              }

        plot_img(MYSTIC_AM.img, **kwargs)
        plt.show()
        exit()

    elif PlotType=='Images':

        # Do plotting outside pyplis to get more control of plot.
        kwargs = {'title': 'a) On-band radiance (mW m$^{-2}$ nm$^{-1}$)',
                  'xlabel': 'Horizontal pixel number',
                  'ylabel': 'Vertical pixel number',
                  'cmap': cmaps.Greys_r,
                  'flip_x': False,
                  'log_img': False,
                  'threshold':  tau_A.threshold,
#                  'extent': [0,MYSTIC_AM.meta["pix_width"],MYSTIC_AM.meta["pix_heigth"],0]
              }

        plot_img(MYSTIC_AM.img, **kwargs)

        kwargs['title']= 'b) Off-band radiance (mW m$^{-2}$ nm$^{-1}$)'
        plot_img(MYSTIC_BM.img, **kwargs)

        kwargs['title']='c) Apparent absorbance, linear scale'
        kwargs['cmap']= cmaps.Greys
        #        kwargs['vmin']= 0.0
        #        kwargs['vmax']= 0.2
        plot_img(tau_A.img, **kwargs)

        kwargs['log_img']= True
        kwargs['title']='d) Apparent absorbance, log$_{10}$ scale'
        plot_img(tau_A.img, **kwargs)

        plt.show()
        exit()

    elif PlotType=='radiances_and_absorbances':
        kwargs1 = {'title': 'a) On-band radiance (mW m$^{-2}$ nm$^{-1}$)',
                   'xlabel': 'Horizontal pixel number',
                   'ylabel': 'Vertical pixel number',
                   'cmap': cmaps.Greys_r,
                   'log_img': False,
                   'threshold':  tau_A.threshold,
                   'extent': [0,MYSTIC_AM.meta["pix_width"],0,MYSTIC_AM.meta["pix_heigth"]],
                   'data': MYSTIC_AM.img
              }
        kwargs2 = {'title': 'b) Off-band radiance (mW m$^{-2}$ nm$^{-1}$)',
                   'xlabel': 'Horizontal pixel number',
                   'ylabel': 'Vertical pixel number',
                   'cmap': cmaps.Greys_r,
                   'log_img': False,
                   'threshold':  tau_A.threshold,
                   'extent': [0,MYSTIC_BM.meta["pix_width"],0,MYSTIC_BM.meta["pix_heigth"]],
                   'data': MYSTIC_BM.img
              }
        kwargs3 = {'title': 'c) Apparent absorbance, linear scale',
                   'xlabel': 'Horizontal pixel number',
                   'ylabel': 'Vertical pixel number',
                   'cmap': cmaps.Greys,
                   'log_img': False,
                   'vmin' :  0.0,
                   'vmax':  1.0,
                   'threshold':  tau_A.threshold,
                   'extent': [0,MYSTIC_BM.meta["pix_width"],0,MYSTIC_BM.meta["pix_heigth"]],
                   'data': tau_A.img
              }
        kwargs4 = {'title': 'd) Apparent absorbance, log$_{10}$ scale',
                   'xlabel': 'Horizontal pixel number',
                   'ylabel': 'Vertical pixel number',
                   'cmap': cmaps.Greys,
                   'log_img': True,
                   'threshold':  tau_A.threshold,
                   'extent': [0,MYSTIC_BM.meta["pix_width"],0,MYSTIC_BM.meta["pix_heigth"]],
                   'data': tau_A.img
              }

        kwargsall = {'kwargs1': kwargs1,
                     'kwargs2': kwargs2,
                     'kwargs3': kwargs3,
                     'kwargs4': kwargs4,
                 }
        plot_images_column(**kwargsall)

        if pngfolder != '':
            pngfile=pngfolder+experimenttype+timestartstr+'_rad_abs_V1.png'
            print("pngfile", pngfile)
            plt.savefig(pngfile)
        else:
            plt.show(block=True)

    elif PlotType=='radiances_and_absorbances_NOANNO':
        fig, ax = plt.subplots(1, 1, figsize=(15,4.5))
        cmap = cmaps.Greys
        vmin = 0.0
        vmax = 1.0
        img = tau_A.img
        extent = [0,MYSTIC_BM.meta["pix_width"],0,MYSTIC_BM.meta["pix_heigth"]]
        print extent
#        im = ax.imshow(img,cmap=cmap,extent=extent,origin='lower',vmin=vmin,vmax=vmax,interpolation='bilinear' )
        im = ax.imshow(img,cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,interpolation='bilinear' )
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.axis('off')

        if pngfolder != '':
            pngfile=pngfolder+experimenttype+timestartstr+'_rad_abs_NOANNO_V1.png'
            print("pngfile", pngfile)
            plt.savefig(pngfile)
        else:
            plt.show(block=True)
        exit()
