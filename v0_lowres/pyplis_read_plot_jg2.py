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
import matplotlib.cm as cmaps
from matplotlib.legend import rcParams

import argparse
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
    import os
    plt.close("all")
    
    OUT_DIR = './output/'

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
            COL_NUM1L = 250
            COL_NUM1R = 250
            COL_NUM2 = 290
            COL_NUM3 = 300
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
        

        flow = FLOW

        ImagesAA = load_MYSTIC_images(filebase, wavelengthA,
                                      loc, timestart, timeend, timestep, experimenttype,
                                      wavelengthB, filebaseBG,
                                      timestampBG, experimenttypeBG, flip_y=True, Filter=Filter,
                                           parser = argparse.ArgumentParser(description='Plot camera images and apparant absorbance. '+
#                                  'Script requires that both foreground and background simulations'+
#                                  'are available in folder specified by InputFolder.')
# 
# 
#     parser.add_argument('InputFolder', type=str,
#                         help='Full folder name holding experiment')
#     parser.add_argument('--PlotType', type=str, default='Images',
#                         help='Type of plot/analysis to be made')
#     parser.add_argument('--InputFolderBG', type=str, default='',
#                         help='Full folder name holding experiment')
#     parser.add_argument('Camera', type=str,
#                         help='Camera location to be plotted')
#     parser.add_argument("--timestart", help="Start time stamp of LES simulation", default='')
#     parser.add_argument("--timeend", help="End time stamp of LES simulation", default='')
#     parser.add_argument("--timestep", help="Time step in stamp units of LES simulation", default='1')
#     parser.add_argument("--timestampBG", help="Time stamp of BG simulation", default='')
#     parser.add_argument("--wavelengthA", help="wavelength name of camera A", default='310')
#     parser.add_argument("--wavelengthB", help="wavelength name of camera B", default='330')
#     parser.add_argument('--experimenttype', type=str, default='',
#                         help='String defining type of experiment. See code for more info.')
#     parser.add_argument('--experimenttypeBG', type=str, default='-999',
#                         help='String defining type of BG experiment. See code for more info.')
#     parser.add_argument('--experimenttypesensitivity', type=str, default='',
#                         help='String defining type of sensitivity experiment. See code for more info.')
#     parser.add_argument('--experimenttypesensitivityBG', type=str, default='',
#                         help='String defining type of BG sensitivity experiment. See code for more info.')
#     parser.add_argument("--pngfolder", help="Name of png file folder", default='')
# 
#     args = parser.parse_args()
#     filebase      = args.InputFolder
#     filebaseBG    = args.InputFolderBG
#     PlotType      = args.PlotType
#     loc           = args.Camera
#     pngfolder     = args.pngfolder
#     wavelengthA   = args.wavelengthA
#     wavelengthB   = args.wavelengthB
#     timestart     = args.timestart
#     timestep      = args.timestep
#     timeend      = args.timeend
#     Filter='Gaussian'#AddRows=True)
        print('No. of images: {}'.format(len(ImagesAA)))
        
        NOISEAMPS =  [0, .4, .8, 2, 4]
        FLOW_WINS = np.arange(10, 111, 25)
        FLOW_ITER = [5, 10]
        FLOW_LEVS = [4, 6, 10]
        
        noiseamp = NOISEAMPS[0]
        outdir = OUT_DIR
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        files = glob.glob('{}*.png'.format(outdir))
        
        ts_pcs1 = []
        ts_pcs2 = []
        for i, im in enumerate(ImagesAA):
            im.img = im.img  + noiseamp * im.img.std() * np.random.random(im.img.shape)
            # noisy = f + 0.4 * f.std() * np.random.random(f.shape)
            #            print im
            ax = im.show_img()
            ax.set_title("AA image {}".format(i+1))
            
            pcs1.plot_line_on_grid(ax=ax)
            pcs2.plot_line_on_grid(ax=ax)
            ax.legend()
            ax.figure.savefig('{}aa_img{}_N{}.png'.format(outdir, i,
                              noiseamp))
            
            
            
            ts_pcs1.append(sum(pcs1.get_line_profile(im.img)))
            
            ts_pcs2.append(sum(pcs2.get_line_profile(im.img)))
            
            
            im = im.pyr_down(3)
            ax = im.show_img()
            ax.set_title("AA image (small) {}".format(i+1))
            
            pcs1.plot_line_on_grid(ax=ax)
            pcs2.plot_line_on_grid(ax=ax)
            ax.legend()
            ax.figure.savefig('{}small_aa_img{}_N{}.png'.format(outdir, i,
                              noiseamp))
        fig, ax = plt.subplots(1,1)
        ax.plot(ts_pcs1, label=pcs1.line_id, color=pcs1.color)
        ax.plot(ts_pcs2, label=pcs2.line_id, color=pcs2.color)   
        ax.legend()
        ax.set_ylabel('ICA (sum)')
        ax.set_xlabel('Img number')
        fig.savefig('{}PCS_timeseries.png'.format(outdir))        