import math
import os
import sys
import numpy as np
import multiprocessing
from subprocess import Popen,PIPE, STDOUT, call
from glob import glob
import time

# Where are we
myhost = ''#os.uname()[1]
if myhost=='ak':
    home = os.environ['HOME']+'/'
else:
    home = '/xnilu_wrk/aky/'

class UVspec:
    def __init__(self):
        # Set some uvspec input that most likely will stay the same
        self.inp = {
            "data_files_path"    : home+'develop/libRadtran/data/',
            "quiet"              : '',
            }
        self.inp["rte_solver"]  = 'montecarlo'
        self.inp["mc_backward"] = ''
        self.inp["atmosphere_file"] = home+'develop/libRadtran/data/atmmod/afglms.dat'
        self.inp["source"] = 'solar '+home+'develop/libRadtran/data/solar_flux/atlas_plus_modtran'
        self.inp["mc_vroom"] = 'on'
        self.inp["albedo"] = '0.0'

        return

    def WriteInputFile(self, InputFile=None, verbose=False):
        self.InputFile=InputFile
        if verbose:
            print("Writing uvspec input file", InputFile)
        try:
            f = open(InputFile,'w')
        except:
            print("UVspec.WriteInputFile: No uvspec input file name given.")
            exit()
        for key in self.inp:
            if verbose:
                sys.stdout.write( key + ' ' + str(self.inp[key]) + '\n')
            f.write( key + ' ' + str(self.inp[key]) + '\n')
        f.close()
        return

    def add_mc_basename_to_input_file(self,mc_basename,fn):
        f = open(fn,'a')
        f.write('{0:s}\n'.format('mc_basename '+mc_basename))
        f.close()

    def worker(self, input,output):
        """thread worker function"""
        verbose = 0
        self.SingleRun(input,output,verbose)
        return

    def SingleRun(self,inp, out, verbose):
        if verbose:
            print("Running uvspec with input file: ", inp)
            print("Output to file                : ", out)

        cmd = home+'/develop/libRadtran/bin/uvspec '+  ' < ' + inp  +  ' > ' + out
        p   = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
        return

    def Run(self,inp, out, verbose=False, n_processes=1):
        debug=False # True
        if verbose:
            print("Running uvspec with input file: ", inp)
            print("Output to file                : ", out)
            print("Number of processors          : ", n_processes)

        tmp_out_base = 'tmp_mystic.out_'
        tmp_inp_base = 'tmp_mystic.inp_'
        # Remove all old files
        for filename in glob(tmp_out_base+"*"):
            if not debug:
                os.remove(filename)
        for filename in glob(tmp_inp_base+"*"):
            if not debug:
                os.remove(filename)

        jobs = []
        for i in range(n_processes):
            # Copy input file to temporary input file to be able to add different
            # mc_basenames to the file without destroying the input file
            tmp_inp = tmp_inp_base+str(i)
            cmd = 'cp '+inp+' '+tmp_inp
            Popen([r"cp",inp, tmp_inp]).wait()
            mc_basename = tmp_out_base+'NP_'+str(i)
            self.add_mc_basename_to_input_file(mc_basename,tmp_inp)
            tmp_out = tmp_out_base+str(i)
            if verbose:
                print('Starting process:',i,' inp:',tmp_inp,' out:',tmp_out)
            if not debug:
                p = multiprocessing.Process(target=self.worker, 
                                            args=(tmp_inp,tmp_out))
                jobs.append(p)
                p.start()
        for j in jobs:
            j.join()

        if verbose:
            print('All processes done. Read output, average and calculate std.')
        InputFiles = tmp_out_base+'NP_'+'*'+'.rad.spc'
        OutputFile= out
        Average_spc_Files(InputFiles, OutputFile, verbose=True)


def zenith(lat, lon, year, month, day,hour,min=0, sec=0, stdlong=0,output=1, uvspecpath=''):

    cmd = uvspecpath+'zenith '+str(day)+' '+str(month)+' '+str(hour)+' '+str(min)+' '+str(sec)+\
          ' '+'-y '+str(year)+' -a '+str(lat)+' -o '+str(lon)+' -s '+str(stdlong)+' -q'
    res = Popen(cmd,shell=True,stdout=PIPE)
    res.wait()
    vals = res.communicate()
    vals = vals[0].split()
    sza = float(vals[1])
    return sza


def get_vals(fn,option):
    """ Returns the values for option in an input file.

        Usage:

        values = get_vals(input_filename,optionname)

        Input:
           filename    uvspec input file
           optionname  name of uvspec option

        Output:
           values      list of option values

        Author: Arve Kylling
        Date:   2011-05-23
    """

    f  = open(fn,'r')
    vals = ''
    for line in f:
        l = line.split()
# This does not work with the new input options.....
#        if ( l[0] == option ):
#            vals = l[1:len(l)]
#        print(l, option
        if option in line:
            nopts = len(option.split())
            vals = l[nopts:len(l)]
#            print(l, option, nopts, vals

    f.close()
    return vals


def Average_spc_Files(InputFiles, OutputFile, verbose=False):
    # First check that all files have the same number of lines. If not
    # the files are surely different.
    i = 0
    for fn in glob(InputFiles):
        with open(fn) as fp:
            nlin = sum(1 for line in fp)
            if i==0:
                nlin0=nlin
            else:
                if nlin != nlin0:
                    print('nlin: ' + str(nlin) + ', not equal nlin0: ' 
                          + str(nlin0))
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
    for f in  glob(InputFiles):
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
            print(sys.stderr, l, arg, s0, s1, s2[l])
            arg = 0.0
        std       = (1.0/s0)*math.sqrt(arg)
        f.write('{0:8.2f} {1:3d} {2:3d} {3:3d} {4:9.4f} {5:9.4f}\n'.format(wvl[0,l], ix[0,l], iy[0,l], iz[0,l], s1/s0, std))
        l = l + 1
    f.close()
    return

def read_rad_spc(fn, STD=False, verbose=False):
    # Read MYSTIC mc.rad.spc file
    if verbose:
        print("Reading MYSTIC mc.rad.spc file: ", fn)
    if STD:
        wvl,ix,iy,iz,rad, std = np.loadtxt(fn, unpack=True)
        return (wvl,ix,iy,iz,rad,std)
    else:
        wvl,ix,iy,iz,rad = np.loadtxt(fn, unpack=True)
        return (wvl,ix,iy,iz,rad)
