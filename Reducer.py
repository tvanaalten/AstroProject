from __future__ import print_function
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.io.fits import HDUList
from astropy.io.fits import PrimaryHDU
from astropy.io.fits import getval,setval

import os
import Constants
    
class Reducer:
    
    objdir = None
    mydir = None
    image_names = None
    
    fil = None
        
    bias = None
    flatfield = None
    set_size = None
    n_sets = 0
    
    
    def __init__(self, dir, filter, image_names, bias_file, flat_file):
        self.objdir = dir
        self.mydir = os.path.expanduser('~/')
        self.image_names = image_names
        self.fil = filter
        self.get_bias_and_flatfield(bias_file, flat_file)
        
    def get_bias_and_flatfield(self, bias_file, flat_file):
        
        self.bias=getdata(self.objdir+bias_file)

        self.flat=getdata(self.objdir+flat_file)

        
        
    def reduce(self, has_sets):
        newdir = self.objdir + Constants.working_directory
        
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        
        newdir = newdir + Constants.image_directory
        
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        
        strlen=len(self.image_names)
        
        
        i = 0        
        
        for file in os.listdir(self.objdir):
            print(file)
            if file[:strlen]==self.image_names:
                filter=getval(self.objdir+file,'FILTER',ignore_missing_end=True)
                print(filter, self.fil)
                if filter[:1]==self.fil or filter == self.fil:
                    print(file,filter[:1]) 
                    data=(getdata(self.objdir+file,ignore_missing_end=True) - self.bias) / self.flat
                    head=getheader(self.objdir+file,ignore_missing_end=True)
                    hdu = PrimaryHDU(data, head)
                    hdul = HDUList([hdu], None)
                    
                    
                    filepath = newdir + Constants.reduced_prefix + self.image_names
                    
                    if(has_sets):
                        
                        array = file.split("_")
                        set = array[1]
                        i = array[2].split(".")[0]
                        
                        if i > self.set_size:
                            self.set_size = i
                        
                        if set > n_sets:
                            n_sets = set
                        
                        print(set)
                        print(i)
                        
                        filepath += "_" + set + "_" + i
                    
                    else:
                        i+=1
                        n_length = len(str(i))
                        
                        for j in range(3-n_length):
                            filepath += "0"
                        
                        filepath += i
                            
                    filepath += Constants.fits_extension


                    
                    
                    print(newdir)
                    print(filepath)
                    
                    hdul.writeto(filepath, overwrite=True)
                    print(i)
                    
    def get_set_info(self):
        return self.set_size, self,n_sets
                    
            
        
    
    
