from photutils import CircularAperture
from photutils import DAOStarFinder
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import os
import pandas
from tabulate import tabulate
import Constants
import Utilities
from astropy.table import Table


class Cataloguer:
    
    filesdir = None 
    image_names = None 
    set_size = None
    has_sets = None
    n_sets = None
    
    def __init__(self, dir, image_names, has_sets, set_size, n_sets):

        self.filesdir = dir
        self.image_names = image_names
        self.has_sets = has_sets
        self.n_sets = n_sets
        self.set_size = has_sets
        
    def catalogue(self):
                
        imagedir = self.filesdir + Constants.working_directory + Constants.image_directory 
            
        file = None
        
        if not self.has_sets:
            file = Constants.reduced_prefix + self.image_names + "0001" + Constants.fits_extension 
        else:
            file = Constants.reduced_prefix + self.image_names + "_1_001" + Constants.fits_extension
        
        if not self.has_sets:
            self.n_sets = 1
        
        for i in range(1, self.n_sets+1):
            
            image_data = fits.getdata(imagedir + file, ext=0)
            mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)    
            
            sources = self.find_stars(image_data, std)
            
            self.convert_to_ra_and_dec(sources, imagedir + file)
            if self.has_sets:
                filepath = self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_" + str(i) + Constants.standard_file_extension
            else:
                filepath = self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension

            sources.write(filepath, format = Constants.table_format, overwrite=True)
            self.make_reg_file(sources)
            
            if not self.has_sets:
                file = Constants.reduced_prefix + self.image_names + "000" + str(i) + Constants.fits_extension 
            else:
                file = Constants.reduced_prefix + self.image_names + "_" + str(i) + "_001" + Constants.fits_extension
        
        
        times = self.filesdir + "workspace/" + Constants.time_file
        if(os.path.exists(times)):
            open(times, "w").close()


    
        file = imagedir + file
        i = 2
        set = 1
        
        while (os.path.exists(file)):
            print(times)
            self.add_times(times, fits.getheader(file))
            print(i)
            if not self.has_sets:
                file = imagedir + Constants.reduced_prefix + self.image_names + "000" + str(i) + ".fits"
            
            else:
                
                if i > self.set_size:
                    i = 1
                    set+= 1
                
                file = imagedir + Constants.reduced_prefix + self.image_names + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension
            
            i+= 1
            

    
    def find_stars(self, image, std):
        daofind = DAOStarFinder(fwhm=8, threshold=3*std)   # find sources with FWHMs of ~8 pixels and peaks ~3-sigma above the background
        sources = daofind(image)

        for col in sources.colnames:    
            sources[col].info.format = '%.8g'  # for consistent table output
           
        return sources
   
    def convert_to_ra_and_dec(self, sources, image_file):
        
        # find the wcs assosiated with the fits image using astropy and the header
        wcs = WCS(fits.open(image_file)[0].header)

        # make two new coloums of 0's
        sources['RA'] = sources['xcentroid'] * 0
        sources['DEC'] = sources['xcentroid'] * 0

        # replace the 0's with ra and dec
        for x in range(0,len(sources)):
            ra, dec = wcs.all_pix2world(sources['xcentroid'][x], sources['ycentroid'][x], 0) 
            sources['RA'][x] = ra
            sources['DEC'][x] = dec
            
    def add_times(self, file, header):
                
        f = open(file, "a+")
        print(file)

        f.write(str(header['DATE-OBS']) + "\r\n")
        
    def make_reg_file(self, table):
        
        f = open(self.filesdir + self.image_names + ".reg", "a+")
        
        xs = table['xcentroid']
        ys = table['ycentroid']
        
        for i in range(len(xs)):
            f.write("point " + str(xs[i]) + " " + str(ys[i]) + " # point=circle 4 \r\n")

    def map_star_ids(self, image_size):  
        
        previous_cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_2" + Constants.standard_file_extension, format = Constants.table_format)

        for i in range(3, self.n_sets):
            cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_" + str(i) + Constants.standard_file_extension, format = Constants.table_format)
            shifts = self.find_shift_between_catalogues(previous_cat, cat, image_size)
            previous_cat = cat
            print(shifts)
            
            
    def find_shift_between_catalogues(self, catalogue1, catalogue2, image_size):
        
        x1s = catalogue1['xcentroid']
        y1s = catalogue1['ycentroid']
        
        x2s = catalogue2['xcentroid']
        y2s = catalogue2['ycentroid']
        
        fluxes = catalogue1["flux"]
        
        max = 0
        
        for i in range(len(fluxes)):
            if float(fluxes[i]) > fluxes[max] and x1s[i] > image_size - 200 and x1s[i] < image_size + 200 and y1s[i] > image_size - 200 and y1s[i] < image_size + 200:
                max = i
            
        x = x1s[max]
        y = y1s[max]
        
        print(max)
        distances = []
        
        for i in range(len(x1s)):
            
            if i != max:
                shifts = [0, 0]
            
                shifts[0] = x1s[i] - x
                shifts[1] = y1s[i] - y
            
                distances.append(shifts)
            
        for i in range(3176, len(x2s)):
            matches = 0
            matched = {}


            for j in range(len(x2s)):
                if i != j:
                    shift = [x2s[j] - x2s[i], y2s[j] - y2s[i]]
                    
                    for k in range(len(distances)):
                        if distances[k][0]*0.999 < shift[0] and distances[k][0]*1.001 > shift[0] and distances[k][1]*0.999 < shift[1] and distances[k][1]*1.001 > shift[1] and not k in matched:
                            matches += 1
                            matched[k] = True
                    
            print(matches, len(distances))
            if matches > len(distances) * 0.01:
                return x2s[i] - x1s[max], y2s[i] - y1s[max]
                    
                    
            
            
        
                
                    

                    
                    
                    
                    
                    
                      
                            
                            
                
                            
                    
                            
            
        
                    
                        
            
            
            
        
        
        

        
        
            
            
        
        
        

        
            
        
        
        
        
        
    
    