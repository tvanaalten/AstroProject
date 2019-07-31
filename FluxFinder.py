import astropy
import os
from astropy.table import Table
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Constants
import Utilities

class FluxFinder:
    

    directory = None
    flux_dir = None
    light_curve_dir = None
    image_names = None
    x_shifts = None
    y_shifts = None
    catalogues = []
    set_size = None
    n_sets = None
    has_sets = None
    obs_start_time = []
    
    def __init__(self, directory, image_names, has_sets, n_sets, set_size):
        self.directory = directory
        self.image_names = image_names
        self.n_sets = n_sets
        self.set_size = set_size
        self.has_sets = has_sets
        if not has_sets:
            self.catalogues.append(Table.read(directory + Constants.working_directory + Constants.catalogue_prefix + image_names + Constants.standard_file_extension, format = Constants.table_format))
        else:
            for i in range(1, n_sets+1):
                self.catalogues.append(Table.read(directory + Constants.working_directory + Constants.catalogue_prefix + image_names + "_" + str(i) + Constants.standard_file_extension, format = Constants.table_format))

        self.get_shifts()
        
        self.flux_dir = directory + Constants.working_directory + Constants.flux_directory
        if not os.path.exists(self.flux_dir):
            os.mkdir(self.flux_dir)
        
        self.light_curve_dir = directory + Constants.working_directory  + Constants.light_curve_directory
        if not os.path.exists(self.light_curve_dir):
            os.mkdir(self.light_curve_dir)
        
    def find_all_fluxes(self):
        
        file = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names 
        
        if not self.has_sets:
            file += "_0001" 
        else:
            file +="_1_001"
        
        file += Constants.fits_extension

        i = 1
        set = 1
        print(file)
        while (os.path.exists(file)):
            self.find_fluxes(i, set)
            
            i+=1
            
            if self.has_sets and i > self.set_size:
                set += 1
                i = 1
           
            file = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names 
            
            if not self.has_sets:
                file += "000" + str(i) 
            else:
                file += "_" + str(set) + "_" + Utilities.format_index(i)

            file += Constants.fits_extension
            
            
        
    def get_shifts(self):
        
        shifts_path = self.directory + Constants.working_directory + Constants.shift_file
        t = Table.read(shifts_path, format = Constants.table_format)
        self.x_shifts = t['xshifts']
        self.y_shifts = t['yshifts']
        
    
    def get_total_shift(self, image_number, set):
    #get total shift for ith image, starting from 1 (shift at 1 = 0)
    
        total_x = 0
        total_y = 0
        for j in range((set-1)*(self.set_size-1)+1, (set-1)*(self.set_size-1) + image_number-1):
            total_x += self.x_shifts[j-1]
            total_y += self.y_shifts[j-1]
        
        print(total_x, total_y)
        return total_x, total_y
    
    def find_fluxes(self, image_number, set):

        x_shift, y_shift = self.get_total_shift(image_number, set)
        positions = (self.catalogues[set-1]['xcentroid'] + x_shift, self.catalogues[set-1]['ycentroid'] + y_shift) # x,y positions of the aperture centers
        apertures = CircularAperture(positions, r=9) # Shape and size of aperture object
        
        data_path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names
        if not self.has_sets:
            data_path += "000" + str(image_number)
        else:
            data_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
        
        data_path += Constants.fits_extension
            
        # find the counts in each aperture

        image_data = fits.getdata(data_path, ext = 0)
        phot_table = aperture_photometry(image_data, apertures) 
        
        # Local background subtraction

        # Define size of background aperture
        annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)
        apertures = CircularAperture(positions, r=9)
        apers = [apertures, annulus_apertures]

        # find counts in each aperture and annulus
        phot_table2 = aperture_photometry(image_data, apers)

        for col in phot_table2.colnames:
            phot_table2[col].info.format = '%.8g'  # for consistent table output
    
        # MEAN
        # Background sub using mean
        
        # Calc the mean background in the second aperture ring
        bkg_mean = phot_table2['aperture_sum_1'] / annulus_apertures.area()
        phot_table2['mean'] = bkg_mean
        
        # Calc background level in each main aperture and subtract
        bkg_sum = bkg_mean * apertures.area()
        final_sum = phot_table2['aperture_sum_0'] - bkg_sum
        phot_table2['residual_aperture_sum_mean'] = final_sum
        phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
                
        # MEDIAN

        # Background sub using median
        phot_table2['median'] = phot_table2['aperture_sum_1']*0
        
        for q in range(0,len(phot_table2)):
            x = positions[0][q]
            y = positions[1][q]   
            xypos = (x, y)
            annulus = CircularAnnulus(xypos, r_in=10, r_out=15)
            ann_mask = annulus.to_mask(method='center')[0]
            weighted_data = ann_mask.multiply(image_data)
            phot_table2['median'][q] = np.median(weighted_data[weighted_data != 0])
            
        # Calc the median background in the second aperture ring
        bkg_med = phot_table2['median']
        
        # Calc background level in each main aperture and subtract
        bkg_sum = bkg_med * apertures.area()
        final_sum = phot_table2['aperture_sum_0'] - bkg_sum
        
        phot_table2['residual_aperture_sum_med'] = final_sum
        phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
        
        out_path = self.flux_dir + Constants.flux_prefix + self.image_names
        
        if not self.has_sets:
            out_path += "000" + str(image_number)
        else:
            out_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
            
        out_path += Constants.standard_file_extension
        phot_table2.write(out_path,  format = Constants.table_format, overwrite = True)
        
    def make_light_curves(self):
        
        times_file = self.directory + Constants.working_directory + Constants.time_file
        
        times = [line.rstrip('\n') for line in open(times_file)]

        file = self.flux_dir + Constants.flux_prefix + self.image_names
        if not self.has_sets:
            file += "0001" 
        else:
            file += "_1_001"
        file += Constants.standard_file_extension
        
        i = 1
        set = 1
        
        all_light_curves = []
        light_curves = []
        while (os.path.exists(file)):
            t = Table.read(file, format = Constants.table_format)

            i+= 1
            
                
            if self.has_sets and self.set_size < i:
                i = 1
                set+= 1
                all_light_curves.append(light_curves)
                light_curves = []
                break
            
            if set > self.n_sets:
                break
           
            file = self.flux_dir + Constants.flux_prefix + self.image_names
            if not self.has_sets:
                file += "000" + str(i) 
            else:
                file += "_" + str(set) + "_" + Utilities.format_index(i)
            file += Constants.standard_file_extension
            
            print(file)
            
            for j in range(len(t['id'])):
                if (i == 2 and set == 1) or i == 1:
                    light_curves.append(Table(names = ('time','counts')))
                
                date_and_time = times[i-1]
                time = date_and_time.split("T")[1]
                
                
                time_split = time.split(":")
                time_split[0] = int(time_split[0])
                time_split[1] = int(time_split[1])
                time_split[2] = int(time_split[2])
                
                if len(self.obs_start_time) == 0:
                    self.obs_start_time = time_split
                    
                time_elapsed_hms = Utilities.diff(time_split, self.obs_start_time)
                
                
                if time_elapsed_hms[0] >= 0:
                
                    time_elapsed = time_elapsed_hms[0] * 3600 + time_elapsed_hms[1] * 60 + time_elapsed_hms[2]
                
                else:
                    
                    time_elapsed_hms_till_midnight = Utilities.diff([23, 59, 60], self.obs_start_time)
                    time_elapsed = time_elapsed_hms_till_midnight[0] * 3600 + time_elapsed_hms_till_midnight[1] * 60 + time_elapsed_hms_till_midnight[2]
                    time_elapsed += time_split[0] * 3600 + time_split[1] * 60 + time_split[0]
                    
                
                                
                
                light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_mean'][j])))

            
            

        
        i = 0
        #for i in range(len(all_light_curves)):
        for j in range(len(all_light_curves[i])):
            file = self.light_curve_dir + self.image_names + Constants.identifier + str(j+1) + Constants.standard_file_extension
           
            table = all_light_curves[i][j]
            table.write(file, format = Constants.table_format, overwrite=True)
        
    def plot_light_curve(self, id):
        file = self.light_curve_dir + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension
        table = Table.read(file, format=Constants.table_format)
        times = table['time']
        fluxes = table['counts']
        
        plt.plot(times, fluxes)
        plt.xlabel("time (s)")
        plt.ylabel("flux")
        
        
        
            
            
        
            
    
    
    
    
            