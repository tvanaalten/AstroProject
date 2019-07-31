from astropy.table import Table
from astropy.io import fits
from photutils import centroid_2dg
import matplotlib.pyplot as plt
import os
import Constants
import Utilities

class ShiftFinder:
    
    directory = None
    image_names = None
    has_sets = None
    
    def __init__(self, dir, names, has_sets):
        self.directory = dir
        self.image_names = names
        self.has_sets = has_sets
    
    def get_reference_coordinates(self, set):
        
        if not self.has_sets:
            t = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension, format=Constants.table_format)
        else:
            t = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_" + str(set) + Constants.standard_file_extension, format=Constants.table_format)

        xs = t["xcentroid"]
        ys = t["ycentroid"]
        
        fluxes = t["flux"]
        max = 0
        
        for i in range(len(fluxes)):
            if float(fluxes[i]) > fluxes[max]:
                max = i
            
        x = xs[max]
        y = ys[max]
        
        size = 10
        
        if not self.has_sets:
            image = fits.open(self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names + "0001" + Constants.fits_extension)
        else:
            image = fits.open(self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names + "_" + str(set)+ "_001" + Constants.fits_extension)
        
        data_array = image[0].data
        data_square = data_array[int(y)-size:int(y)+size, int(x)-size:int(x)+size] # note your x,y coords need to be an int
        plt.imshow(data_square, origin='lower',  cmap='viridis')
        plt.plot(10, 10, 'r*', label='new center')
        plt.show()



                
        return x, y
    
    def find_shift(self, previous_x, previous_y, image_path):

        # choose how big a box to draw around your star
        size = 10  # 10 will make a 20x20 pixel box
        
        # open the second image and load it as an array
        image = fits.open(image_path)  # replace 'FullComb.fits' with the namne of your second image file
        data_array = image[0].data
        # Select just the box around the star
        data_square = data_array[int(previous_y)-size:int(previous_y)+size, int(previous_x)-size:int(previous_x)+size] # note your x,y coords need to be an int
        plt.imshow(data_square, origin='lower',  cmap='viridis')

        print(previous_x, previous_y)

        
        x, y = centroid_2dg(data_square)
        
        plt.plot(x, y, 'r*', label='new center')
        plt.show()


                
        x_shift = ((x - size) - (previous_x - int(previous_x))) # also correcting for the shift due to using int
        y_shift = ((y - size) - (previous_y - int(previous_y)))
        
        return x_shift, y_shift
    
        
    def get_all_shifts(self, set_size, n_sets):
                
        shift_file = self.directory + Constants.working_directory + Constants.shift_file
        
        if(os.path.exists(shift_file)):
            open(shift_file, "w").close()
        
        f = open(shift_file, "a+")

        x_shifts = []
        y_shifts = []
        
        image = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names 
    
        if not self.has_sets:
            image += "0002"
        else:
            image += "_1_002"
        
        image += Constants.fits_extension

        i = 3
        set = 1
        
        x, y = self.get_reference_coordinates(set)
                
        while (os.path.exists(image)):
                        
            x_shift, y_shift = self.find_shift(x, y, image)
            
            x_shifts.append(x_shift)
            y_shifts.append(y_shift)
            
            x += x_shift
            y += y_shift
            
            image = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names
            
            if not self.has_sets:            
                image += "000" + str(i)
            else:
                image += "_" + str(set) + "_" + Utilities.format_index(i)
            
            image += Constants.fits_extension
                
            print(set, i)

            i+= 1
            
            
            if self.has_sets and i > set_size:
                i = 2
                set += 1
                              
                if set > n_sets:
                    break
                
                
                x, y = self.get_reference_coordinates(set)

            
            
        table = Table([x_shifts, y_shifts], names = ('xshifts','yshifts'))
        
        table.write(shift_file, format = Constants.table_format, overwrite=True)

        
        


        
        
    
    