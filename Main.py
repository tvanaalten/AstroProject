import AstroProject.Reducer
import AstroProject.Cataloguer
import AstroProject.ShiftFinder
import AstroProject.FluxFinder

r = AstroProject.Reducer.Reducer("/Users/Thomas/Documents/Thomas_test/","No filter","l198", "bias_001.fit", "flat-1_001.fit")
r.reduce(True)
n = r.get_set_size()

c = AstroProject.Cataloguer.Cataloguer("/Users/Thomas/Documents/Thomas_test/", "l198", True, 50, 7)
c.catalogue()

sf = AstroProject.ShiftFinder.ShiftFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True)
sf.get_all_shifts(50, 7)

ff = AstroProject.FluxFinder.FluxFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)
ff.find_all_fluxes()
ff.make_light_curves()
#c.map_star_ids(1500)

    



