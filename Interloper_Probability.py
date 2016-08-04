from __future__ import division
import numpy as np
import matplotlib.pyplot as pp
from astropy.io import ascii
from scipy.spatial import cKDTree
from scipy.stats import beta


# how far to search for objects
_dist_bound = 1.0 / 0.06
# how many sample points to try
_n_samps = 500

def calculate_probability(folder_name):

    f125w_file = '/Users/Jenna/Documents/Matt Mechtley Quasar Research/all_quasar_images/{folder_name}/F125W/{folder_name}_F125W_phot.cat'.format(
        folder_name=folder_name)
    f160w_file = '/Users/Jenna/Documents/Matt Mechtley Quasar Research/all_quasar_images/{folder_name}/F160W/{folder_name}_F160W_phot.cat'.format(
        folder_name=folder_name)

    cat_f160w = ascii.read(f160w_file,format = 'sextractor')
    cat_f125w = ascii.read(f125w_file,format = 'sextractor')

    xy_image = np.column_stack([cat_f160w['X_IMAGE'], cat_f160w['Y_IMAGE']])
    kd = cKDTree(xy_image)

    extents = [xy_image.min(axis=0), xy_image.max(axis=0)]

    test_x_image = np.random.uniform(low=extents[0][0], high=extents[1][0], size=_n_samps)
    test_y_image = np.random.uniform(low=extents[0][1], high=extents[1][1], size=_n_samps)

    test_pts = np.column_stack([test_x_image, test_y_image])

    dists, inds = kd.query(test_pts, k=1, distance_upper_bound=_dist_bound)

    # make cut based on those with found neighbors
    found_cat_f160w = cat_f160w[inds[dists < np.inf]]
    found_cat_f125w = cat_f125w[inds[dists < np.inf]]
    successes = len(found_cat_f160w)
    failures = test_pts.shape[0] - successes

    # check the definition of beta distribution, the +1 is basically so the math works even in the case where successes or failures is 0
    found_dist = beta(successes + 1, failures + 1)
    # this gives us a 2-tuple of the 68% confidence interval, i.e. there is a 68% probability the value falls between (min,max)
    # 68% is the usually-quoted "1-sigma" statistical uncertainty
    found_frac = found_dist.interval(0.68)

    print 'Without Color Criteria'
    print 'Samples: {:d} Matched: {:d} Frac: {:0.3g}-{:0.3g}' .format(successes+failures, successes, found_frac[0], found_frac[1])

    # now make cuts based on photometry
    phot_mask = (found_cat_f160w['MAG_AUTO'] < 26.4)

    found_cat_f160w = found_cat_f160w[phot_mask]
    found_cat_f125w = found_cat_f125w[phot_mask]
    successes = len(found_cat_f160w)
    failures = test_pts.shape[0] - successes

    found_dist = beta(successes + 1, failures + 1)
    found_frac = found_dist.interval(0.68)

    print 'With Mag Criteria'
    print 'Sample: {:d} Matched: {:d} Frac: {:0.3g}-{:0.3g}' .format(successes+failures, successes, found_frac[0], found_frac[1])

    #now make color cuts
    jh_color = found_cat_f125w['MAG_AUTO'] - found_cat_f160w['MAG_AUTO']
    color_mask = jh_color < 0.5
    found_cat_f160w = found_cat_f160w[color_mask]
    found_cat_f125w = found_cat_f125w[color_mask]

    successes = len(found_cat_f160w)
    failures = test_pts.shape[0] - successes

    found_dist = beta(successes + 1, failures + 1)
    found_frac = found_dist.interval(0.68)

    print 'With Color Criteria'
    print 'Sample: {:d} Matched: {:d} Frac: {:0.3g}-{:0.3g}'.format(successes + failures, successes, found_frac[0],
                                                                    found_frac[1])

    pp.scatter(test_x_image, test_y_image, c='Orange')
    pp.scatter(found_cat_f160w['X_IMAGE'], found_cat_f160w['Y_IMAGE'], c='Purple')
    pp.xlabel('X_IMAGE')
    pp.ylabel('Y_IMAGE')
    pp.title(folder_name)
    pp.show()


folder_names = ['NDWFS-J142516.30+325409.0', 'CFHQS-J003311.40-012524.9', 'QSO-J0005-0006',
                'SDSS-J012958.51-003539.7', 'SDSS-J020332.39+001229.3', 'SDSS-J205406.42-000514.8']

for folder_name in folder_names:
    try:
        calculate_probability(folder_name)
    except IOError:
        print('Could not find catalogs for {folder}' .format(folder=folder_name))
