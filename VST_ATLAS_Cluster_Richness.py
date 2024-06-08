import os
import numpy as np
import pandas as pd
from scipy import stats
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from astropy.coordinates import SkyCoord

def calc_kcor(filter_name, redshift, colour_name, colour_value):
    """
    K-corrections calculator in Python.
    Reference: http://kcor.sai.msu.ru
    @param filter_name: Name of the filter to calculate K-correction for.
    @param redshift: Redshift of a galaxy.
    @param colour_name: Name of the colour, e.g., 'g - r'.
    @param colour_value: Value of the galaxy's colour.
    @return: K-correction in the specified filter for the given redshift and colour.
    """
    
    coeff = {
        'B_BRc': [[0, 0, 0, 0], [-1.99412, 3.45377, 0.818214, -0.630543], ...],
        # Add the remaining coefficients here...
    }

    c = coeff[filter_name + '_' + colour_name.replace(' - ', '')]
    kcor = 0.0

    for x, a in enumerate(c):
        for y, b in enumerate(c[x]):
            kcor += c[x][y] * redshift**x * colour_value**y
                
    return kcor

if __name__ == "__main__":
    import doctest
    doctest.testmod()

################################################################################################

# Constants
h = 0.678  # Planck 2018

# Load data
df = pd.read_csv('./catalogue/ATLAS_new_29_nulls_v2_no_assoc_p2_gsgt5lt300_nagtive_ra_merged_subsample.csv', sep=",")
df = df[df['clusterID_new'].groupby(df['clusterID_new']).transform('size') >= 5]

# 3 sigma clipping function to remove outliers
def is_outlier(s):
    lower_limit = s.mean() - (s.std() * 3)
    upper_limit = s.mean() + (s.std() * 3)
    return ~s.between(lower_limit, upper_limit)

df = df[~df.groupby('clusterID_new')['col4'].apply(is_outlier)]

# Calculations for div_1, div_4, div_7 and their errors
df['div_1'] = df.apply(lambda row: row.col1 / row.col3**2, axis=1)
df['div_1_err'] = df.apply(lambda row: 1 / row.col3**2, axis=1)
df['div_4'] = df.apply(lambda row: row.col4 / row.col5**2, axis=1)
df['div_4_err'] = df.apply(lambda row: 1 / row.col5**2, axis=1)
df['div_7'] = df.apply(lambda row: row.col7 / row.col8**2, axis=1)
df['div_7_err'] = df.apply(lambda row: 1 / row.col8**2, axis=1)

# Calculate g - i color
df['gi'] = df['gmag'] - df['imag']

# Group data by cluster ID
groupedp = df.groupby('clusterID_new')

# Calculate size, mean RA and Dec for each group
sizep = df.groupby('clusterID_new').size()  
ra_meanp = groupedp['ra'].mean()
dec_meanp = groupedp['dec'].mean()

# Calculate weighted means
wmp_1 = (groupedp.apply(lambda g: g.div_1.sum())) / groupedp.apply(lambda g: g.div_1_err.sum())    
wmp_4 = (groupedp.apply(lambda g: g.div_4.sum())) / groupedp.apply(lambda g: g.div_4_err.sum())    
wmp_7 = (groupedp.apply(lambda g: g.div_7.sum())) / groupedp.apply(lambda g: g.div_7_err.sum())    

# Calculate standard error of the mean
wm_errp_1 = groupedp['col1'].sem()
wm_errp_4 = groupedp['col4'].sem()
wm_errp_7 = groupedp['col7'].sem()

# Calculate angular diameter distance
D_A = cosmo.angular_diameter_distance(wmp_4).value

# Convert 1 h^-1 Mpc to degrees
size = (1000.0 / (h * cosmo.kpc_proper_per_arcmin(wmp_4).value)) / 60.0   
size_Abell = size.copy()
factor = cosmo.kpc_proper_per_arcmin(wmp_4).value
max_sep = size * u.deg

# Initialize arrays
n = np.zeros(groupedp.ngroups)
nr = np.zeros(groupedp.ngroups)
n200 = np.zeros(groupedp.ngroups)
B = np.zeros(groupedp.ngroups)
sep_arr = np.zeros(groupedp.ngroups)
sep_arr_deg = np.zeros(groupedp.ngroups)

# Calculate luminosity distance
d_L = cosmo.luminosity_distance(wmp_4).value

# Calculate L_s and M_s
L_s_i = 2.08 * (10**10) * (0.7**(-2))
M_s_i2 = -2.5 * np.log10((0.4 * L_s_i))

# Calculate K-correction
gi = groupedp.apply(lambda g: g.gi.mean())
K_cor = calc_kcor('i', wmp_4, 'g - i', gi)

# Calculate i-band magnitude limit and BCG magnitude
i_mag_lm = M_s_i2 + 25 + 5 * np.log10(d_L) + (K_cor.values)
M_bcg = (groupedp['imag'].max()) - 25 - (5 * np.log10(d_L)) - (K_cor.values)

# Calculate L_bcg
L_bcg = 10**(0.4 * (4.74 - M_bcg))

# Iterate over each group to calculate separations and counts
for i, (clusterID_new, grp) in enumerate(groupedp):
    ra_mean = grp['ra'].mean()
    dec_mean = grp['dec'].mean()
    c = SkyCoord(ra_mean, dec_mean, unit=u.deg, frame='fk5')
    A = grp['imag'].max()
    B[i] = i_mag_lm[i]    

    grp = grp.loc[(grp['imag'] < A) & (grp['imag'] > B[i])]
    catalog = SkyCoord(grp['ra'], grp['dec'], unit=u.deg, frame='fk5')
    
    try:
        sep_arr[i] = np.max(c.separation(catalog).arcminute) * factor[i] / 1000.0
        sep_arr_deg[i] = (sep_arr[i] * 180.0) / (np.pi * D_A[i])
    except ValueError:
        pass

    for j in range(len(grp.index)):
        if c.separation(catalog[j]) < max_sep[i]:   
            n[i] += 1

idx = np.where([wmp_4 > 0.1])[1]
n[idx] = n[idx] * (220 / ((571 * wmp_4.values[idx]**2) - 819 * wmp_4.values[idx] + 292.8))

# Calculate R200
R200 = 0.142 * (n**0.6) * 180.0 / (np.pi * D_A) * u.deg

# Calculate nr within R200
for i, (clusterID_new, grp) in enumerate(groupedp):
    ra_mean = grp['ra'].mean()
    dec_mean = grp['dec'].mean()
    c = SkyCoord(ra_mean, dec_mean, unit=u.deg, frame='fk5')
    A = grp['imag'].max()
    B[i] = i_mag_lm[i]    

    grp = grp.loc[(grp['imag'] < A) & (grp['imag'] > B[i])]
    catalog = SkyCoord(grp['ra'], grp['dec'], unit=u.deg, frame='fk5')

    for j in range(len(grp.index)):
        if c.separation(catalog[j]) < R200[i]:   
            nr[i] += 1

nr_ns = nr.copy()
nr[idx] = nr[idx] * (220 / ((571.7 * wmp_4.values[idx]**2) - 819 * wmp_4.values[idx] + 292.8))

# Save results to file
np.savetxt('ATLAS_new_29_nulls_v2_no_assoc_gsgt5lt300_nagtive_ra_richness_post_viva_3sig_v2_subsample.txt',
           np.c_[np.asarray(sizep.index, dtype=np.int64), ra_meanp, dec_meanp, sizep, wmp_1, wm_errp_1, wmp_4, wm_errp_4, wmp_7, wm_errp_7, n, nr, max_sep, R200, K_cor, sep_arr, sep_arr_deg, factor, size_Abell],
           header='clusterID RA DEC size wm1 wm1err wm4 wm4err wm7 wm7err n nr200 radius_1mpc radius K_cor sep_arr sep_arr_deg factor Abell_radius',
           fmt=['%i', '%f', '%f', '%i', '%f', '%f', '%f', '%f', '%f', '%f', '%i', '%i', '%f', '%f', '%f', '%f', '%f', '%f', '%f'])
