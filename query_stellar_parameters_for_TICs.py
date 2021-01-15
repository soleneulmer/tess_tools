from astropy import units as u
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
import pandas as pd

def get_gaia_data(tic):
    """
    Input: tic = float
    Ouputs: params = list containing Right Ascension, Declination, 
    Proper motions, and radial velocity of the star
    """
    try:
        # Query GAIA catalog
        r = Gaia.query_object('TIC'+str(tic), radius=5*u.arcsec)

        # Reading parameters of interest
        rv = round(r['radial_velocity'].data.data[0], 2)
        pmra = round(r['pmra'].data.data[0] * 10**(-3), 4)
        pmdec = round(r['pmdec'].data.data[0] * 10**(-3), 4)
        ra = round(r['ra'].data.data[0], 6)
        dec = round(r['dec'].data.data[0], 6)
        params = [ra, dec, pmra, pmdec, rv]
    except:
        params = [nan, nan, nan, nan, nan]
    return params

def get_exofop_data(tic):
    """
    Inputs: tic = float
    Ouputs: top_params = list containing Right Ascension, Declination, 
    Proper motions, Distance, V magnitude and B-V magnitude difference,
    and the stellar effective temperature
    btm_params = list containing stellar radius, mass, density, metallicity and logg
    """
    try:
        # Query TICv8 catalog on MAST
        r = Catalogs.query_object("TIC"+ str(tic), radius=.01, catalog="TIC")

        # Reading parameters of interest
        ra = round(r['ra'].data.data[0], 6)        # deg
        dec = round(r['dec'].data.data[0], 6)      # deg
        pmra = round(r['pmRA'].data.data[0], 6)    # mas/yr
        pmdec = round(r['pmDEC'].data.data[0], 6)  # mas/yr
        dist = round(r['d'].data.data[0], 4)       # pc

        Vmag = round(r['Vmag'].data.data[0], 4)
        Bmag = round(r['Bmag'].data.data[0], 4)
        BV = Bmag - Vmag

        Teff = round(r['Teff'].data.data[0], 0)   # Kelvins

        top_params = [ra, dec, pmra, pmdec, dist, Vmag, BV, Teff]

        radius = round(r['rad'].data.data[0], 4)          # R_sun
        mass = round(r['mass'].data.data[0], 4)           # M_sun
        density = round(r['rho'].data.data[0], 4) * 1.41  # g/cm3
        metallicity = round(r['MH'].data.data[0], 4)
        logg = round(r['logg'].data.data[0], 4)           # cm/s2


        btm_params = [radius, mass, density, metallicity, logg]
    except:
        top_params = [nan, nan, nan, nan, nan, nan, nan, nan]
        btm_params = [nan, nan, nan, nan, nan]
    return top_params, btm_params


def quadratic_ld_coefficients(teff, logg, path_to_file='/home/solene/transits/limb_darkening_Claret/table5_Claret2018.dat'):
    """
    Quadratic Limb darkening coefficients from Claret 2018
    Vizier repository: http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A+A/618/A20/
    Inputs: teff = float, stellar effective temperature
    logg = float, stellar logg
    path_to_file = string, path to the table 5 containing quadratic limb darkening coefficients for TESS
    """
    # Rounding to match the values in the table
    teff_ld = round(teff, -2)
    logg_ld = round(logg * 2) / 2

    # Loading table
    ld_df = pd.read_csv(path_to_file, delimiter='|', skiprows=5, skipfooter=1,
                        names=['logg', 'Teff', 'Z', 'LHP', 'a', 'b',
                           'mu', 'chi2', 'od', 'Sys'])
    # Selecting both LD coefficients based on Teff and logg
    a_ld = ld_df[(ld_df['logg'] == logg_ld) & (ld_df['Teff'] == teff_ld)]['a'].values[0]
    b_ld = ld_df[(ld_df['logg'] == logg_ld) & (ld_df['Teff'] == teff_ld)]['b'].values[0]
    return [a_ld, b_ld]
