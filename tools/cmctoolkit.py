import numpy as np
import pandas as pd
import h5py
import scipy.optimize
import scipy.interpolate
import gzip
import time
import os

################################################################################
# DEFINE CGS CONSTANTS
################################################################################

#Universal constants
h, c, k, hbar = 6.6260755e-27, 2.99792458e10, 1.380658e-16, 1.05457266e-27
G, sigma = 6.67259e-8, 5.67051e-5

# Subsets
startype_all = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
startype_star = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
startype_ms = np.array([0, 1])
startype_giant = np.array([2, 3, 4, 5, 6, 7, 8, 9])
startype_wd = np.array([10, 11, 12])
startype_other = np.array([7])
startype_remnant = np.array([10, 11, 12, 13, 14])
startype_bh = np.array([14])

def make_unitdict(convfile):
    """
    Helper function which converts a conversion ratio file into a unit dictionary.
    An example of a conversion file is 'initial.conv.sh'. Function is called in
    Snapshot class.

    Parameters
    ----------
    convfile: list
        specially formatted convfile list wherein each line is an element
    """
    base_dict = {'code': 1, # code units

                 # Fundamental
                    'g': float(convfile[5][12:]), # grams, massunitcgs
                 'msun': float(convfile[7][13:]), # msun, massunitmsun
                   'cm': float(convfile[13][14:]), # cm, lengthunitcgs
                   'pc': float(convfile[15][17:]), # pc, lengthunitparsec
                    's': float(convfile[17][12:]), # s, timeunitcgs
                  'myr': float(convfile[19][13:]), # myr, timeunitsmyr

                 # Stellar
                  's_g': float(convfile[9][13:]), # g stellar, mstarunitcgs
               's_msun': float(convfile[11][14:]), # msun stellar, mstarunitmsun

                 # N-body
                 'nb_s': float(convfile[21][14:]), # s N-body s, nbtimeunitcgs
               'nb_myr': float(convfile[23][15:]), # myr N-body Myr, nbtimeunitsmyr
                }

    custom_dict = {
                 # Custom
                 'nb_km/s': 1e-5 * base_dict['cm'] / base_dict['nb_s'], # km/s
                     'erg': base_dict['g'] * base_dict['cm'] ** 2 / base_dict['s'] ** 2, # erg
                   'erg.s': base_dict['g'] * base_dict['cm'] ** 2 / base_dict['s'], # erg*s (angular momentum)
                   'erg/s': base_dict['g'] * base_dict['cm'] ** 2 / base_dict['s'] ** 3, # erg/s
      'erg/s/cm2/angstrom': base_dict['g'] / (1e8 * base_dict['cm'] * base_dict['s'] ** 3), # erg/s/cm2/angstrom
               'erg/s/cm3': base_dict['g'] / (base_dict['cm'] * base_dict['s'] ** 3), # erg/s/cm3
                    'lsun': 2.599e-34 * base_dict['g'] * base_dict['cm'] ** 2 / base_dict['s'] ** 3, # lsun
                    'rsun': (1 / 6.96e10) * base_dict['cm'], # rsun
                'angstrom': 1e8 * base_dict['cm'], # angstrom
                     'kpc': 1e-3 * base_dict['pc'], # kpc
                    'g/s2': base_dict['g'] / base_dict['s'] ** 2, # g/s^2 (spectral flux unit)
                      'jy': 1e23 * base_dict['g'] / base_dict['s'] ** 2, # jansky
                     'gyr': 1e-3 * base_dict['myr'], # gyr
                      'kg': 1e-3 * base_dict['g'], #kg
                      'Hz': base_dict['s'] ** -1, # Hz (frequency)
                    '1/yr': 31556926.0153 * base_dict['s'] ** -1, # yr^-1 (frequency)
                   'nb_Hz': base_dict['nb_s'] ** -1, # n-body Hz (frequency)

                 # Angular
                     'rad': 1, # radians
                     'deg': 180 / np.pi, # degrees
                  'arcmin': 60 * 180 / np.pi, # arcmin
                  'arcsec': 3600 * 180 / np.pi, # arcsec
                     'mas': 1e3 * 3600 * 180 / np.pi, # milliarcsecond
                  }

    unitdict = {**base_dict, **custom_dict}

    return unitdict

def load_filter(fname):
    """
    Filter function which loads ascii filter functions with two columns:
       1. wavelength[ANGSTROM]
       2. transmission
    Function assumes no header. Also assumes that the last line is empty.

    Parameters
    ----------
    fname: str
        path of filter function

    Returns
    -------
    filt: pd.DataFrame
        filter function table
    """
    # Read filter function
    f = open(fname, 'r')
    text = f.read().split('\n')[:-1] # Assumes last line is empty
    f.close()

    # Convert to a pandas table
    wavelength = np.array([text[ii].split()[0] for ii in range(len(text))])
    transmission = np.array([text[ii].split()[1] for ii in range(len(text))])

    filt = {'wavelength[ANGSTROM]': wavelength.astype(float),
            'transmission': transmission.astype(float)}
    filt = pd.DataFrame(filt)

    return filt

def load_filtertable(fname):
    """
    Filter function which loads ascii filter functions with five columns:
       1. filtname
       2. path
       3. zp_spectralflux[ERG/S/CM2/ANGSTROM]
    Function assumes no header. Also assumes that the last line is empty.

    Parameters
    ----------
    fname: str
        path of filter table file

    Returns
    -------
    filttable: pd.DataFrame
        table containing information about filter functions
    """
    # Read filter table
    f = open(fname, 'r')
    text = f.read().split('\n')
    f.close()

    # Convert to a pandas table
    filtname = np.array([text[ii].split()[0] for ii in range(len(text))])
    path = np.array([text[ii].split()[1] for ii in range(len(text))])
    zp_spectralflux = np.array([text[ii].split()[2] for ii in range(len(text))])
    
    filttable = {'filtname': filtname,
            'path': path,
            'zp_spectralflux[ERG/S/CM2/ANGSTROM]': zp_spectralflux.astype(float)}
    filttable = pd.DataFrame(filttable)

    return filttable

def smooth_filter(filtfunc, seglength=100):
    """
    Takes a filter function with many points and smooth it over
    """
    filttable = load_filter(filtfunc)
    
    # Open new file
    f = open(filtfunc.replace('.dat', '') + '_abridged.dat', 'w')
    
    # Partition the filter function into seglength point segments (except for the last
    # one)
    for ii in range(int(np.ceil(len(filttable) / seglength))):
        segment = filttable[seglength * ii: seglength * (ii + 1)]
        
        # Calculate area under the curve of the filter function
        wavelength_arr = np.array(segment['wavelength[ANGSTROM]'])
        transmission_arr = np.array(segment['transmission'])
        area = np.sum(0.5 * (transmission_arr[1:] + transmission_arr[:-1])
                          * (wavelength_arr[1:] - wavelength_arr[:-1]))
                                  
        wavelength = (wavelength_arr[-1] + wavelength_arr[0]) / 2
        trans = area / (wavelength_arr[-1] - wavelength_arr[0])
                
        f.write('{} {}\n'.format(wavelength, trans))
        
    f.close()

def add_mags(mag1, mag2):
    """
    Helper function which adds two magnitudes together along luminosity lines.

    Parameters
    ----------
    mag1: float
        first magnitude being added

    mag2: float
        second magnitude being added

    Returns
    -------
    tot_mag: float
        sum of the magnitudes
    """
    tot_mag = -2.5 * np.log10( 10 ** (-mag1 / 2.5) + 10 ** (-mag2 / 2.5) )
    return tot_mag

def find_t_ms(z, m):
    """
    Helper function for find_MS_TO()
    """
    eta = np.log10(z/0.02)
    a1 = 1.593890e3+2.053038e3*eta+1.231226e3*eta**2.+2.327785e2*eta**3.
    a2 = 2.706708e3+ 1.483131e3*eta+ 5.772723e2*eta**2.+ 7.411230e1*eta**3.
    a3 = 1.466143e2 - 1.048442e2*eta - 6.795374e1*eta**2. - 1.391127e1*eta**3.
    a4 = 4.141960e-2 + 4.564888e-2*eta + 2.958542e-2*eta**2 + 5.571483e-3*eta**3.
    a5 = 3.426349e-1
    a6 = 1.949814e1 + 1.758178*eta - 6.008212*eta**2. - 4.470533*eta**3.
    a7 = 4.903830
    a8 = 5.212154e-2 + 3.166411e-2*eta - 2.750074e-3*eta**2. - 2.271549e-3*eta**3.
    a9 = 1.312179 - 3.294936e-1*eta + 9.231860e-2*eta**2. + 2.610989e-2*eta**3.
    a10 = 8.073972e-1

    m_hook = 1.0185 + 0.16015*eta + 0.0892*eta**2.
    m_HeF = 1.995 + 0.25*eta + 0.087*eta**2.
    m_FGB = 13.048*(z/0.02)**0.06/(1+0.0012*(0.02/z)**1.27)

    t_BGB = (a1+a2*m**4.+a3*m**5.5+m**7.)/(a4*m**2.+a5*m**7.)
    x = np.max([0.95,np.min([0.95-0.03*(eta+0.30103)]),0.99])
    mu = np.max([0.5, 1.0-0.01*np.max([a6/(m**a7), a8+a9/m**a10])])
    t_hook = mu*t_BGB

    t_MS = np.max([t_hook, x*t_BGB])

    return (t_MS)

def SSE_MS_get_L_and_R(M, Z, t):
    """
    Use SSE (Hurley et al. 2000) to get L and R as a function of M, Z, t. M should
    be given as an array and Z, t should be single numbers.
    Main sequence implemented here only.
    
    M: array of masses given in solar masses
    Z: float absolute metallicity (Zsun = 0.02)
    t: age (Gyr)
    """
    # Helper functions  
    def get_zdep_param(Z, alpha, beta, gamma, eta, mu):
        """
        Helper function to retrieve a parameter at a certain SSE/BSE metallicity
        parametezed as in Hurley et al. 2000 and Tout et al. 1996 as:

        param = alpha + beta * k + gamma * k^2 + eta * k^3 + mu * k^4

        where

        zeta = log(Z / 0.02)
        """
        zeta = np.log10(Z / 0.02)

        param = alpha + beta * zeta + gamma * zeta ** 2 + eta * zeta ** 3 + mu * zeta ** 4

        return param
    
    zeta = np.log10(Z / 0.02) # called zeta in text
    t *= 1000 # convert Gyr to Myr (SSE's units)
    
    # Retrieve ZAMS L and R from Tout et al. 1996
    alpha_tout = get_zdep_param(Z, 0.39704170, -0.32913574, 0.34776688, 0.37470851, 0.09011915)
    beta_tout = get_zdep_param(Z, 8.52762600, -24.41225973, 56.43597107, 37.06152575, 5.45624060)
    gamma_tout = get_zdep_param(Z, 0.00025546, -0.00123461, -0.00023246, 0.00045519, 0.00016176)
    delta_tout = get_zdep_param(Z, 5.43288900, -8.62157806, 13.44202049, 14.51584135, 3.39793084)
    epsilon_tout = get_zdep_param(Z, 5.56357900, -10.32345224, 19.44322980, 18.97361347, 4.16903097)
    zeta_tout = get_zdep_param(Z, 0.78866060, -2.90870942, 6.54713531, 4.05606657, 0.53287322)
    eta_tout = get_zdep_param(Z, 0.00586685, -0.01704237, 0.03872348, 0.02570041, 0.00383376)

    theta_tout = get_zdep_param(Z, 1.71535900, 0.62246212, -0.92557761, -1.16996966, -0.30631491)
    iota_tout = get_zdep_param(Z, 6.59778800, -0.42450044, -12.13339427, -10.73509484, -2.51487077)
    kappa_tout = get_zdep_param(Z, 10.08855000, -7.11727086, -31.67119479, -24.24848322, -5.33608972)
    lambda_tout = get_zdep_param(Z, 1.01249500, 0.32699690, -0.00923418, -0.03876858, -0.00412750)
    mu_tout = get_zdep_param(Z, 0.07490166, 0.02410413, 0.07233664, 0.03040467, 0.00197741)
    nu_tout = get_zdep_param(Z, 0.01077422, 0, 0, 0, 0)
    xi_tout = get_zdep_param(Z, 3.08223400, 0.94472050, -2.15200882, -2.49219496, -0.63848738)
    omicron_tout = get_zdep_param(Z, 17.84778000, -7.45345690, -48.96066856, -40.05386135, -9.09331816)
    pi_tout = get_zdep_param(Z, 0.00022582, -0.00186899, 0.00388783, 0.00142402, -0.00007671)

    Lzams = (alpha_tout * M ** 5.5 + beta_tout * M ** 11) / (gamma_tout + M ** 3 + delta_tout * M ** 5 + epsilon_tout * M ** 7 + zeta_tout * M ** 8 + eta_tout * M ** 9.5)
    Rzams = (theta_tout * M ** 2.5 + iota_tout * M ** 6.5 + kappa_tout * M ** 11 + lambda_tout * M ** 19 + mu_tout * M ** 19.5) / (nu_tout + xi_tout * M ** 2 + omicron_tout * M ** 8.5 + M ** 18.5 + pi_tout * M ** 19.5)

    # Model parameters
    a1 = get_zdep_param(Z, 1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2, 0)
    a2 = get_zdep_param(Z, 2.706708e3, 1.483131e3, 5.772723e2, 7.411230e1, 0)
    a3 = get_zdep_param(Z, 1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1, 0)
    a4 = get_zdep_param(Z, 4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3, 0)
    a5 = get_zdep_param(Z, 3.426349e-1, 0, 0, 0, 0)
    a6 = get_zdep_param(Z, 1.949814e1, 1.758178e0, -6.008212e0, -4.470533e0, 0)
    a7 = get_zdep_param(Z, 4.903830e0, 0, 0, 0, 0)
    a8 = get_zdep_param(Z, 5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3, 0)
    a9 = get_zdep_param(Z, 1.312179e0, -3.294936e-1, 9.231860e-2, 2.610989e-2, 0)
    a10 = get_zdep_param(Z, 8.073972e-1, 0, 0, 0, 0)
    a11p = get_zdep_param(Z, 1.031538e0, -2.434480e-1, 7.732821e0, 6.460705e0, 1.374484e0)
    a12p = get_zdep_param(Z, 1.043715e0, -1.577474e0, -5.168234e0, -5.596506e0, -1.299394e0)
    a13 = get_zdep_param(Z, 7.859573e2, -8.542048e0, -2.642511e1, -9.585707e0, 0)
    a14 = get_zdep_param(Z, 3.858911e3, 2.459681e3, -7.630093e1, -3.486057e2, -4.861703e1)
    a15 = get_zdep_param(Z, 2.888720e2, 2.952979e2, 1.850341e2, 3.797254e1, 0)
    a16 = get_zdep_param(Z, 7.196580e0, 5.613746e-1, 3.805871e-1, 8.398728e-2, 0)

    a11 = a11p * a14
    a12 = a12p * a14

    sigma = np.log10(Z)
    log_a17 = np.max([0.097 - 0.1072 * (sigma + 3), np.max([0.097, np.min([0.1461, 0.1461 + 0.1237 * (sigma + 2)])])])
    a17 = 10 ** log_a17

    a18p = get_zdep_param(Z, 2.187715e-1, -2.154437e0, -3.768678e0, -1.975518e0, -3.021475e-1)
    a19p = get_zdep_param(Z, 1.466440e0, 1.839725e0, 6.442199e0, 4.023635e0, 6.957529e-1)
    a20 = get_zdep_param(Z, 2.652091e1, 8.178458e1, 1.156058e2, 7.633811e1, 1.950698e1)
    a21 = get_zdep_param(Z, 1.472103e0, -2.947609e0, -3.312828e0, -9.945065e-1, 0)
    a22 = get_zdep_param(Z, 3.071048e0, -5.679941e0, -9.745523e0, -3.594543e0, 0)
    a23 = get_zdep_param(Z, 2.617890e0, 1.019135e0, -3.292551e-2, -7.445123e-2, 0)
    a24 = get_zdep_param(Z, 1.075567e-2, 1.773287e-2, 9.610479e-3, 1.732469e-3, 0)
    a25 = get_zdep_param(Z, 1.476246e0, 1.899331e0, 1.195010e0, 3.035051e-1, 0)
    a26 = get_zdep_param(Z, 5.502535e0, -6.601663e-2, 9.968707e-2, 3.599801e-2, 0)

    a18 = a18p * a20
    a19 = a19p * a20
    
    a27 = get_zdep_param(Z, 9.511033e1, 6.819618e1, -1.045625e1, -1.474939e1, 0)
    a28 = get_zdep_param(Z, 3.113458e1, 1.012033e1, -4.650511e0, -2.463185e0, 0)
    a29p = get_zdep_param(Z, 1.413057e0, 4.578814e-1, -6.850581e-2, -5.588658e-2, 0)
    a30 = get_zdep_param(Z, 3.910862e1, 5.196646e1, 2.264970e1, 2.873680e0, 0)
    a31 = get_zdep_param(Z, 4.597479e0, -2.855179e-1, 2.709724e-1, 0, 0)
    a32 = get_zdep_param(Z, 6.682518e0, 2.827718e-1, -7.294429e-2, 0, 0)
    
    a29 = a29p ** a32
    
    a33 = np.min([1.4, 1.5135 + 0.3769 * zeta]) # sic
    a33 = np.max([0.6355 - 0.4192 * zeta, np.max([1.25, a33])])
    
    a34 = get_zdep_param(Z, 1.910302e-1, 1.158624e-1, 3.348990e-2, 2.599706e-3, 0)
    a35 = get_zdep_param(Z, 3.931056e-1, 7.277637e-2, -1.366593e-1, -4.508946e-2, 0)
    a36 = get_zdep_param(Z, 3.267776e-1, 1.204424e-1, 9.988332e-2, 2.455361e-2, 0)
    a37 = get_zdep_param(Z, 5.990212e-1, 5.570264e-2, 6.207626e-2, 1.777283e-2, 0)
    
    a38 = get_zdep_param(Z, 7.330122e-1, 5.192827e-1, 2.316416e-1, 8.346941e-3, 0)
    a39 = get_zdep_param(Z, 1.172768e0, -1.209262e-1, -1.193023e-1, -2.859837e-2, 0)
    a40 = get_zdep_param(Z, 3.982622e-1, -2.296279e-1, -2.262539e-1, -5.219837e-2, 0)
    a41 = get_zdep_param(Z, 3.571038e0, -2.223625e-2, -2.611794e-2, -6.359648e-3, 0)
    a42 = get_zdep_param(Z, 1.9848e0, 1.1386e0, 3.5640e-1, 0, 0)
    a43 = get_zdep_param(Z, 6.300e-2, 4.810e-2, 9.840e-3, 0, 0)
    a44 = get_zdep_param(Z, 1.200e0, 2.450e0, 0, 0, 0)
    
    a42 = np.min([1.25, np.max([1.1, a42])]) # sic
    a44 = np.min([1.3, np.max([0.45, a44])])
    
    a45 = get_zdep_param(Z, 2.321400e-1, 1.828075e-3, -2.232007e-2, -3.378734e-3, 0)
    a46 = get_zdep_param(Z, 1.163659e-2, 3.427682e-3, 1.421393e-3, -3.710666e-3, 0)
    a47 = get_zdep_param(Z, 1.048020e-2, -1.231921e-2, -1.686860e-2, -4.234354e-3, 0)
    a48 = get_zdep_param(Z, 1.555590e0, -3.223927e-1, -5.197429e-1, -1.066441e-1, 0)
    a49 = get_zdep_param(Z, 9.7700e-2, -2.3100e-1, -7.5300e-2, 0, 0)
    a50 = get_zdep_param(Z, 2.4000e-1, 1.8000e-1, 5.9500e-1, 0, 0)
    a51 = get_zdep_param(Z, 3.3000e-1, 1.3200e-1, 2.1800e-1, 0, 0)
    a52 = get_zdep_param(Z, 1.1064e0, 4.1500e-1, 1.8000e-1, 0, 0)
    a53 = get_zdep_param(Z, 1.1900e0, 3.7700e-1, 1.7600e-1, 0, 0)
    
    a49 = np.max([a49, 0.145]) # sic
    a50 = np.min([a50, 0.306 + 0.053 * zeta])
    a51 = np.min([a51, 0.3625 + 0.062 * zeta])
    a52 = np.max([a52, 0.9])
    a53 = np.max([a53, 1.0])
    if Z > 0.01:
        a52 = np.min([a52, 1.0])
        a53 = np.min([a53, 1.1])
    
    a54 = get_zdep_param(Z, 3.855707e-1, -6.104166e-1, 5.676742e0, 1.060894e1, 5.284014e0)
    a55 = get_zdep_param(Z, 3.579064e-1, -6.442936e-1, 5.494644e0, 1.054952e1, 5.280991e0)
    a56 = get_zdep_param(Z, 9.587587e-1, 8.777464e-1, 2.017321e-1, 0, 0)
    
    a57 = np.min([1.4, 1.5135 + 0.3769 * zeta]) # sic
    a57 = np.max([0.6355 - 0.4192 * zeta, np.max([1.25, a57])])
    
    a58 = get_zdep_param(Z, 4.907546e-1, -1.683928e-1, -3.108742e-1, -7.202918e-2, 0)
    a59 = get_zdep_param(Z, 4.537070e0, -4.465455e0, -1.612690e0, -1.623246e0, 0)
    a60 = get_zdep_param(Z, 1.796220e0, 2.814020e-1, 1.423325e0, 3.421036e-1, 0)
    a61 = get_zdep_param(Z, 2.256216e0, 3.773400e-1, 1.537867e0, 4.396373e-1, 0)
    a62 = get_zdep_param(Z, 8.4300e-2, -4.7500e-2, -3.5200e-2, 0, 0)
    a63 = get_zdep_param(Z, 7.3600e-2, 7.4900e-2, 4.4260e-2, 0, 0)
    a64 = get_zdep_param(Z, 1.3600e-1, 3.5200e-2, 0, 0, 0)
    a65 = get_zdep_param(Z, 1.564231e-3, 1.653042e-3, -4.439786e-3, -4.951011e-3, -1.216530e-3)
    a66 = get_zdep_param(Z, 1.4770e0, 2.9600e-1, 0, 0, 0)
    a67 = get_zdep_param(Z, 5.210157e0, -4.143695e0, -2.120870e0, 0, 0)
    a68 = get_zdep_param(Z, 1.1160e0, 1.6600e-1, 0, 0, 0)
    
    a62 = np.max([0.065, a62]) # sic
    if Z < 0.004:
        a63 = np.min([0.055, a63])
    a64 = np.max([0.091, np.min([0.121, a64])])
    a66 = np.max([a66, np.min([1.6, -0.308 - 1.046 * zeta])])
    a66 = np.max([0.8, np.min([0.8 - 2.0 * zeta, a66])])
    a68 = np.max([0.9, np.min([a68, 1.0])])
    if a68 > a66:
        a64 = (a58 * a66 ** a60) / (a59 + a66 ** a61)
    a68 = np.min([a68, a66])
    
    a69 = get_zdep_param(Z, 1.071489e0, -1.164852e-1, -8.623831e-2, -1.582349e-2, 0)
    a70 = get_zdep_param(Z, 7.108492e-1, 7.935927e-1, 3.926983e-1, 3.622146e-2, 0)
    a71 = get_zdep_param(Z, 3.478514e0, -2.585474e-2, -1.512955e-2, -2.833691e-3, 0)
    a72 = get_zdep_param(Z, 9.132108e-1, -1.653695e-1, 0, 3.636784e-2, 0)
    a73 = get_zdep_param(Z, 3.969331e-3, 4.539076e-3, 1.720906e-3, 1.897857e-4, 0)
    a74 = get_zdep_param(Z, 1.600e0, 7.640e-1, 3.322e-1, 0, 0)
    
    if Z > 0.01: # sic
        a72 = np.max([a72, 0.95])
    a74 = np.max([1.4, np.min([a74, 1.6])])
    
    a75 = get_zdep_param(Z, 8.109e-1, -6.282e-1, 0, 0, 0)
    a76 = get_zdep_param(Z, 1.192334e-2, 1.083057e-2, 1.230969e0, 1.551656e0, 0)
    a77 = get_zdep_param(Z, -1.668868e-1, 5.818123e-1, -1.105027e1, -1.668070e1, 0)
    a78 = get_zdep_param(Z, 7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2, 0)
    a79 = get_zdep_param(Z, 9.409838e0, 1.522928e0, 0, 0, 0)
    a80 = get_zdep_param(Z, -2.7110e-1, -5.7560e-1, -8.3800e-2, 0, 0)
    a81 = get_zdep_param(Z, 2.4930e0, 1.1475e0, 0, 0, 0)
    
    a75 = np.max([1.0, np.min([a75, 1.27])]) # sic
    a75 = np.max([a75, 0.6355 - 0.4192 * zeta])
    a76 = np.max([a76, -0.1015564 - 0.2161264 * zeta - 0.05182516 * zeta ** 2])
    a77 = np.max([-0.3868776 - 0.5457078 * zeta - 0.1463472 * zeta ** 2, np.min([0.0, a77])])
    a78 = np.max([0.0, np.min([a78, 7.454 + 9.046 * zeta])])
    a79 = np.min([a79, np.max([2.0, -13.3 - 18.6 * zeta])])
    a80 = np.max([0.0585542, a80])
    a81 = np.min([1.5, np.max([0.4, a81])])

    c1 = -8.67073e-2
    
    eta = 10 * np.ones(len(M))
    cond = (Z <= 0.0009) & (M >= 1.1)
    eta[cond] = 20
    cond = (Z <= 0.0009) & (M < 1.1) & (M > 1.0)
    eta[cond] = (10 * (1.1 - M[cond]) + 20 * (M[cond] - 1.0)) / 0.1
    
    Mhook = 1.0185 + 0.16015 * zeta + 0.0892 * zeta ** 2
    
    # Timescales
    tBGB = (a1 + a2 * M ** 4 + a3 * M ** 5.5 + M ** 7) / (a4 * M ** 2 + a5 * M ** 7)
    mu = np.max([0.5 * np.ones(len(M)), 1.0 - 0.01 * np.max([a6 / M ** a7, a8 + a9 / M ** a10], axis=0)], axis=0)
    thook = mu * tBGB
    x = np.max([0.95, np.min([0.95 - 0.03 * (zeta + 0.30103), 0.99])])
    tMS = np.max([thook, x * tBGB], axis=0)
    epsilon = 0.01
    tau1 = t / thook
    tau1[tau1 > 1.0] = 1.0
    tau2 = (t - (1.0 - epsilon) * thook) / (epsilon * thook)
    tau2[tau2 > 1.0] = 1.0
    tau2[tau2 < 0.0] = 0.0

    # Quantities at end of the MS
    Ltms = (a11 * M ** 3 + a12 * M ** 4 + a13 * M ** (a16 + 1.8)) / (a14 + a15 * M ** 5 + M ** a16)
    Rtms = np.zeros(len(M))
    
    cond = (M <= a17)
    Rtms[cond] = (a18 + a19 * M[cond] ** a21) / (a20 + M[cond] ** a22)
    cond = (M >= a17 + 0.1)
    Rtms[cond] = (c1 * M[cond] ** 3 + a23 * M[cond] ** a26 + a24 * M[cond] ** (a26 + 1.5)) / (a25 + M[cond] ** 5)
    Rtms_lower = (a18 + a19 * a17 ** a21) / (a20 + a17 ** a22)
    Rtms_upper = (c1 * (a17 + 0.1) ** 3 + a23 * (a17 + 0.1) ** a26 + a24 * (a17 + 0.1) ** (a26 + 1.5)) / (a25 + (a17 + 0.1) ** 5)
    cond = (M > a17) & (M < a17 + 0.1)
    Rtms[cond] = (Rtms_lower * (a17 + 0.1 - M[cond]) + Rtms_upper * (M[cond] - a17)) / 0.1
    cond = (M < 0.5)
    Rtms[cond] = np.max([Rtms[cond], 1.5 * Rzams[cond]], axis=0)
    
    # Define luminosity and radius alpha/beta parameters, and gamma
    alphaL = np.zeros(len(M))
    B = (a45 + a46 * 2.0 ** a48) / (2.0 ** 0.4 + a47 * 2.0 ** 1.9)
    cond = (M < 0.5)
    alphaL[cond] = a49
    cond = (M >= 0.5) & (M < 0.7)
    alphaL[cond] = a49 + 5.0 * (0.3 - a49) * (M[cond] - 0.5)
    cond = (M >= 0.7) & (M < a52)
    alphaL[cond] = 0.3 + (a50 - 0.3) * (M[cond] - 0.7) / (a52 - 0.7)
    cond = (M >= a52) & (M < a53)
    alphaL[cond] = a50 + (a51 - a50) * (M[cond] - a52) / (a53 - a52)
    cond = (M >= a53) & (M < 2.0)
    alphaL[cond] = a51 + (B - a51) * (M[cond] - a53) / (2.0 - a53)
    cond = (M >= 2.0)
    alphaL[cond] = (a45 + a46 * M[cond] ** a48) / (M[cond] ** 0.4 + a47 * M[cond] ** 1.9)
    
    betaL = a54 - a55 * M ** a56
    betaL[betaL < 0] = 0
    
    alphaR = np.zeros(len(M))
    B = (a58 * a66 ** a60) / (a59 + a66 ** a61)
    C = (a58 * a67 ** a60) / (a59 + a67 ** a61)
    cond = (M < 0.5)
    alphaR[cond] = a62
    cond = (M >= 0.5) & (M < 0.65)
    alphaR[cond] = a62 + (a63 - a62) * (M[cond] - 0.5) / 0.15
    cond = (M >= 0.65) & (M < a68)
    alphaR[cond] = a63 + (a64 - a63) * (M[cond] - 0.65) / (a68 - 0.65)
    cond = (M >= a68) & (M < a66)
    alphaR[cond] = a64 + (B - a64) * (M[cond] - a68) / (a66 - a68)
    cond = (M >= a66) & (M < a67)
    alphaR[cond] = (a58 * M[cond] ** a60) / (a59 + M[cond] ** a61)
    cond = (M > a67)
    alphaR[cond] = C + a65 * (M[cond] - a67)
    
    betaRp = np.zeros(len(M))
    B = a69 * 2.0 ** 3.5 / (a70 + 2.0 ** a71)
    C = a69 * 16.0 ** 3.5 / (a70 + 16.0 ** a71)
    cond = (M <= 1.0)
    betaRp[cond] = 1.06
    cond = (M > 1.0) & (M < a74)
    betaRp[cond] = 1.06 + (a72 - 1.06) * (M[cond] - 1.0) / (a74 - 1.06)
    cond = (M >= a74) & (M < 2.0)
    betaRp[cond] = a72 + (B - a72) * (M[cond] - a74) / (2.0 - a74)
    cond = (M >= 2.0) & (M <= 16.0)
    betaRp[cond] = a69 * M[cond] ** 3.5 / (a70 + M[cond] ** a71)
    cond = (M > 16.0)
    betaRp[cond] = C + a73 * (M[cond] - 16.0)
    betaR = betaRp - 1
    
    gamma = np.zeros(len(M))
    B = a76 + a77 * (1.0 - a78) ** a79
    cond = (M <= 1.0)
    gamma[cond] = a76 + a77 * (M[cond] - a78) ** a79
    
    if a75 == 1.0:
        C = B
    else:
        C = a80
        cond = (M > 1.0) & (M <= a75)
        gamma[cond] = B + (a80 - B) * ((M[cond] - 1.0) / (a75 - 1.0)) ** a81
    
    cond = (M > a75) & (M < a75 + 0.1)
    gamma[cond] = C - 10.0 * (M[cond] - a75) * C
    cond = (M >= a75 + 0.1)
    gamma[cond] = 0
    
    # delta_L and delta_R
    delta_L = np.zeros(len(M))
    B = np.min([a34 / a33 ** a35, a36 / a33 ** a37], axis=0)
    cond = (M <= Mhook)
    delta_L[cond] = 0.0
    cond = (M > Mhook) & (M < a33)
    delta_L[cond] = B * ((M[cond] - Mhook) / (a33 - Mhook)) ** 0.4
    cond = (M >= a33)
    delta_L[cond] = np.min([a34 / M[cond] ** a35, a36 / M[cond] ** a37], axis=0)
    
    delta_R = np.zeros(len(M))
    B = (a38 + a39 * 2.0 ** 3.5) / (a40 * 2.0 ** 3 + 2.0 ** a41) - 1.0
    cond = (M <= Mhook)
    delta_R[cond] = 0.0
    cond = (M > Mhook) & (M <= a42)
    delta_R[cond] = a43 * ((M[cond] - Mhook) / (a42 - Mhook)) ** 0.5
    cond = (M > a42) & (M < 2.0)
    delta_R[cond] = a43 + (B - a43) * ((M[cond] - a42) / (2.0 - a42)) ** a44
    cond = (M >= 2.0)
    delta_R[cond] = (a38 + a39 * M[cond] ** 3.5) / (a40 * M[cond] ** 3 + M[cond] ** a41) - 1.0
    
    # Calculate L and R finally
    tau = t / tMS
    
    log_LMS_LZAMS = alphaL * tau + betaL * tau ** eta
    log_LMS_LZAMS += (np.log10(Ltms / Lzams) - alphaL - betaL) * tau ** 2
    log_LMS_LZAMS += -delta_L * (tau1 ** 2 - tau2 ** 2)
    
    log_RMS_RZAMS = alphaR * tau + betaR * tau ** 10 + gamma * tau ** 40
    log_RMS_RZAMS += (np.log10(Rtms / Rzams) - alphaR - betaR - gamma) * tau ** 3
    log_RMS_RZAMS += -delta_R * (tau1 ** 3 - tau2 ** 3)
    
    L = Lzams * 10 ** log_LMS_LZAMS
    R = Rzams * 10 ** log_RMS_RZAMS
    
    # Modify radius in accordance with Tout et al. 1997 (low-mass MS stars can be "substantially degenerate" in regime below)
    X = 0.76 - 3.0 * Z
    R = np.max([R, 0.0258 * (1.0 + X) ** (5 / 3) * M ** (-1 / 3)], axis=0)
    
    return L, R

def SSE_MS_get_flux(M, Z, t, filttable):
    """
    Function which calculates main sequence magnitudes given a filttable
    
    M: array of masses given in solar masses
    Z: float absolute metallicity (Zsun = 0.02)
    t: age (Gyr)
    
    Returns a dictionary with keys = filter names, elements = magnitude arrays
    """
    lum_lsun, rad_rsun = SSE_MS_get_L_and_R(M, Z, t)
    
    lum = 3.848e33 * lum_lsun # in erg/s
    rad = 6.957e10 * rad_rsun # cm
    
    sigma = 5.67051e-5
    Teff_K = (lum / (4 * np.pi * sigma * rad ** 2)) ** 0.25
    
    if type(filttable) == str:
        filttable = load_filtertable(filttable)
    elif type(filttable) != pd.DataFrame:
        raise ValueError('filttable must be either str or pd.DataFrame')

    # Read filter files
    filtnames = np.array(filttable['filtname'])
    filtfuncs = [load_filter(filttable.loc[ii,'path']) for ii in range(len(filttable))]

    zp_spectralflux = 1e8 * np.array(filttable['zp_spectralflux[ERG/S/CM2/ANGSTROM]'])

    # Calculate magnitudes
    filtdict = {} # at end, format as {key: mag array, ...}
    
    for ii in range(len(filtnames)):
        mag_arr = np.nan * np.ones(len(M))
        
        # Get filter function information
        wavelength_cm = 1e-8 * np.array(filtfuncs[ii]['wavelength[ANGSTROM]']) # cm
        transmission = np.array(filtfuncs[ii]['transmission'])

        passband_wl = np.sum(0.5 * (transmission[1:] + transmission[:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]))

        # Use trapezoid rule to evaluate integral of filtfunc * Planck distribution
        planck = 2 * h * c ** 2 / (wavelength_cm.reshape((1, wavelength_cm.size)) ** 5 * (np.exp(h * c / (k * np.outer(Teff_K, wavelength_cm))) - 1))

        planck_weighted = planck * transmission.reshape((1, transmission.size))
        integrated_planck_weighted = np.sum(0.5 * (planck_weighted[:,1:] + planck_weighted[:,:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]), axis=1)

        luminosity_cgs = 4 * np.pi ** 2 * rad ** 2 * integrated_planck_weighted

        spectral_lum = luminosity_cgs / (4 * np.pi * 3.086e19 ** 2 * passband_wl)

        # Calculate magnitudes (exclude black holes)
        mag_arr = -2.5 * np.log10(spectral_lum / zp_spectralflux[ii])
        
        filtdict[filtnames[ii]] = mag_arr

    return filtdict
    
def find_MS_TO(t, z):
    """
    Iteratively calculate main sequence mass turnoff as a function of cluster age
    and metallicity.
    
    Parameters
    ----------
    t: float
        age of cluster (in Gyr)
    
    z: float
        cluster metallicity
        
    Returns
    -------
    mto: float
        turn-off mass
    """
    t *= 1000
    
    # Make a grid of masses
    m_arr = np.logspace(-2, 3, 1000)
    t_arr = np.array([find_t_ms(z, m) for m in m_arr])
    
    # Interpolate t as a function of m to find turn-off mass
    interp = scipy.interpolate.interp1d(t_arr, m_arr)
    mto = float(interp(t))

    return mto

class Snapshot:
    """
    Snapshot class for snapshot file, usually something like 'initial.snap0137.dat.gz' 
    or 'initial.snapshots.h5', paired alongside conversion file and, preferably, 
    distance and metallicity.

    Parameters
    ----------
    fname: str
        filename of snapshot

    conv: str, dict, or pd.DataFrame
        if str, filename of unitfile (e.g., initial.conv.sh)
        if dict, dictionary of unit conversion factors
        if pd.DataFrame, table corresponding to initial.conv.sh

    snapshot_name: str
        key name for h5 snapshots; if unspecified, defaults to the last snapshot

    dist: float (default: None)
        distance to cluster in kpc
        
    z: float (default: None)
        cluster metallicity

    Attributes
    ----------
    data: pd.DataFrame
        snapshot table

    unitdict: dict
        Dictionary containing unit conversion information

    dist: float (None)
        distance to cluster in kpc

    age: float (None)
        age of cluster in Gyr 

    filtertable: pd.DataFrame
        table containing information about all filter for which photometry exists
    """
    def __init__(self, fname, conv, snapshot_name=None, dist=None, z=None):
        self.dist = dist
        self.z = z

        # Can either load old gzip snapshots or new hdf5 snapshots
        if '.dat.gz' in fname:
            # For gzip, read in snapshot as long list where each line is a str
            f = gzip.open(fname,'rb')
            text = str(f.read()).split('\\n')
            f.close()

            # Extract column names
            colrow = text[1].replace('#',' ').replace('  ',' ').split()
            colnames = [ colrow[ii][len(str(ii+1))+1:].replace(':','') for ii in range(len(colrow)) ]

            # Extract snapshot time
            t_snapshot = float(text[0].split('#')[1].split('=')[1].split()[0])

            # Make a list of lists each of which contains the contents of each object
            text = text[2:-1]
            rows = np.array([ np.array(row.split())[:len(colnames)] for row in text ])

            rows[np.where(rows == 'na')] = 'nan'
            rows = rows.astype(float)

            # Build a dictionary and cast to pandas DataFrame object
            table = {}
            for ii in range(len(colnames)):
                table[colnames[ii]] = rows[:, ii]

            self.data = pd.DataFrame(table)

        elif '.h5' in fname:
            # Take the last snapshot from the file if not specified
            if snapshot_name == None:
                snapshot_name = list(h5py.File(fname,'r').keys())[-1]

            # New version of CMC prints out pandas DataFrame formatted hdf5 files
            self.data = pd.read_hdf(fname,key=snapshot_name)

            # Extract snapshot time
            t_snapshot = float(snapshot_name.split('=')[1])

        else:
            raise ValueError('unsupported snapshot file type')

        # Now, build conversion dictionary
        if (type(conv) == str) or (type(conv) == pd.DataFrame):
            if type(conv) == str:
                f = open(conv, 'r')
                convfile = f.read().split('\n')
                f.close()

            # Produce unit dictionary
            self.unitdict = make_unitdict(convfile)

        elif type(conv) == dict:
            self.unitdict = conv

        else:
            raise ValueError('convfile must be either str or pd.DataFrame')

        # Also read in the time of the snapshot (code units) and convert to gyr
        self.age = self.convert_units(t_snapshot, 'code', 'gyr')

        self.filtertable = pd.DataFrame({'filtname': [],
                                             'path': [],
                              'zp_spectralflux[JY]': []})

    def convert_units(self, arr, in_unit='code', out_unit='code'):
        """
        Converts an array from CODE units to 'unit' using conversion factors specified
        in unitfile.

        Note: 's_' preceding an out_unit refers to 'stellar' quantities. 'nb_' refers
        to n-body units. Without these tags, it is presumed otherwise.

        Parameters
        ----------
        arr: array-like
            array to be converted

        in_unit: str (default: 'code')
            unit from which arr is to be converted

        out_unit: str (default: 'code')
            unit to which arr is to be converted

        Returns
        -------
        converted: array-like
            converted array
        """
        # Make sure both specified units are good
        if in_unit in self.unitdict.keys():
            ValueError('{} is not a recognized unit.'.format(in_unit))
        elif out_unit in self.unitdict.keys():
            ValueError('{} is not a recognized unit.'.format(out_unit))

        # Converted array
        converted = self.unitdict[out_unit] * arr / self.unitdict[in_unit]

        return converted

    def make_cuts(self, min_mass=None, max_mass=None, dmin=None, dmax=None, max_lum=None, fluxdict=None):
        """
        Helper method to return a boolean array where a given set of cuts are
        satisfied.

        Parameters
        ----------
        min_mass: float (default: None)
            If specified, only include stars above this mass

        min_mass: float (default: None)
            If specified, only include stars below this mass

        dmin: float (default: None)
            If specified, only include stars outside this projected radius

        dmax: float (default: None)
            If specified, only include stars inside this projected radius

        max_lum: float (default: None)
            IF specified, only include stars below this luminosity [LSUN]

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        Returns
        -------
        good: array-like of bool
            boolean array specifying where cuts are satisfied
        """
        # At the beginning, nothing is cut
        good = np.ones(len(self.data)).astype(bool)

        single = (self.data['binflag'] != 1)
        binary = (self.data['binflag'] == 1)

        # Mass cuts
        if min_mass is not None: # Pretend binaries are a single star with double mass
            good = good & ( ( (self.data['m[MSUN]'] > min_mass)                          & single ) |
                            ( (self.data['m0[MSUN]'] + self.data['m1[MSUN]'] > min_mass) & binary ) )
        if max_mass is not None:
            good = good & ( ( (self.data['m[MSUN]'] < max_mass)                          & single ) |
                            ( (self.data['m0[MSUN]'] + self.data['m1[MSUN]'] < max_mass) & binary ) )

        # Cuts on projected radius (in d)
        if (dmin is not None) | (dmax is not None):
            if 'd[PC]' not in self.data.keys():
                self.make_2d_projection()

            d_pc_arr = np.array(self.data['d[PC]'])

            if dmin is not None:
                good = good & ( d_pc_arr > dmin )
            if dmax is not None:
                good = good & ( d_pc_arr < dmax )
                
        # Cut on luminosity
        if max_lum is not None:
            good = good & ( ( single & (self.data['luminosity[LSUN]'] < max_lum) ) |
                            ( binary & (self.data['bin_star_lum0[LSUN]'] + self.data['bin_star_lum1[LSUN]'] < max_lum) ) )

        # Make sure all of the filters are actually there
        if fluxdict is not None:
            if not np.in1d(np.array(list(fluxdict.keys())), self.filtertable['filtname']).all():
                raise ValueError('One or more filters specified do not have photometry in this table.')

            # Cut on (observed) magnitudes
            for filtname in fluxdict.keys():
                faint_cut = fluxdict[filtname][0]
                bright_cut = fluxdict[filtname][1]

                colname = 'obsMag_{}'.format(filtname)
                bincolname0 = 'bin_obsMag0_{}'.format(filtname)
                bincolname1 = 'bin_obsMag1_{}'.format(filtname)

                if faint_cut is not None:
                    good = good & ( ( (self.data[colname] < faint_cut)                                       & single ) |
                                    ( (add_mags(self.data[bincolname0], self.data[bincolname1]) < faint_cut) & binary ) )

                if bright_cut is not None:
                    good = good & ( ( (self.data[colname] > bright_cut)                                       & single ) |
                                    ( (add_mags(self.data[bincolname0], self.data[bincolname1]) > bright_cut) & binary ) )

        return good

    def add_photometry(self, filttable):
        """
        Function which, assuming black-body behavior, assigns observed magnitudes
        to stars in desired filters

        For each filter, adds the following columns (# = filter name):
        absMag_#: absolute magnitude in filter # for single (np.nan for binary or black hole)
        bin_absMag0_#: absolute magnitude in filter # for first star in binary (np.nan for single or black hole)
        bin_absMag1_#: absolute magnitude in filter # for second star in binary (np.nan for single or black hole)
        tot_absMag_#: total magnitude in filter #, same as absMag_# for singles and is the magnitude sum of a binary pair if binary

        If distance is given, also add:
        obsMag_#: observed magnitude in filter # for single (np.nan for binary or black hole)
        bin_obsMag0_#: observed magnitude in filter # for first star in binary (np.nan for single or black hole)
        bin_obsMag1_#: observed magnitude in filter # for second star in binary (np.nan for single or black hole)
        tot_obsMag_#: total observed magnitude in filter #, same as absMag_# for singles and is the magnitude sum of a binary pair if binary
        
        Code assumes calculation in energy units, photometry in energy units.

        Parameters
        ----------
        filttable: str or pd.DataFrame
            if str, path to filter table
            if pd.DataFrame, table containing information about filters (see function: load_filtertable)

        Returns
        -------
        none
        """
        # If Teff has not been calculated, calculate it
        if 'Teff[K]' not in self.data.keys():
            self.calc_Teff()

        if type(filttable) == str:
            filttable = load_filtertable(filttable)
        elif type(filttable) != pd.DataFrame:
            raise ValueError('filttable must be either str or pd.DataFrame')

        # Read filter files
        filtnames = np.array(filttable['filtname'])
        filtfuncs = [load_filter(filttable.loc[ii,'path']) for ii in range(len(filttable))]
        
        zp_spectralflux = self.convert_units(np.array(filttable['zp_spectralflux[ERG/S/CM2/ANGSTROM]']), 'erg/s/cm2/angstrom', 'erg/s/cm3')

        if self.dist is not None:
            distance_modulus = 5 * np.log10(self.dist / 0.01)

        # Calculate magnitudes
        for ii in range(len(filtnames)):
            self.data['absMag_' + filtnames[ii]] = np.nan * np.ones(len(self.data))
            self.data['bin_absMag0_' + filtnames[ii]] = np.nan * np.ones(len(self.data))
            self.data['bin_absMag1_' + filtnames[ii]] = np.nan * np.ones(len(self.data))

            # Get filter function information
            wavelength_cm = np.array(self.convert_units(filtfuncs[ii]['wavelength[ANGSTROM]'], 'angstrom', 'cm'))
            transmission = np.array(filtfuncs[ii]['transmission'])
            
            # Get normalization integral
            passband_wl = np.sum(0.5 * (transmission[1:] + transmission[:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]))

            # Use trapezoid rule to evaluate integral of filtfunc * Planck distribution
            Teff_K = self.data.loc[(self.data['binflag'] != 1) & (self.data['startype'] != 14), 'Teff[K]']
            planck = 2 * h * c ** 2 / (wavelength_cm.reshape((1, wavelength_cm.size)) ** 5 * (np.exp(h * c / (k * np.outer(Teff_K, wavelength_cm))) - 1))

            planck_weighted = planck * transmission.reshape((1, transmission.size))
            integrated_planck_weighted = np.sum(0.5 * (planck_weighted[:,1:] + planck_weighted[:,:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]), axis=1)

            rad_rsun = self.data.loc[(self.data['binflag'] != 1) & (self.data['startype'] != 14), 'radius[RSUN]']
            rad_cm = self.convert_units(rad_rsun, 'rsun', 'cm')
            luminosity_cgs = 4 * np.pi ** 2 * rad_cm ** 2 * integrated_planck_weighted

            spectral_lum = luminosity_cgs / (4 * np.pi * self.convert_units(10, 'pc', 'cm') ** 2 * passband_wl)

            # Calculate magnitudes (exclude black holes)
            self.data.loc[(self.data['binflag'] != 1) & (self.data['startype'] != 14), 'absMag_' + filtnames[ii]] = -2.5 * np.log10(spectral_lum / zp_spectralflux[ii])

            if self.dist is not None:
                self.data['obsMag_' + filtnames[ii]] = self.data['absMag_' + filtnames[ii]]
                self.data['obsMag_' + filtnames[ii]] += distance_modulus

            # Repeat this process for the first star in each binary
            Teff0_K = self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14), 'bin_Teff0[K]']
            planck0 = 2 * h * c ** 2 / (wavelength_cm.reshape((1, wavelength_cm.size)) ** 5 * (np.exp(h * c / (k * np.outer(Teff0_K, wavelength_cm))) - 1))

            planck_weighted0 = planck0 * transmission.reshape((1, transmission.size))
            integrated_planck_weighted0 = np.sum(0.5 * (planck_weighted0[:,1:] + planck_weighted0[:,:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]), axis=1)

            rad0_rsun = self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14), 'bin_star_radius0[RSUN]']
            rad0_cm = self.convert_units(rad0_rsun, 'rsun', 'cm')
            luminosity0_cgs = 4 * np.pi ** 2 * rad0_cm ** 2 * integrated_planck_weighted0

            spectral_lum0 = luminosity0_cgs / (4 * np.pi * self.convert_units(10, 'pc', 'cm') ** 2 * passband_wl)

            # Calculate magnitudes
            self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14), 'bin_absMag0_' + filtnames[ii]] = -2.5 * np.log10(spectral_lum0 / zp_spectralflux[ii])

            if self.dist is not None:
                self.data['bin_obsMag0_' + filtnames[ii]] = self.data['bin_absMag0_' + filtnames[ii]]
                self.data['bin_obsMag0_' + filtnames[ii]] += distance_modulus

            # Repeat this process for the second star in each binary
            Teff1_K = self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype1'] != 14), 'bin_Teff1[K]']
            planck1 = 2 * h * c ** 2 / (wavelength_cm.reshape((1, wavelength_cm.size)) ** 5 * (np.exp(h * c / (k * np.outer(Teff1_K, wavelength_cm))) - 1))

            planck_weighted1 = planck1 * transmission.reshape((1, transmission.size))
            integrated_planck_weighted1 = np.sum(0.5 * (planck_weighted1[:,1:] + planck_weighted1[:,:-1]) * (wavelength_cm[1:] - wavelength_cm[:-1]), axis=1)

            rad1_rsun = self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype1'] != 14), 'bin_star_radius1[RSUN]']
            rad1_cm = self.convert_units(rad1_rsun, 'rsun', 'cm')
            luminosity1_cgs = 4 * np.pi ** 2 * rad1_cm ** 2 * integrated_planck_weighted1

            spectral_lum1 = luminosity1_cgs / (4 * np.pi * self.convert_units(10, 'pc', 'cm') ** 2 * passband_wl)

            # Calculate magnitudes
            self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype1'] != 14), 'bin_absMag1_' + filtnames[ii]] = -2.5 * np.log10(spectral_lum1 / zp_spectralflux[ii])

            if self.dist is not None:
                self.data['bin_obsMag1_' + filtnames[ii]] = self.data['bin_absMag1_' + filtnames[ii]]
                self.data['bin_obsMag1_' + filtnames[ii]] += distance_modulus
                
            # Add total magnitude columns together
            self.data['tot_absMag_' + filtnames[ii]] = np.nan * np.ones(len(self.data))
            
            good_single = (self.data['binflag'] != 1) & (self.data['startype'] != 14)
            self.data.loc[good_single, 'tot_absMag_' + filtnames[ii]] = self.data.loc[good_single, 'absMag_' + filtnames[ii]]
            
            good_binary = (self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14) & (self.data['bin_startype1'] != 14)
            self.data.loc[good_binary, 'tot_absMag_' + filtnames[ii]] = add_mags(self.data.loc[good_binary, 'bin_absMag0_' + filtnames[ii]],
                                                                                 self.data.loc[good_binary, 'bin_absMag1_' + filtnames[ii]])
                                                                                 
            good0_bad1 = (self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14) & (self.data['bin_startype1'] == 14)
            self.data.loc[good0_bad1, 'tot_absMag_' + filtnames[ii]] = self.data.loc[good0_bad1, 'bin_absMag0_' + filtnames[ii]]
            
            good1_bad0 = (self.data['binflag'] == 1) & (self.data['bin_startype0'] == 14) & (self.data['bin_startype1'] != 14)
            self.data.loc[good1_bad0, 'tot_absMag_' + filtnames[ii]] = self.data.loc[good1_bad0, 'bin_absMag1_' + filtnames[ii]]
            
            if self.dist is not None:
                self.data['tot_obsMag_' + filtnames[ii]] = self.data['tot_absMag_' + filtnames[ii]]
                self.data['tot_obsMag_' + filtnames[ii]] += distance_modulus

            # Add filter to filtertable
            filterrow = pd.DataFrame({'filtname': [filtnames[ii]],
                                          'path': [filttable.loc[ii,'path']],
           'zp_spectralflux[ERG/S/CM2/ANGSTROM]': [filttable.loc[ii,'zp_spectralflux[ERG/S/CM2/ANGSTROM]']],
                                     })

            self.filtertable = self.filtertable.append(filterrow, ignore_index=True)

    def make_2d_projection(self, seed=0):
        """
        Appends to a snapshot table a column projecting stellar positions onto the
        2-dimensional plane of the sky (with randomly generated angles). This adds
        a column 'd' which is the projected radius.

        Adds the following columns:
        d : Projected radial distance (code units)
        x, y, z: Projected x, y, z coordinates (code units)
        vx, vy, vz: Projected x, y, z components of velocity (code units)
        vd: Projected radial component of velocity (code units)
        va: Projected azimuthal component of velocity (code units)
        r[PC] : Radial distance (pc)
        d[PC] : Projected radial distance (pc)
        x[PC], y[PC], z[PC]: Projected x, y, z coordinates (pc)
        vx[KM/S], vy[KM/S], vz[KM/S]: Projected x, y, z components of velocity (km/s)
        vd[KM/S]: Projected radial component of velocity (km/s)
        va[KM/S]: Projected azimuthal component of velocity (km/s)

        Parameters
        ----------
        seed: int (default: None)
            random seed, if None then don't set a seed
        """
        r_arr = self.data['r']
        vr_arr = self.data['vr']
        vt_arr = self.data['vt']

        # Generate needed random numbers
        if seed is not None:
            np.random.seed(seed)
        cosTheta = np.random.uniform(-1, 1, len(r_arr)) # Cosine of polar angle
        sinTheta = np.sqrt(1 - cosTheta ** 2)
        Phi = np.random.uniform(0, 2*np.pi, len(r_arr)) # Azimuthal angle
        vAngle = np.random.uniform(0, 2*np.pi, len(r_arr)) # Additional angle for tangential velocities

        # Project positions and convert to pc
        d_arr = r_arr * sinTheta
        x_arr = d_arr * np.cos(Phi)
        y_arr = d_arr * np.sin(Phi)
        z_arr = r_arr * cosTheta

        r_pc_arr = self.convert_units(r_arr, 'code', 'pc')
        d_pc_arr = self.convert_units(d_arr, 'code', 'pc')
        x_pc_arr = self.convert_units(x_arr, 'code', 'pc')
        y_pc_arr = self.convert_units(y_arr, 'code', 'pc')
        z_pc_arr = self.convert_units(z_arr, 'code', 'pc')

        # Project velocities and convert to km/s
        abs_v = np.hypot(vr_arr, vt_arr)
        dotTheta = np.cos(vAngle) * vt_arr / r_arr
        dotPhi = np.sin(vAngle) * vt_arr / (r_arr * sinTheta)

        vx_arr = vr_arr * sinTheta * np.cos(Phi) - r_arr * dotPhi * sinTheta * np.sin(Phi) - r_arr * dotTheta * cosTheta * np.cos(Phi)
        vy_arr = vr_arr * sinTheta * np.sin(Phi) + r_arr * dotPhi * sinTheta * np.cos(Phi) - r_arr * dotTheta * cosTheta * np.sin(Phi)
        vz_arr = vr_arr * cosTheta + r_arr * dotTheta * sinTheta
        vd_arr = vx_arr * np.cos(Phi) + vy_arr * np.sin(Phi) # projected radial component of velocity
        va_arr = -vx_arr * np.sin(Phi) + vy_arr * np.cos(Phi) # azimuthal component of velocity (line-of-sight tangential)

        vr_kms_arr = self.convert_units(vr_arr, 'code', 'nb_km/s')
        vx_kms_arr = self.convert_units(vx_arr, 'code', 'nb_km/s')
        vy_kms_arr = self.convert_units(vy_arr, 'code', 'nb_km/s')
        vz_kms_arr = self.convert_units(vz_arr, 'code', 'nb_km/s')
        vd_kms_arr = self.convert_units(vd_arr, 'code', 'nb_km/s')
        va_kms_arr = self.convert_units(va_arr, 'code', 'nb_km/s')

        # Add new columns to table
        self.data['d'] = d_arr
        self.data['x'], self.data['y'], self.data['z'] = x_arr, y_arr, z_arr
        self.data['vx'], self.data['vy'], self.data['vz'] = vx_arr, vy_arr, vz_arr
        self.data['vd'], self.data['va'] = vd_arr, va_arr

        self.data['r[PC]'] = r_pc_arr
        self.data['d[PC]'] = d_pc_arr
        self.data['x[PC]'], self.data['y[PC]'], self.data['z[PC]'] = x_pc_arr, y_pc_arr, z_pc_arr
        self.data['vx[KM/S]'], self.data['vy[KM/S]'], self.data['vz[KM/S]'] = vx_kms_arr, vy_kms_arr, vz_kms_arr
        self.data['vd[KM/S]'], self.data['va[KM/S]'] = vd_kms_arr, va_kms_arr

    def calc_Teff(self):
        """
        Calculate the effective temperature Teff for every (non-BH) star in the catalog.

        Adds the following columns:
        Teff[K] : Effective temperature of a single star, np.nan for a binary or black hole
        bin_Teff0[K] : Effective temperature for first star in a binary, np.nan for a single or black hole
        bin_Teff1[K] : Effective temperature for second star in a binary, np.nan for a single or black hole

        Parameters
        ----------
        none
        """
        # Add columns for Teff
        self.data['Teff[K]'] = np.nan * np.ones(len(self.data))
        self.data['bin_Teff0[K]'] = np.nan * np.ones(len(self.data))
        self.data['bin_Teff1[K]'] = np.nan * np.ones(len(self.data))

        # First, calculate Teff for non-BH singles
        single = (self.data['binflag'] != 1) & (self.data['startype'] != 14)
        lum_lsun = self.data.loc[single, 'luminosity[LSUN]']
        rad_rsun = self.data.loc[single, 'radius[RSUN]']

        lum = self.convert_units(lum_lsun, 'lsun', 'erg/s')
        rad = self.convert_units(rad_rsun, 'rsun', 'cm')

        self.data.loc[single, 'Teff[K]'] = (lum / (4 * np.pi * rad ** 2 * sigma)) ** 0.25

        # Next, calculate Teff for non-BH doubles... start with the first of the pair
        binary0 = (self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14)
        lum0_lsun = self.data.loc[binary0, 'bin_star_lum0[LSUN]']
        rad0_rsun = self.data.loc[binary0, 'bin_star_radius0[RSUN]']

        lum0 = self.convert_units(lum0_lsun, 'lsun', 'erg/s')
        rad0 = self.convert_units(rad0_rsun, 'rsun', 'cm')

        self.data.loc[(self.data['binflag'] == 1) & (self.data['bin_startype0'] != 14), 'bin_Teff0[K]'] = (lum0 / (4 * np.pi * rad0 ** 2 * sigma)) ** 0.25

        # Same as above but for the second star in each binary pair
        binary1 = (self.data['binflag'] == 1) & (self.data['bin_startype1'] != 14)
        lum1_lsun = self.data.loc[binary1, 'bin_star_lum1[LSUN]']
        rad1_rsun = self.data.loc[binary1, 'bin_star_radius1[RSUN]']

        lum1 = self.convert_units(lum1_lsun, 'lsun', 'erg/s')
        rad1 = self.convert_units(rad1_rsun, 'rsun', 'cm')

        self.data.loc[binary1, 'bin_Teff1[K]'] = (lum1 / (4 * np.pi * rad1 ** 2 * sigma)) ** 0.25

    def calc_surface_gravity(self):
        """
        Adds surface gravity columns to snapshot table.

        Adds the following columns:
        g[CM/S] : Surface gravity of a single star, np.nan for a binary
        bin_g0[CM/S] : Surface gravity for first star in a binary, np.nan for a single
        bin_g1[CM/S] : Surface gravity for second star in a binary, np.nan for a single

        Parameters
        ----------
        none

        Returns
        -------
        none
        """
        # Add columns for g
        self.data['g[CM/S2]'] = np.nan * np.ones(len(self.data))
        self.data['bin_g0[CM/S2]'] = np.nan * np.ones(len(self.data))
        self.data['bin_g1[CM/S2]'] = np.nan * np.ones(len(self.data))

        # Add surface gravity for single
        mass_msun = self.data.loc[self.data['binflag'] != 1, 'm[MSUN]']
        mass_g = self.convert_units(mass_msun, 'msun', 'g')
        rad_rsun = self.data.loc[self.data['binflag'] != 1, 'radius[RSUN]']
        rad_cm = self.convert_units(rad_rsun, 'rsun', 'cm')

        self.data.loc[self.data['binflag'] != 1, 'g[CM/S2]'] = G * mass_g / rad_cm ** 2

        # Add surface gravity for first binary
        mass0_msun = self.data.loc[self.data['binflag'] != 1, 'm0[MSUN]']
        mass0_g = self.convert_units(mass0_msun, 'msun', 'g')
        rad0_rsun = self.data.loc[self.data['binflag'] != 1, 'bin_star_radius0[RSUN]']
        rad0_cm = self.convert_units(rad0_rsun, 'rsun', 'cm')

        self.data.loc[self.data['binflag'] == 1, 'g0[CM/S2]'] = G * mass0_g / rad0_cm ** 2

        # Add surface gravity for second binary
        mass1_msun = self.data.loc[self.data['binflag'] != 1, 'm1[MSUN]']
        mass1_g = self.convert_units(mass1_msun, 'msun', 'g')
        rad1_rsun = self.data.loc[self.data['binflag'] != 1, 'bin_star_radius1[RSUN]']
        rad1_cm = self.convert_units(rad1_rsun, 'rsun', 'cm')

        self.data.loc[self.data['binflag'] == 1, 'g1[CM/S2]'] = G * mass1_g / rad1_cm ** 2

    def make_spatial_density_profile(self, bin_edges=None, min_mass=None, max_mass=None, fluxdict=None, startypes=startype_star):
        """
        Creates a spatial density profile

        Parameters
        ----------
        bin_edges: array-like
            Bin edges of radial density profile (if None, make 100 bins from min to max, logarithmic spacing)

        min_mass: float (default: None)
            If specified, only include stars above this mass

        min_mass: float (default: None)
            If specified, only include stars below this mass

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        bin_edges: array-like
            Bin edges

        profile: array-like
            Number density in each annulus

        e_profile: array-like
            Poisson error in number density in each annulus
        """
        if 'd[PC]' not in self.data.keys():
            self.make_2d_projection()

        # If bin_edges is not specified, use default
        if bin_edges is None:
            bin_edges = np.logspace( np.log10(np.min(self.data['d[PC]'])), np.log10(np.max(self.data['d[PC]'])), 100 )

        # Read columns and make cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict) # make cuts
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        mass_arr = np.array(self.data.loc[good, 'm[MSUN]'])
        d_pc_arr = np.array(self.data.loc[good, 'd[PC]'])

        count, _ = np.histogram(d_pc_arr, bin_edges)

        # Create spatial density profile by binning and dividing by annulus area
        area = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
        profile = count / area
        e_profile = np.sqrt(count) / area

        return bin_edges, profile, e_profile

    def velocity_dispersion(self, min_mass=None, max_mass=None, dmin=None, dmax=None, startypes=startype_star):
        """
        Calculate velocity dispersion in km/s.

        Parameters
        ----------
        min_mass: float (default: None)
            If specified, only include stars above this mass

        min_mass: float (default: None)
            If specified, only include stars below this mass

        dmin: float (default: None)
            If specified, only include stars outside this projected radius

        dmax: float (default: None)
            If specified, only include stars inside this projected radius

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        veldisp: float
            Velocity dispersion (km/s)

        e_veldisp: float
            Uncertainty in velocity dispersion (km/s)
        """
        if dmin == None:
            dmin = 0
        if dmax == None:
            dmax = np.inf

        _, veldisp, e_veldisp = self.make_velocity_dispersion_profile(bin_edges=[dmin, dmax], min_mass=min_mass, max_mass=max_mass, startypes=startypes)

        veldisp = veldisp[0]
        e_veldisp = e_veldisp[0]

        return veldisp, e_veldisp

    def make_velocity_dispersion_profile(self, bin_edges=None, min_mass=None, max_mass=None, fluxdict=None, startypes=startype_star):
        """
        Calculate velocity dispersion in km/s.

        Parameters
        ----------
        bin_edges: array-like (default: None)
            Bin edges of mass function (if None, make 100 bins from min to max, logarithmic spacing)

        min_mass: float (default: None)
            If specified, only include stars above this mass

        min_mass: float (default: None)
            If specified, only include stars below this mass

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        bin_edges: array-like
            Bin edges of mass function (pc)

        veldisp_profile: array-like
            Velocity dispersion profile (km/s)

        e_veldisp_profile: array-like
            Uncertainty in velocity dispersion profile (km/s)
        """
        if 'd[PC]' not in self.data.keys():
            self.make_2d_projection()

        # Define cut and find relevant arrays, only using relevant startypes
        # As long as one of a binary pair is a good startype, include it
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        mass_arr = np.array(self.data.loc[good, 'm[MSUN]'])
        d_pc_arr = np.array(self.data.loc[good, 'd[PC]'])
        v_kms_arr = np.array(np.hypot(self.data.loc[good, 'vx[KM/S]'], np.hypot(self.data.loc[good, 'vy[KM/S]'], self.data.loc[good, 'vz[KM/S]'])))

        # If bin_edges is not specified, use default
        if bin_edges is None:
            bin_edges = np.logspace( np.log10(np.min(self.data['d[PC]'])), np.log10(np.max(self.data['d[PC]'])), 100 )

        # Calculate velocity dispersion profile
        hist, _ = np.histogram(d_pc_arr, bin_edges)
        veldisp_profile = np.array([ np.sqrt(np.average(v_kms_arr[np.where( (d_pc_arr >= bin_edges[ii]) & (d_pc_arr < bin_edges[ii+1]) )] ** 2))
                                  for ii in range(len(bin_edges) - 1) ]) / np.sqrt(3)

        e_veldisp_profile = veldisp_profile / np.sqrt(2 * hist)

        return bin_edges, veldisp_profile, e_veldisp_profile

    def make_mass_function(self, bin_edges=None, dmin=None, dmax=None, fluxdict=None, startypes=startype_star):
        """
        Creates a mass function

        Parameters
        ----------
        bin_edges: array-like
            Bin edges of mass function (if None, make 100 bins from min to max, logarithmic spacing)

        dmin: float (default: None)
            If specified, only include stars outside this projected radius

        dmax: float (default: None)
            If specified, only include stars inside this projected radius

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        bin_edges: array-like
            Bin edges

        mf: array-like
            Mass function

        e_mf: array-like
            Mass function uncertainty
        """
        # If bin_edges is not specified, use default
        if bin_edges is None:
            bin_edges = np.logspace( np.log10(np.min(self.data['m[MSUN]'])), np.log10(np.max(self.data['m[MSUN]'])), 100 )

        # Make relevant cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(dmin=dmin, dmax=dmax, fluxdict=fluxdict)
        
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )

        mass_arr = np.array(self.data['m[MSUN]'])

        binary = np.array(self.data['binflag'] == 1)
        good_both_bin = np.array(np.in1d(bin_startype0_arr, startypes) & np.in1d(bin_startype1_arr, startypes))
        good_bin0 = np.array(np.in1d(bin_startype0_arr, startypes) & ~np.in1d(bin_startype1_arr, startypes))
        good_bin1 = np.array(~np.in1d(bin_startype0_arr, startypes) & np.in1d(bin_startype1_arr, startypes))
        
        mass_arr[good & binary & good_both_bin] = np.array(self.data.loc[good & binary & good_both_bin, 'm0[MSUN]']) + np.array(self.data.loc[good & binary & good_both_bin, 'm1[MSUN]'])
        mass_arr[good & binary & good_bin0] = np.array(self.data.loc[good & binary & good_bin0, 'm0[MSUN]'])
        mass_arr[good & binary & good_bin1] = np.array(self.data.loc[good & binary & good_bin1, 'm1[MSUN]'])
        
        mass_arr = mass_arr[good]

        # Obtain mass function
        count, _ = np.histogram(mass_arr, bin_edges)
        dm = bin_edges[1:] - bin_edges[:-1]

        mf = count / dm
        e_mf = np.sqrt(count) / dm

        return mf, e_mf

    def fit_mass_function_slope(self, init_guess=2.35, min_mass=None, max_mass=None, dmin=None, dmax=None, fluxdict=None, startypes=startype_star): # SCIPY DEPENDENCY
        """
        Fits the mass function slope in a bin-independent way

        Parameters
        ----------
        init_guess: float (default: 2.35)
            initial guess for the power law slope (if unspecified, guess Salpeter)

        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass

        dmin: float (default: None)
            If specified, only include stars outside this projected radius

        dmax: float (default: None)
            If specified, only include stars inside this projected radius

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        alpha: float
            Converged value for mass function slope (if failed fit, returns np.nan)
        """
        # Make relevant cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, dmin=dmin, dmax=dmax, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        mass_arr = np.array(self.data.loc[good, 'm[MSUN]'])

        # Define -1 * loglikelihood function, we wish to minimize this likelihood function
        # (in other words, minimize -1 * loglikelihood)
        if min_mass == None:
            min_mass = np.min(mass_arr)
        if max_mass == None:
            max_mass = np.inf

        N = len(mass_arr)
        minusloglike = lambda alpha: -1 * ( N * np.log((1 - alpha) / (max_mass ** (1-alpha) - min_mass ** (1-alpha))) - alpha * np.sum(np.log(mass_arr)) )

        # Perform fit
        res = scipy.optimize.minimize(minusloglike, init_guess)

        if res.success:
            alpha = res.x[0]
        else:
            alpha = np.nan

        return alpha

    def make_smoothed_number_profile(self, bins=80, min_mass=None, max_mass=None, fluxdict=None, startypes=startype_star, min_logr=-1.5):
        """
        Creates smoothed number density profile by smearing out stars probabilistically.

        Parameters
        ----------
        bins: int (default: 80)
            number of bins used for number density profile

        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        min_logr: float (default: -1.5)
            Minimum logarithmic radius in parsec

        Returns
        -------
        bin_center: array-like
            radial points at which profile is evaluated (in pc)

        profile: array-like
            array of number densities

        e_profile: array-like
            uncertainty in number density
        """
        # Define cut and find relevant arrays, only using relevant startypes
        # As long as one of a binary pair is a good startype, include it
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )

        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)

        # Probabilistically count stars in each bin
        bin_edges = np.logspace( min_logr, np.log10(np.max(r_pc_arr)), bins+1 )
        bin_edges = bin_edges.reshape(1, len(bin_edges))

        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )
        weight = np.sum(weight, axis=0)

        # Make profile
        bin_center = (bin_edges[0,:-1] + bin_edges[0,1:]) / 2

        area = np.pi * (bin_edges[0,1:] ** 2 - bin_edges[0,:-1] ** 2)
        profile = weight / area
        e_profile = np.sqrt(weight) / area

        return bin_center, profile, e_profile

    def make_smoothed_brightness_profile(self, filtname, bins=80, min_mass=None, max_mass=None, max_lum=None, fluxdict=None, startypes=startype_star, min_logr=-1.5):
        """
        Creates smoothed surface brightness profile by smearing out stars probabilistically.

        Parameters
        ----------
        filtname: str
            filter name

        bins: int (default: 80)
            number of bins used for number density profile

        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass
            
        max_lum: float (default: None)
            IF specified, only include stars below this luminosity [LSUN]

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        min_logr: float (default: -1.5)
            Minimum logarithmic radius in parsec

        Returns
        -------
        bin_center: array-like
            radial points at which profile is evaluated (in arcsec)

        profile: array-like
            array of surface brightness values (mag arcsec^-2)
        """
        assert self.dist is not None
        
        # Make relevant cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict, max_lum=max_lum)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)

        # Probabilistically count stars in each bin
        bin_edges = np.logspace( min_logr, np.log10(np.max(r_pc_arr)), bins+1 )
        bin_edges = bin_edges.reshape(1, len(bin_edges))

        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )

        # Pull magnitudes from specified filter (making sure binaries
        # are handled correctly)
        if filtname not in np.array(self.filtertable['filtname']): # make sure the filter actually has photometry
            raise ValueError('Given filter does not have photometry')

        mag = self.data.loc[good, 'obsMag_{}'.format(filtname)]
        binmag0 = self.data.loc[good, 'bin_obsMag0_{}'.format(filtname)]
        binmag1 = self.data.loc[good, 'bin_obsMag1_{}'.format(filtname)]
        startype0 = self.data.loc[good, 'bin_startype0']
        startype1 = self.data.loc[good, 'bin_startype1']
        binflag = self.data.loc[good, 'binflag']

        good_bin01 = (binflag == 1) & np.in1d(startype0, startypes) & np.in1d(startype1, startypes)
        good_bin0 = (binflag == 1) & np.in1d(startype0, startypes) & ~np.in1d(startype1, startypes)
        good_bin1 = (binflag == 1) & ~np.in1d(startype0, startypes) & np.in1d(startype1, startypes)
        mag[good_bin01] = add_mags(binmag0[good_bin01], binmag1[good_bin01])
        mag[good_bin0] = binmag0[good_bin0]
        mag[good_bin1] = binmag1[good_bin1]
        mag = mag.values.reshape(len(mag), 1)

        mag_weight = np.sum(10 ** (-mag / 2.5) * weight, axis=0)

        # Make profile (mag / arcsec^2)
        area = np.pi * (bin_edges[0,1:] ** 2 - bin_edges[0,:-1] ** 2)
        arcsec2 = self.convert_units(1, 'arcsec', 'rad') ** 2
        dist_pc = self.convert_units(self.dist, 'kpc', 'pc')
        angular_area = area / dist_pc ** 2

        bin_center = 206265 * ( (bin_edges[0,:-1] + bin_edges[0,1:]) / 2 ) / dist_pc
        profile = -2.5 * np.log10( (arcsec2 / angular_area) * mag_weight )

        return bin_center, profile

    def make_smoothed_veldisp_profile(self, bins=80, min_mass=None, max_mass=None, dmax=None, fluxdict=None, startypes=startype_star, min_logr=-1.5):
        """
        Creates smoothed velocity dispersion profile by smearing out stars probabilistically.

        Parameters
        ----------
        bins: int (default: 80)
            number of bins used for number density profile

        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass
            
        dmax: float (default: None)
            If specified, the outermost bin boundary is this value (otherwise, it's the max radial value)

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list
        
        min_logr: float (default: -1.5)
            Minimum logarithmic radius in parsec

        Returns
        -------
        bin_center: array-like
            radial points at which profile is evaluated (in pc)

        veldisp_profile: array-like
            array of velocity dispersions (in km/s)

        e_veldisp_profile: array-like
            uncertainty in velocity dispersions (in km/s)
        """
        # Make relevant cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)

        # Probabilistically count stars in each bin
        if dmax is not None:
            bin_edges = np.logspace( min_logr, np.log10(dmax), bins+1 )
        else:
            bin_edges = np.logspace( min_logr, np.log10(np.max(r_pc_arr)), bins+1 )
        bin_edges = bin_edges.reshape(1, len(bin_edges))

        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )

        # Read in velocities
        bin_center = (bin_edges[0,:-1] + bin_edges[0,1:]) / 2

        v_arr = np.array(np.hypot(self.data.loc[good, 'vt'], self.data.loc[good, 'vr']))
        v_kms_arr = self.convert_units(v_arr, 'code', 'nb_km/s')
        v_kms_arr = v_kms_arr.reshape(len(v_kms_arr), 1)

        # Calculate velocity dispersion & uncertainty
        veldisp_profile = np.sqrt(np.sum(weight * v_kms_arr ** 2, axis=0) / np.sum(3 * weight, axis=0))
        weight = np.sum(weight, axis=0)
        e_veldisp_profile = veldisp_profile / np.sqrt(2 * weight)

        return bin_center, veldisp_profile, e_veldisp_profile

    def calculate_renclosed(self, enclosed_frac=0.5, qty='mass', bins=200, min_mass=None, max_mass=None, fluxdict=None, startypes=startype_star):
        """
        Calculate radius enclosing some percentage of mass/light by probabilistically binning stars and interpolating the CDF

        Parameters
        ----------
        enclosed_frac: float between 0 and 1, exclusive (default: 0.5)
            fraction of enclosed mass at which radius is desired
            
        qty: str, either 'mass' or 'light' (default: mass)
            depending on option, either calculates half-mass or half-light radius
        
        bins: int (default: 80)
            number of bins used for number density profile

        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list

        Returns
        -------
        rhm: float
            half-mass radius
        """
        if (enclosed_frac >= 1) or (enclosed_frac <= 0):
            raise ValueError('enclosed_frac must be between 0 and 1, exclusive')
        
        # Find startypes and make cuts
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']
        
        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        
        # Calculate weighted mass distribution
        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)

        binflag = np.array(self.data.loc[good, 'binflag'])
        startype0 = np.array(self.data.loc[good, 'bin_startype0'])
        startype1 = np.array(self.data.loc[good, 'bin_startype1'])
        
        binary = (binflag == 1)
        startype0_ok = (np.in1d(startype0, startypes))
        startype1_ok = (np.in1d(startype1, startypes))
        
        if qty == 'mass':
            mass_arr = np.array(self.data.loc[good, 'm[MSUN]'])
            mass0_arr = np.array(self.data.loc[good, 'm0[MSUN]'])
            mass1_arr = np.array(self.data.loc[good, 'm1[MSUN]'])
            
            mass_arr[binary & startype0_ok & startype1_ok] = mass0_arr[binary & startype0_ok & startype1_ok] + mass1_arr[binary & startype0_ok & startype1_ok]
            mass_arr[binary & startype0_ok & ~startype1_ok] = mass0_arr[binary & startype0_ok & ~startype1_ok]
            mass_arr[binary & ~startype0_ok & startype1_ok] = mass1_arr[binary & ~startype0_ok & startype1_ok]
        elif qty == 'light':
            lum_arr = np.array(self.data.loc[good,'luminosity[LSUN]'])
            lum0_arr = np.array(self.data.loc[good,'bin_star_lum0[LSUN]'])
            lum1_arr = np.array(self.data.loc[good,'bin_star_lum1[LSUN]'])
            
            lum_arr[binary & startype0_ok & startype1_ok] = lum0_arr[binary & startype0_ok & startype1_ok] + lum1_arr[binary & startype0_ok & startype1_ok]
            lum_arr[binary & startype0_ok & ~startype1_ok] = lum0_arr[binary & startype0_ok & ~startype1_ok]
            lum_arr[binary & ~startype0_ok & startype1_ok] = lum1_arr[binary & ~startype0_ok & startype1_ok]
        else:
            raise ValueError("qty must be either 'mass' or 'light'")
            
        # Probabilistically count stars in each bin and calculate CDF
        bin_edges = np.logspace( -5, np.log10(np.max(r_pc_arr)), bins+1 )
        bin_edges = bin_edges.reshape(1, len(bin_edges))
        bin_center = (bin_edges[0,:-1] + bin_edges[0,1:]) / 2
    
        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )
        
        if qty == 'mass':
            bin_mass = np.sum(mass_arr.reshape(len(mass_arr), 1) * weight, axis=0)
            cumdist = np.cumsum(bin_mass) / np.sum(bin_mass)
        elif qty == 'light':
            bin_light = np.sum(lum_arr.reshape(len(lum_arr), 1) * weight, axis=0)
            cumdist = np.cumsum(bin_light) / np.sum(bin_light)
        
        # Interpolate bin_center AS A FUNCTION OF cumdist and calculate rhm
        interp = scipy.interpolate.interp1d(cumdist, bin_center)
        renclosed = float(interp(enclosed_frac))
        
        return renclosed

    def make_paramdict(self):
        """
        Create a dictionary with some useful quantities about the cluster.
        
        Parameters
        ----------
        none
        
        Returns
        -------
        paramdict: dict
            parameter dictionary
        """
        # Initialize parameter dictionary
        paramdict = {}

        # Get relevant columns
        binary = (self.data['binflag'] == 1)
        startype = np.array(self.data['startype'])
        bin_startype0 = np.array(self.data['bin_startype0'])
        bin_startype1 = np.array(self.data['bin_startype1'])
        ospin = np.array(self.data['ospin'])
        ospin0 = np.array(self.data['ospin0'])
        ospin1 = np.array(self.data['ospin1'])

        mass = np.array(self.data['m[MSUN]'])
        mass0 = np.array(self.data['m0[MSUN]'])
        mass1 = np.array(self.data['m1[MSUN]'])

        lum = np.array(self.data['luminosity[LSUN]'])
        lum0 = np.array(self.data['bin_star_lum0[LSUN]'])
        lum1 = np.array(self.data['bin_star_lum1[LSUN]'])

        dmdt0 = np.array(self.data['dmdt0'])
        dmdt1 = np.array(self.data['dmdt1'])

        # Get count values
        paramdict['Nsys'] = len(startype)
        paramdict['Nobj'] = np.sum(~binary) + 2 * np.sum(binary)
        paramdict['Nlum'] = np.sum(~binary & np.in1d(startype, startype_star)) + np.sum(binary & (np.in1d(bin_startype0, startype_star) | np.in1d(bin_startype1, startype_star)))
        paramdict['Nstar'] = np.sum(np.in1d(startype, startype_star)) + np.sum(np.in1d(bin_startype0, startype_star)) + np.sum(np.in1d(bin_startype1, startype_star))

        for ii in range(16):
            paramdict['N_{}'.format(ii)] = np.sum(np.in1d(startype, ii)) + np.sum(np.in1d(bin_startype0, ii)) + np.sum(np.in1d(bin_startype1, ii))

        # Counting binaries
        paramdict['Nbin'] = np.sum(binary)
        paramdict['Nbinl'] = np.sum(np.in1d(bin_startype0, startype_star) | np.in1d(bin_startype1, startype_star))
        paramdict['Nbins'] = np.sum(np.in1d(bin_startype0, startype_star) & np.in1d(bin_startype1, startype_star))
        paramdict['Nbinr'] = np.sum(np.in1d(bin_startype0, startype_remnant) | np.in1d(bin_startype1, startype_remnant))
        paramdict['Nbint'] = np.sum( (dmdt0 != 0) & ~np.isnan(dmdt0) )

        for ii in range(16):
            for jj in range(ii+1):
                if ii == jj:
                    paramdict['Nbin_{}_{}'.format(ii, jj)] = np.sum(np.in1d(bin_startype0, ii) & np.in1d(bin_startype1, ii))
                else:
                    paramdict['Nbin_{}_{}'.format(ii, jj)] = np.sum(np.in1d(bin_startype0, ii) & np.in1d(bin_startype1, jj)) + np.sum(np.in1d(bin_startype0, jj) & np.in1d(bin_startype1, ii))

        # Mass transferring binaries
        for ii in range(16):
            for jj in range(16):
                if ii == jj:
                    paramdict['Nbint_{}_on_{}'.format(ii, jj)] = np.sum(np.in1d(bin_startype0, ii) & np.in1d(bin_startype1, ii) & (dmdt0 != 0) & ~np.isnan(dmdt0))
                else:
                    paramdict['Nbint_{}_on_{}'.format(ii, jj)] = np.sum(np.in1d(bin_startype0, ii) & np.in1d(bin_startype1, jj) & (dmdt1 > 0) & ~np.isnan(dmdt0)) + np.sum(np.in1d(bin_startype0, jj) & np.in1d(bin_startype1, ii) & (dmdt0 > 0) & ~np.isnan(dmdt0))

        # Get mass quantities
        paramdict['Mtot'] = np.sum(mass * ~binary) + np.sum(mass0 * binary) + np.sum(mass1 * binary)
        paramdict['Mstar'] = np.sum(mass * np.in1d(startype, startype_star)) + np.sum(mass0 * np.in1d(bin_startype0, startype_star)) + np.sum(mass1 * np.in1d(bin_startype1, startype_star))
        paramdict['Mbh'] = np.sum(mass * np.in1d(startype, 14)) + np.sum(mass0 * np.in1d(bin_startype0, 14)) + np.sum(mass1 * np.in1d(bin_startype1, 14))

        paramdict['Ltot'] = np.sum(lum * np.in1d(startype, startype_star)) + np.sum(lum0 * np.in1d(bin_startype0, startype_star)) + np.sum(lum1 * np.in1d(bin_startype1, startype_star))

        # Count other exotica

        # Millisecond pulsars
        per = 2 * np.pi / self.convert_units(ospin[np.in1d(startype, 13)], '1/yr', 'Hz')
        per0 = 2 * np.pi / self.convert_units(ospin0[np.in1d(bin_startype0, 13)], '1/yr', 'Hz')
        per1 = 2 * np.pi / self.convert_units(ospin1[np.in1d(bin_startype1, 13)], '1/yr', 'Hz')
        per_list = list(per) + list(per0) + list(per1)

        if len(per_list) > 0:
            paramdict['min_ns_per'] = np.min(per_list)
            paramdict['Nmsp'] = np.sum(np.array(per_list) <= 0.01)
        else:
            paramdict['min_ns_per'] = np.nan
            paramdict['Nmsp'] = 0

        # Turn-off mass, but only calculate this for cluster ages > 10 Myr
        if (self.z is not None) and (self.age > 0.01):
            mto = find_MS_TO(self.age, self.z)
            paramdict['mto'] = mto
            paramdict['Nbs'] = np.sum((mass > mto) & np.in1d(startype, startype_ms)) + np.sum((mass0 > mto) & np.in1d(bin_startype0, startype_ms)) + np.sum((mass1 > mto) & np.in1d(bin_startype1, startype_ms))
        else:
            paramdict['mto'] = np.nan
            paramdict['Nbs'] = np.nan
            
        # Most massive object of each stellar type
        for ii in range(16):
            if paramdict['N_{}'.format(ii)] != 0:
                paramdict['mmax_{}'.format(ii)] = np.max(np.append(mass, np.append(mass0, mass1))[np.append(startype, np.append(bin_startype0, bin_startype1)) == ii])
                paramdict['mmin_{}'.format(ii)] = np.min(np.append(mass, np.append(mass0, mass1))[np.append(startype, np.append(bin_startype0, bin_startype1)) == ii])
            else:
                paramdict['mmax_{}'.format(ii)] = np.nan
                paramdict['mmin_{}'.format(ii)] = np.nan
        
        # Half-mass & half-light radii (calculated using only stars)
        paramdict['rhm'] = self.calculate_renclosed(enclosed_frac=0.5, qty='mass')
        paramdict['rhl'] = self.calculate_renclosed(enclosed_frac=0.5, qty='light')
        
        # Velocity dispersion for all startypes (velocity dispersion of any combination of startypes can be deduced)
        # & central velocity dispersion
        paramdict['sig'] = self.make_smoothed_veldisp_profile(bins=1)[1][0]
        paramdict['sighm'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhm'])[1][0]
        paramdict['sighl'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhl'])[1][0]
        
        paramdict['sig_star'] = self.make_smoothed_veldisp_profile(bins=1, startypes=startype_star)[1][0]
        paramdict['sighm_star'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhm'], startypes=startype_star)[1][0]
        paramdict['sighl_star'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhl'], startypes=startype_star)[1][0]
        
        paramdict['sig_ms'] = self.make_smoothed_veldisp_profile(bins=1, startypes=startype_ms)[1][0]
        paramdict['sighm_ms'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhm'], startypes=startype_ms)[1][0]
        paramdict['sighl_ms'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhl'], startypes=startype_ms)[1][0]
        
        if np.sum(np.in1d(startype, startype_giant)) + np.sum(np.in1d(bin_startype0, startype_giant)) + np.sum(np.in1d(bin_startype1, startype_giant)) != 0:
            paramdict['sig_ms'] = self.make_smoothed_veldisp_profile(bins=1, startypes=startype_giant)[1][0]
            paramdict['sighm_ms'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhm'], startypes=startype_giant)[1][0]
            paramdict['sighl_ms'] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhl'], startypes=startype_giant)[1][0]
        else:
            paramdict['sig_ms'] = np.nan
            paramdict['sighm_ms'] = np.nan
            paramdict['sighl_ms'] = np.nan
        
        for ii in range(16):
            if paramdict['N_{}'.format(ii)] != 0:
                paramdict['sig_{}'.format(ii)] = self.make_smoothed_veldisp_profile(bins=1, startypes=ii)[1][0]
                paramdict['sighm_{}'.format(ii)] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhm'], startypes=ii)[1][0]
                paramdict['sighl_{}'.format(ii)] = self.make_smoothed_veldisp_profile(bins=1, dmax=paramdict['rhl'], startypes=ii)[1][0]
            else:
                paramdict['sig_{}'.format(ii)] = np.nan
                paramdict['sighm_{}'.format(ii)] = np.nan
                paramdict['sighl_{}'.format(ii)] = np.nan
        
        paramdict['age'] = self.age
        
        return paramdict

    def binary_fraction(self, min_q=None, max_q=None, min_mass=None, max_mass=None, dmin=None, dmax=None, fluxdict=None, startypes=startype_star, bin_startypes=startype_star):
        """
        Calculate the binary fraction subject to some cuts.

        Parameters
        ----------
        min_q: float (default: None)
            If specified, only report binary fraction for binaries with q > qmin.
        
        max_q: float (default: None)
            If specified, only report binary fraction for binaries with q < qmax.
        
        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass
        
        dmin: float (default: None)
            If specified, only include stars outside this projected radius. This
            is not done by random projection.

        dmax: float (default: None)
            If specified, only include stars inside this projected radius. This
            is not done by random projection.

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None

        startypes: array-like (default: startype_star)
            If specified, only include startypes in this list for consideration
        
        bin_startypes: array-like (default: startype_star)
            If specified, only include startypes in this list as binaries

        Returns
        -------
        binfrac: float
            Binary fraction, defined as the fraction of detectable point sources which are
            actually binary systems (!= fraction of stars in binaries).
        
        e_binfrac: float
            Error in the binary fraction, assuming Poissoninan counting error.
        """
        # Define cut and find relevant arrays, only using relevant startypes
        # As long as one of a binary pair is a good startype, include it
        startype_arr = self.data['startype']
        bin_startype0_arr = self.data['bin_startype0']
        bin_startype1_arr = self.data['bin_startype1']

        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        good = good & ( ((self.data['binflag'] != 1) & np.in1d(startype_arr, startypes))
                      | ((self.data['binflag'] == 1) & (np.in1d(bin_startype0_arr, startypes) | np.in1d(bin_startype1_arr, startypes)) ) )
        
        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)
        
        # Probabilistically count stars in each bin
        if dmin is None:
            dmin = 0
        if dmax is None:
            dmax = np.max(r_pc_arr)
        
        bin_edges = np.array([dmin, dmax])
        bin_edges = bin_edges.reshape(1, len(bin_edges))
        
        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )
        
        # Sum up all weights, and all binary weights, and calculate the binary
        # fraction
        all_weight = np.sum(weight)
        
        bin_good = (self.data.loc[good, 'binflag'] == 1)
        bintype_good = np.in1d(bin_startype0_arr[good], bin_startypes) & np.in1d(bin_startype1_arr[good], bin_startypes)
        
        if (min_q is not None) or (max_q is not None):
            q_arr = np.array(self.data.loc[good, 'm1[MSUN]'] / self.data.loc[good, 'm0[MSUN]']) # mass ratio
            q_arr[q_arr > 1] **= -1
            
            if (min_q is not None) and (max_q is not None):
                good_q = (q_arr > min_q) & (q_arr < max_q)
            elif (min_q is not None):
                good_q = (q_arr > min_q)
            elif (max_q is not None):
                good_q = (q_arr < max_q)
            
            bin_weight = np.sum(weight[good_q & bintype_good & bin_good])
        else:    
            bin_weight = np.sum(weight[bintype_good & bin_good])
        
        binfrac = bin_weight / all_weight
        e_binfrac = np.sqrt(binfrac * (binfrac + 1) / all_weight)
        
        return binfrac, e_binfrac
    
    def binary_fraction_photometric(self, filttable, mag_filter, blue_filter, red_filter, min_q, max_q=None, min_mass=None, max_mass=None,
                                          dmin=None, dmax=None, fluxdict=None, primary_flux_faint=None, primary_flux_brite=None, color_pad=0.1):
        """
        Calculate the binary fraction subject to some cuts. Use cuts on the CMD.

        Parameters
        ----------
        filttable: pd.DataFrame
            table containing information about filter functions

        mag_filter: str (default: None)
            Filter name for the "mag filter" (y-axis of CMD)
        
        blue_filter: str (default: None)
            Filter name for the "blue filter"
        
        red_filter: str (default: None)
            Filter name for the "red filter"

        min_q: float (default: None)
            If specified, only report binary fraction for binaries with q > qmin.
        
        max_q: float (default: None)
            If specified, only report binary fraction for binaries with q < qmax.
        
        min_mass: float (default: None)
            If specified, only include stars above this mass

        max_mass: float (default: None)
            If specified, only include stars below this mass
        
        dmin: float (default: None)
            If specified, only include stars outside this projected radius. This
            is not done by random projection.

        dmax: float (default: None)
            If specified, only include stars inside this projected radius. This
            is not done by random projection.

        fluxdict: dict (default: None)
            If specified, makes upper and lower (observed) magnitude cuts in certain filters
            Should be in the following format:
            {'filtname1': [faint1, brite1], 'filtname2': [faint1, brite1], ...}

            If you don't have a cut, put None
        
        primary_flux_faint: float (default: None)
            If specified, makes faint magnitude cut in the magnitude filter only to the primary
            star, not the whole system. This is used to more closely match
            the cuts applied to the ACS Globular Cluster Survey in Milone et al. 2012.
            
            Give in absolute magnitude.
        
        primary_flux_brite: float (default: None)
            If specified, makes brite magnitude cut in the magnitude filter only to the primary
            star, not the whole system. This is used to more closely match
            the cuts applied to the ACS Globular Cluster Survey in Milone et al. 2012.
            
            Cuts on mag filter, interpolates wrt. blue_filter-red_filter color, cut between faint_mag
            and brite_mag. Note: You're probably not intending to use both this argument and fluxdict.
            
            Give in absolute magnitude.
        
        color_pad: float (default: 0.1)
            To be considered as a single star or binary of interest, enforce that a source
            is at most this many magnitudes redder than the turnoff. To remove, set to None.

        Returns
        -------
        binfrac: float
            Binary fraction, defined as the fraction of detectable point sources which are
            actually binary systems (!= fraction of stars in binaries).
        
        e_binfrac: float
            Error in the binary fraction, assuming Poissoninan counting error.
        """
        # Add filters to snapshot
        #self.add_photometry(filttable)
        
        # Define cut and find relevant arrays
        good = self.make_cuts(min_mass=min_mass, max_mass=max_mass, fluxdict=fluxdict)
        
        if (primary_flux_faint is not None) or (primary_flux_brite is not None):
            # Generate a mass grid to figure out what masses the specified magnitude
            # cuts refer to; do this by interpolating mass wrt. magnitude.
            # Also calculate the corresponding blue/red mags
            single_bool = (self.data['binflag'] != 1)
            M = np.geomspace(np.min(self.data.loc[single_bool, 'm[MSUN]']), find_MS_TO(self.age, self.z), 1000)
            fdict = SSE_MS_get_flux(M, self.z, self.age, filttable)
            
            mag_to_mass_interp = scipy.interpolate.interp1d(fdict[mag_filter], M)
            mag_to_blue_interp = scipy.interpolate.interp1d(fdict[mag_filter], fdict[blue_filter])
            mag_to_red_interp = scipy.interpolate.interp1d(fdict[mag_filter], fdict[red_filter])
        
            # For each the brite and faint mag, calculate a locus in CMD space corresponding to all
            # possible binaries with that as a primary
            q_arr = np.linspace(1e-2, 1, 1000)
        if primary_flux_faint is not None:
            good_flux = np.ones(len(self.data)).astype(bool)
            
            cut_faint_mass = float(mag_to_mass_interp(primary_flux_faint))
            cut_faint_blue = float(mag_to_blue_interp(primary_flux_faint))
            cut_faint_red = float(mag_to_red_interp(primary_flux_faint))
            
            faint_cut_secondary_dict = SSE_MS_get_flux(q_arr * cut_faint_mass, self.z, self.age, filttable)
            faint_cut_secondary_mag = faint_cut_secondary_dict[mag_filter]
            faint_cut_secondary_blue = faint_cut_secondary_dict[blue_filter]
            faint_cut_secondary_red = faint_cut_secondary_dict[red_filter]
            faint_cut_total_mag = add_mags(primary_flux_faint, faint_cut_secondary_mag)
            faint_cut_total_color = add_mags(cut_faint_blue, faint_cut_secondary_blue) - add_mags(cut_faint_red, faint_cut_secondary_red)
            faint_cut_total_mag = np.append(primary_flux_faint, faint_cut_total_mag)
            faint_cut_total_color = np.append(cut_faint_blue-cut_faint_red, faint_cut_total_color)
            
            faint_interp = scipy.interpolate.interp1d(faint_cut_total_mag, faint_cut_total_color)

            # Remove faint sources
            good_flux[self.data[f'tot_absMag_{mag_filter}'] > faint_cut_total_mag[0]] = False

            faint_locus = (self.data[f'tot_absMag_{mag_filter}'] < faint_cut_total_mag[0]) & (self.data[f'tot_absMag_{mag_filter}'] > faint_cut_total_mag[-1])

            good_flux[faint_locus] = (self.data.loc[faint_locus, f'tot_absMag_{blue_filter}']-self.data.loc[faint_locus, f'tot_absMag_{red_filter}'] < faint_interp(self.data.loc[faint_locus, f'tot_absMag_{mag_filter}']))
            
            good = good & good_flux
        if primary_flux_brite is not None:
            good_flux = np.ones(len(self.data)).astype(bool)
            
            cut_brite_mass = float(mag_to_mass_interp(primary_flux_brite))
            cut_brite_blue = float(mag_to_blue_interp(primary_flux_brite))
            cut_brite_red = float(mag_to_red_interp(primary_flux_brite))
            
            brite_cut_secondary_dict = SSE_MS_get_flux(q_arr * cut_brite_mass, self.z, self.age, filttable)
            brite_cut_secondary_mag = brite_cut_secondary_dict[mag_filter]
            brite_cut_secondary_blue = brite_cut_secondary_dict[blue_filter]
            brite_cut_secondary_red = brite_cut_secondary_dict[red_filter]
            brite_cut_total_mag = add_mags(primary_flux_brite, brite_cut_secondary_mag)
            brite_cut_total_color = add_mags(cut_brite_blue, brite_cut_secondary_blue) - add_mags(cut_brite_red, brite_cut_secondary_red)
            brite_cut_total_mag = np.append(primary_flux_brite, brite_cut_total_mag)
            brite_cut_total_color = np.append(cut_brite_blue-cut_brite_red, brite_cut_total_color)
            
            brite_interp = scipy.interpolate.interp1d(brite_cut_total_mag, brite_cut_total_color)
            
            # Remove brite sources
            good_flux[self.data[f'tot_absMag_{mag_filter}'] < brite_cut_total_mag[-1]] = False
            
            brite_locus = (self.data[f'tot_absMag_{mag_filter}'] < brite_cut_total_mag[0]) & (self.data[f'tot_absMag_{mag_filter}'] > brite_cut_total_mag[-1])

            good_flux[brite_locus] = (self.data.loc[brite_locus, f'tot_absMag_{blue_filter}']-self.data.loc[brite_locus, f'tot_absMag_{red_filter}'] > brite_interp(self.data.loc[brite_locus, f'tot_absMag_{mag_filter}']))

            good = good & good_flux
        
        # Color pad, if specified
        if color_pad is not None:
            single_bool = np.in1d(self.data['startype'], [0, 1])
            M = np.geomspace(np.min(self.data.loc[single_bool, 'm[MSUN]']), find_MS_TO(self.age, self.z), 1000)
            fdict = SSE_MS_get_flux(M, self.z, self.age, filttable)
            mag_arr = fdict[mag_filter]
            blue_arr = fdict[blue_filter]
            red_arr = fdict[red_filter]
            
            # Turn-off is just the bluest point on the model isochrone
            turnoff_color = np.min(blue_arr - red_arr)
            turnoff_mag = mag_arr[blue_arr-red_arr == turnoff_color][0]
            
            good_color = ((self.data[f'tot_absMag_{blue_filter}'] - self.data[f'tot_absMag_{red_filter}']) > turnoff_color - 0.1)
            good = good & good_color
        
        # Use SSE prescriptions
        m_lower = 0.08
        M = np.geomspace(m_lower, find_MS_TO(self.age, self.z), 1000)
        filtdict = SSE_MS_get_flux(M, self.z, self.age, filttable)
    
        if min_q is not None:
            M_min = min_q * M
            filtdict_min = SSE_MS_get_flux(M_min, self.z, self.age, filttable)
            
            mag_min = add_mags(filtdict[mag_filter], filtdict_min[mag_filter])
            blue_min = add_mags(filtdict[blue_filter], filtdict_min[blue_filter])
            red_min = add_mags(filtdict[red_filter], filtdict_min[red_filter])
            
            interp_min = scipy.interpolate.interp1d(mag_min, blue_min-red_min,
                                                    bounds_error=False, kind='linear', fill_value=np.nan) # lin. interp. color vs. mag
            red_good = (self.data.loc[good, f'tot_absMag_{blue_filter}']-self.data.loc[good, f'tot_absMag_{red_filter}'] > interp_min(self.data.loc[good, f'tot_absMag_{mag_filter}']))
        else:
            raise ValueError('Must have min_q > 0. For photometry-independent idenfication, use binary_fraction().')
        
        if max_q is not None:
            M_max = max_q * M
            filtdict_max = SSE_MS_get_flux(M_max, self.z, self.age, filttable)
            
            mag_max = add_mags(filtdict[mag_filter], filtdict_max[mag_filter])
            blue_max = add_mags(filtdict[blue_filter], filtdict_max[blue_filter])
            red_max = add_mags(filtdict[red_filter], filtdict_max[red_filter])
            
            interp_max = scipy.interpolate.interp1d(mag_max, blue_max-red_max,
                                                    bounds_error=False, kind='linear', fill_value=np.nan) # lin. interp. color vs. mag
            blue_good = (self.data.loc[good, f'tot_absMag_{blue_filter}']-self.data.loc[good, f'tot_absMag_{red_filter}'] < interp_max(self.data.loc[good, f'tot_absMag_{mag_filter}']))
        else:
            blue_good = np.ones(len(self.data.loc[good])).astype(bool)
        
        # Probabilistically count stars in each bin
        r_pc_arr = self.convert_units(self.data.loc[good, 'r'], 'code', 'pc')
        r_pc_arr = r_pc_arr.values.reshape(len(r_pc_arr), 1)
        
        if dmin is None:
            dmin = 0
        if dmax is None:
            dmax = np.max(r_pc_arr)
        
        bin_edges = np.array([dmin, dmax])
        bin_edges = bin_edges.reshape(1, len(bin_edges))
        
        inner = r_pc_arr ** 2 - bin_edges[:,:-1] ** 2
        outer = r_pc_arr ** 2 - bin_edges[:,1:] ** 2
        inner[inner < 0] = 0
        outer[outer < 0] = 0

        weight = r_pc_arr ** -1 * ( np.sqrt(inner) - np.sqrt(outer) )
        
        # Identify binaries
        bin_weight = np.sum(weight[red_good & blue_good])
        all_weight = np.sum(weight)
        
        binfrac = bin_weight / all_weight
        e_binfrac = np.sqrt(binfrac * (binfrac + 1) / all_weight)
        
        return binfrac, e_binfrac
    
    def get_blue_stragglers(self, mag_filter, blue_filter, red_filter, filttable, color_pad=0.1):
        """
        A function to identify observationally realistic blue stragglers.
        Procedure is to locate all stars brighter than the turn-off magnitude
        and bluer than the turnoff color optionally padded by some amount.
        Additionally, blue stragglers are restricted to singles which are
        main sequence or binaries containing at least one main sequence star
        and no giant (where such binaries are counted as a single star).
        
        Parameters
        ----------
        mag_filter: str (default: None)
            Filter name for the "mag filter" (y-axis of CMD)
        
        blue_filter: str (default: None)
            Filter name for the "blue filter"
        
        red_filter: str (default: None)
            Filter name for the "red filter"
        
        color_pad: float (default: 0.1)
            To be counted as an "observational" blue straggler, a system must
            be bluer than the turn-off color by at least color_pad.
        
        Returns
        -------
        bs_cat: pd.DataFrame
            Blue straggler catalog
        """
        # Locate the MSTO
        ##mto = find_MS_TO(self.age, self.z)

        ##good = np.where((self.data['m[MSUN]'] < mto) & (self.data['binflag'] != 1) & np.in1d(self.data['startype'], [0, 1]))
        ##ms_singles = self.data.loc[good]
        ##turnoff_mag = np.min(ms_singles[f'tot_obsMag_{mag_filter}'])
        ##turnoff_color = np.min(ms_singles[f'tot_obsMag_{blue_filter}']-ms_singles[f'tot_obsMag_{red_filter}'])
        
        single_bool = np.in1d(self.data['startype'], [0, 1])
        m_lower = 0.08#np.min([np.min(self.data.loc[single_bool, 'm[MSUN]']), 0.08])
        
        M = np.geomspace(m_lower, find_MS_TO(self.age, self.z), 1000)
        fdict = SSE_MS_get_flux(M, self.z, self.age, filttable)
        mag_arr = fdict[mag_filter]
        blue_arr = fdict[blue_filter]
        red_arr = fdict[red_filter]
        
        # Turn-off is just the bluest point on the model isochrone
        turnoff_color = np.min(blue_arr - red_arr)
        turnoff_mag = mag_arr[blue_arr - red_arr == turnoff_color][0]
        
        
        # Identify observational blue stragglers
        color_good = ((self.data[f'tot_absMag_{blue_filter}'] - self.data[f'tot_absMag_{red_filter}']) < turnoff_color - 0.1)
        mag_good = (self.data[f'tot_absMag_{mag_filter}'] < turnoff_mag)
        ms_singles_good = (np.in1d(self.data['startype'], [0, 1]) & (self.data['binflag'] != 1))
        ms_binaries_good = ((np.in1d(self.data['bin_startype0'], [0, 1]) | np.in1d(self.data['bin_startype1'], [0, 1]))
                          & (~np.in1d(self.data['bin_startype0'], [2, 3, 4, 5, 6, 7, 8, 9])
                          & ~np.in1d(self.data['bin_startype1'], [2, 3, 4, 5, 6, 7, 8, 9]))
                          & (self.data['binflag'] == 1))
        
        bs_cat = self.data[color_good & mag_good & (ms_singles_good | ms_binaries_good)]

        return bs_cat
