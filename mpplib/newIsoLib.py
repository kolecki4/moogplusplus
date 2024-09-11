import numpy as np
import sys
import os

from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter
from scipy.stats import linregress
from scipy.integrate import simpson
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

from astropy.io import fits

class stellarParameters:    
    def __init__(self,teff = None,steff = None,logg = None,slogg = None,mass = None,smass = None,rad = None,srad = None,lum = None,slum = None,metal = None, alpha = None):
        self.teff = teff
        self.steff=steff
        self.logg = logg
        self.slogg = slogg
        self.mass = mass
        self.smass = smass
        self.rad = rad
        self.srad = srad
        self.lum = lum
        self.slum = slum
        self.metal = metal
        self.alpha = alpha

def queryPhotometry(star):
    """
    Uses astroquery to return an array containing 
    absolute magnitude photometric data for a target

    Parameters
    ----------
    star : string
        The SIMBAD name of a star.

    Returns
    -------
    starmags : ndarray
        An array which compiles UVBRI, 2MASS, Gaia, and WISE 
        photometric data for the star, with named columns
        according to the band name.

    """
    customSimbad=Simbad()
    customSimbad.add_votable_fields('flux(U)','flux(V)','flux(B)','flux(R)','flux(I)','flux(J)','flux(H)','flux(K)', 'id(GAIA DR2)', 'id(GAIA DR1)', 'parallax', 'plx_error', 'id(WISEA)', 'id(WISE)')
    customSimbad.add_votable_fields('flux_error(U)','flux_error(V)','flux_error(B)','flux_error(R)','flux_error(I)','flux_error(J)','flux_error(H)','flux_error(K)')
    columns = [('U',float),('B',float),('Bp',float),('V',float),
           ('R',float),('G',float),('Rp',float),('I',float),('J',float),('H',float),
           ('K',float),
           ('W1',float),('W2',float),('W3',float),('W4',float)]
    
    starmags = np.zeros(1, dtype = columns)
    magerrs = np.zeros(1, dtype = columns)
    
    #Query Object
    sim=customSimbad.query_object(star)
    #print(sim["ID_GAIA_DR2"])
    if 'Gaia DR2' in sim['ID_GAIA_DR2'][0]:
        gaiaid=sim['ID_GAIA_DR2'][0].replace('Gaia DR2 ', '')
        gaia2Vizier = Vizier(columns = ["Source", "Plx", "e_Plx", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag"])
        gaia2Vizier.ROW_LIMIT = -1
        gai=gaia2Vizier.query_object(star, catalog = 'I/345/gaia2', radius = 120*u.arcsec, coordinate_frame = "J2016")[0]
        gai = gai[np.argwhere(gai["Source"] == int(gaiaid))[0][0]]
        if type(gai['Plx']) == np.float64:
            plx=gai['Plx']
            eplx = gai['e_Plx']

            
        try:
            fuck
            cbj=Vizier.query_object(star, 'I/347/gaia2dis')
            d = cbj['I/347/gaia2dis']['rest'][np.where(cbj[0]['Source'] == int(gaiaid))[0][0]]
            #print('Bailer-Jones+2018 distance found:',d)
        except:
            #print('Using 1/parallax for distance value')
            d = 1/(plx/1000)
            ed = d*eplx/plx
        starmags['G']=gai['Gmag'] - 5*np.log10(d) + 5
        starmags['Bp']=gai['BPmag'] - 5*np.log10(d) + 5
        starmags['Rp']=gai['RPmag'] - 5*np.log10(d) + 5
        magerrs['G']=(gai['e_Gmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Bp']=(gai['e_BPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Rp']=(gai['e_RPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
    
    elif 'Gaia DR3' in sim['ID_GAIA_DR2'][0]:
        gaiaid=sim['ID_GAIA_DR2'][0].replace('Gaia DR3 ', '')
        #print(gaiaid)
        gaia3Vizier = Vizier(columns = ["Source", "Plx", "e_Plx", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag"])
        gaia3Vizier.ROW_LIMIT = -1
        gai=gaia3Vizier.query_object(star, catalog = 'I/355/gaiadr3', radius = 120*u.arcsec, coordinate_frame = "J2016")[0]
        gai = gai[np.argwhere(gai["Source"] == int(gaiaid))[0][0]]
        if type(gai['Plx']) == np.float64:
            plx=gai['Plx']
            eplx = gai['e_Plx']


        try:
            fuck
            cbj=Vizier.query_object(star, 'I/347/gaia2dis')
            d = cbj['I/347/gaia2dis']['rest'][np.where(cbj[0]['Source'] == int(gaiaid))[0][0]]
            #print('Bailer-Jones+2018 distance found:',d)
        except:
            #print('Using 1/parallax for distance value')
            d = 1/(plx/1000)
            ed = d*eplx/plx

        starmags['G']=gai['Gmag'] - 5*np.log10(d) + 5
        starmags['Bp']=gai['BPmag'] - 5*np.log10(d) + 5
        starmags['Rp']=gai['RPmag'] - 5*np.log10(d) + 5
        magerrs['G']=(0.002755**2 + gai['e_Gmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Bp']=(0.00279**2 + gai['e_BPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Rp']=(0.003779**2 + gai['e_RPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
            
    else:
        #print("No GAIA parallax found. Using parallax from SIMBAD sourced from " + sim['PLX_BIBCODE'][0])
        #print('Using 1/parallax for distance value')
        d = 1/(sim['PLX_VALUE'][0]/1000)
        ed = d*sim['PLX_ERROR'][0]/sim['PLX_VALUE'][0]
        
        
    if 'WISEA' in sim['ID_WISEA'][0]:
        Vizier.ROW_LIMIT = -1
        wiseid = sim['ID_WISEA'][0].replace('WISEA ', '')
        wise=Vizier.query_object(star, catalog = 'II/328/allwise', radius = 120*u.arcsec)
        #starmags['W1']=wise['II/328/allwise']['W1mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W2']=wise['II/328/allwise']['W2mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W3']=wise['II/328/allwise']['W3mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W4']=wise['II/328/allwise']['W4mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        
        #magerrs['W1']=(wise['II/328/allwise']['e_W1mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W2']=(wise['II/328/allwise']['e_W2mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W3']=(wise['II/328/allwise']['e_W3mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W4']=(wise['II/328/allwise']['e_W4mag'][ np.where(wise[0]['AllWISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        
    elif 'WISE' in sim['ID_WISE'][0]:
        Vizier.ROW_LIMIT = -1
        wiseid = sim['ID_WISE'][0].replace('WISE ', '')
        wise=Vizier.query_object(star, catalog = 'II/311', radius = 120*u.arcsec)  
        #starmags['W1']=wise['II/311/wise']['W1mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W2']=wise['II/311/wise']['W2mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W3']=wise['II/311/wise']['W3mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
        #starmags['W4']=wise['II/311/wise']['W4mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]] - 5*np.log10(d) + 5
    
        #magerrs['W1']=(wise['II/311/wise']['e_W1mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W2']=(wise['II/311/wise']['e_W2mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W3']=(wise['II/311/wise']['e_W3mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        #magerrs['W4']=(wise['II/311/wise']['e_W4mag'][ np.where(wise[0]['WISE'] == wiseid)[0][0]]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        
    starmags['U']=sim['FLUX_U'][0] - 5*np.log10(d) + 5
    starmags['V']=sim['FLUX_V'][0] - 5*np.log10(d) + 5
    starmags['B']=sim['FLUX_B'][0] - 5*np.log10(d) + 5
    starmags['R']=sim['FLUX_R'][0] - 5*np.log10(d) + 5
    starmags['I']=sim['FLUX_I'][0] - 5*np.log10(d) + 5
    starmags['J']=sim['FLUX_J'][0] - 5*np.log10(d) + 5
    starmags['H']=sim['FLUX_H'][0] - 5*np.log10(d) + 5
    starmags['K']=sim['FLUX_K'][0] - 5*np.log10(d) + 5
    
    magerrs['U']=(sim['FLUX_ERROR_U'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['V']=(sim['FLUX_ERROR_V'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['B']=(sim['FLUX_ERROR_B'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['R']=(sim['FLUX_ERROR_R'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['I']=(sim['FLUX_ERROR_I'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['J']=(sim['FLUX_ERROR_J'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['H']=(sim['FLUX_ERROR_H'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    magerrs['K']=(sim['FLUX_ERROR_K'][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
    
    for band in ['J','H','K']:
        if magerrs[band] == 0:
            magerrs[band] = 0.01
    return starmags, magerrs

def MISTreadiso(age, metall):
    """
    Reads a Dartmouth theoretical isochrone into a numpy array
    
    Parameters
    ----------
    age : float
        The age in Gyr of the isochrone to read. (Rounded to the nearest available value)
    metall : float
        The metallicity ([Fe/H] relative to solar) of the isochrone to read.

    Returns
    -------
    combinediso : ndarray
        The requested isochrone with all synthetic photometry 
        and with all columns labelled.

    """
    isopath = 'isochrones/'
    
    if metall < -1.5:
        metall = -1.5
    if metall > 0.5:
        metall = 0.5
    
    if metall >= 0:                                               
        s1='p'                                                   
    else:                                                        
        s1='m'
    s2 = str(round(10*metall)/10).replace('.','').replace('-','')
    
    mstr = s1+s2
    
    outColumns = [('EEP',float),('M/Mo',float),('LogTeff',float),('LogG',float),
               ('LogL/Lo',float),('U',float),('V',float),('B',float),
               ('R',float),('I',float),('J',float),('H',float),
               ('K',float),('G',float),('Bp',float),('Rp',float),
               ('W1',float),('W2',float),('W3',float),('W4',float)]
    
    availAges = np.arange(5,10.11,0.05)
    
    ### GAIA/2MASS/UVBRI
    with open(isopath + mstr + 'G2U.iso.cmd', 'r') as f:
        logAge = np.log10(age)+9

        isoStartInds = []
        isoLen = -1
        nCols = -1
        
        targetIso = np.argmin(np.abs(availAges - logAge))
        
        j = 0
        while True:
            line = f.readline()
            j+=1
            if "number of isochrones =" in line:
                niso = int(line[-7:-1])
                break
            
        while True:
            line = f.readline()
            ### This denotes the start of an isochrone. Make note of where
            if 'number of EEPs, cols' in line:
                isoStartInds.append(j)
                isoLen = int(line[-9:-3])
                nCols = int(line[-3:-1])
            
            ### This implies we have reached the correct isochrone. Extract it
            if len(isoStartInds) == targetIso+1:
                
                ### Initialize the output array
                combinediso = np.zeros(isoLen, dtype = outColumns)
                ### Skip over unneeded text
                f.readline()
                
                ### Note the names of the columns in the MIST file
                columns = []
                for col in f.readline().split()[1:]:
                    columns.append(col)
                
                ### Line by line, put the data into the combined iso array
                for i in range(isoLen): 
                    line = f.readline().split()
                    combinediso[i] = (line[0],line[3],line[4],line[5],line[6],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[22],line[23],line[24],0,0,0,0)
                break
            j+=1
            
    ### WISE        
    with open(isopath + mstr + 'WISE.iso.cmd', 'r') as f:
    
        isoStartInds = []
        isoLen = -1
        nCols = -1
    
        j = 0
        while True:
            line = f.readline()
            j+=1
            if "number of isochrones =" in line:
                niso = int(line[-7:-1])
                break
            
        while True:
            line = f.readline()
            ### This denotes the start of an isochrone. Make note of where
            if 'number of EEPs, cols' in line:
                isoStartInds.append(j)
                isoLen = int(line[-9:-3])
                nCols = int(line[-3:-1])
            
            ### This implies we have reached the correct isochrone. Extract it
            if len(isoStartInds) == targetIso+1:
                ### Skip over unneeded text
                f.readline()
                
                ### Note the names of the columns in the MIST file
                columns = []
                for col in f.readline().split()[1:]:
                    columns.append(col)
                
                ### Line by line, put the data into the combined iso array
                for i in range(isoLen): 
                    line = f.readline().split()
                    combinediso[i][16],combinediso[i][17],combinediso[i][18],combinediso[i][19] = (line[9],line[10],line[11],line[12])
                break
            j+=1
        
    return combinediso

def gibbs_sampling(image, w_start, samples):

    image_pdf = image / image.sum()
    height, width = image_pdf.shape
    result = []
    w_current = w_start

    for _ in range(samples):
        # sample height
        h_given_w = image_pdf[:, w_current] / image_pdf[:, w_current].sum()
        h_current = np.random.choice(np.array(range(height)), size=1, p=h_given_w)[0]

        # sample width
        w_given_h = image_pdf[h_current, :] / image_pdf[h_current, :].sum()
        w_current = np.random.choice(np.array(range(width)), size=1, p=w_given_h)[0]

        result.append((h_current, w_current))

    return result

def MISTgetTandG(metal, ageRange, target = None, SIMBADphot = None):
    """
    Derives an effective temperature and surface gravity
    for a target by comparing photometry to a MIST 
    theoretical isochrone of a given age and metallicity.
    Also required is an age range which is used to prevent
    inaccuracies
    
    Parameters
    ----------
    metall : float
        The metallicity ([Fe/H] relative to solar) of the isochrone to read.
    ageRange : int
        Flag for MIST/Dartmouth isochrone age range
            0: MIST standard age grid log(Age [yr]) = [5, 7.4]
            1: MIST standard age grid log(Age [yr]) = [7, 9.4]
            2: MIST standard age grid log(Age [yr]) = [9, 10.1]
    target : string
        The SIMBAD name of a star.
    SIMBADphot : ndarray
        The returned value of the querySIMBADphot function.
        A timesaver if you're doing multiple calls of this function
        on the same star.
    
    Returns
    -------
    T : float
        Effective Temperature.
    G : float
        Log(g).
    T_err : float
        The uncertainty on T.
    G_err : float
        The uncertainty on G.
    
    """
    
    if SIMBADphot == None and target != None: 
        SIMBADphot=queryPhotometry(target)
        starmags=SIMBADphot[0]
        magerrs=SIMBADphot[1]
    elif SIMBADphot != None and target == None:
        starmags=SIMBADphot[0]
        magerrs=SIMBADphot[1]
    else:
        print("Exactly one of 'SIMBADphot' and 'target' may be used")
            
    if ageRange == '0':
        availages = 10**(np.arange(5,7.4,0.1)-9)
        linages = np.arange(.0002,.025,.0002)
    elif ageRange == '1':   
        availages = 10**(np.arange(7,9.4,0.1)-9)
        linages = np.arange(.01,2.5,.02)
    elif ageRange == '2':
        availages = 10**(np.arange(9,10.1,0.1)-9)
        linages = np.arange(.1,13,.1)
    else:
        ageRange = ageRange.replace('[','').replace(']','').split(',')
        ageRange = [float(i) for i in ageRange]
        availages = 10**(np.arange(5,10.11,0.05)-9)
        availages = availages[(availages > ageRange[0]) & (availages < ageRange[1])]
        linages = np.linspace(ageRange[0],ageRange[1],125)
    # availages = np.append(np.arange(1,5,0.25),np.arange(5,13.1,0.5))
    availages = 10**(np.arange(9,10.3,0.05)-9)
    # availages = 10**(np.arange(7,9.4,0.1)-9)
    # availages = 10**(np.arange(5,7.4,0.1)-9)
    
    #availages = 10**(np.arange(6,9,0.05)-9)
    
    # Let the program know what data we have that's good to use
    bands = np.array(['U', 'B', 'Bp', 'V', 'R', 'G', 'Rp', 'I', 'J', 'H', 'K', 'W1', 'W2', 'W3', 'W4'])
    waves = np.array([365,445,532,551,658,673,797,806,1250,1650,2150,3400,4600,12000,22000])
    
    
    goodPoints = (np.array([starmags[0][i] != np.zeros_like(starmags)[0][i] for i in range(len(starmags[0]))]) & np.array([magerrs[0][i] != np.zeros_like(magerrs)[0][i] for i in range(len(magerrs[0]))]))
    
    bands = bands[goodPoints]
    waves = waves[goodPoints]
    
    
    
    
    
    
    
    smags = np.array([])
    emags = np.array([])
    imags = np.array([])
    

    for i in bands:
        smags = np.append(smags, starmags[0][i])
        emags = np.append(emags, magerrs[0][i])
    

    
    teffs = np.zeros([len(availages),1500])
    loggs = np.zeros([len(availages),1500])
    mass = np.zeros([len(availages),1500])
    lum = np.zeros([len(availages),1500])
    rad = np.zeros([len(availages),1500])
    
    residuals = np.zeros(shape = [len(availages),1500])
    chi2 = np.zeros(shape = [len(availages),1500])+9e9
    
    for k in range(len(availages)):
        isoc = MISTreadiso(availages[k],metal)
        
        for i in range(len(isoc['EEP'])):
            imags = np.array([])
            
            for j in range(len(bands)):
                imags = np.append(imags,isoc[bands[j]][i])
        
            resids = smags - imags
            chi2[k,i] = np.sum((resids**2)/(emags)**2)/(len(emags) - 2)
            residuals[k,i] = sum(np.abs(resids))
            teffs[k,i] = 10**isoc['LogTeff'][i]
            loggs[k,i] = isoc['LogG'][i]
            mass[k,i] = isoc['M/Mo'][i]
            lum[k,i] = 10**isoc['LogL/Lo'][i]
            rad[k,i] = np.sqrt((6.67408e-11 * isoc['M/Mo'][i]*1.988e30)/((10**isoc['LogG'][i] /100)))/6.957e8
   
    minPoint = np.where(chi2 == np.min(chi2))
    minPoint = (minPoint[0][0], minPoint[1][0])
    
    isoc = MISTreadiso(availages[minPoint[0]],metal)
    
    imags = np.array([])
        
    for j in range(len(bands)):
        imags = np.append(imags,isoc[bands[j]][minPoint[1]])
    
    plt.plot(waves,imags)
    plt.errorbar(waves,smags, yerr = emags)
    plt.figure()
    plt.contourf(range(1500), availages, chi2/np.min(chi2), levels = range(20)*np.min(chi2))
    plt.colorbar()
    plt.xlim(80,200)
    
    avTeff = np.average(teffs[chi2 < 2*np.min(chi2)])
    avMass = np.average(mass[chi2 < 2*np.min(chi2)])
    avLogg = np.average(loggs[chi2 < 2*np.min(chi2)])
    avLum = np.average(lum[chi2 < 2*np.min(chi2)])
    avRad = np.average(rad[chi2 < 2*np.min(chi2)])

    stdTeff = np.std(teffs[chi2 < 2*np.min(chi2)])
    stdLogg = np.std(loggs[chi2 < 2*np.min(chi2)])
    stdMass = np.std(mass[chi2 < 2*np.min(chi2)])
    stdLum = np.std(lum[chi2 < 2*np.min(chi2)])
    stdRad = np.std(rad[chi2 < 2*np.min(chi2)])

    #print("Lowest Chi^2 Value: %.2f" % np.min(chi2))
    
    return stellarParameters(avTeff, stdTeff, avLogg, stdLogg, avMass, stdMass, avRad, stdRad, avLum, stdLum, metal)

def createLDFileName(teff,logg,metal,band):
    
    availmetals = np.array([1,0.5,-0.0,-0.5,-1,-1.5,-2,-3,-4])
    
    if teff > 4900 and teff < 5000:
        teff = 4900
    if teff >= 5000 and teff < 5100:
        teff = 5100


    filename = "PHOENIX/LD/"+band+"/Z"
    #if metal == 0: metal = -0.0
    filename += "%+.1f" % availmetals[np.argmin(np.abs(availmetals - metal))]
    filename += "/lte%05d" % round(teff, -2)
    filename += "%+4.2f" %(-round(logg*2,0)/2)
    filename += "%+3.1f" % availmetals[np.argmin(np.abs(availmetals - metal))]
    filename += ".PHOENIX-ACES-AGSS-COND-LIMBDARK-2011-"
    if band in ["U","B","V","R","I"]:
        filename += "Johnson_"
    filename += band
    filename += ".fits"
    return filename

def fitLDCurve(teff,logg,metal,band):
    ld = fits.open(createLDFileName(teff,logg,metal,band))
    fit = linregress(ld[1].data[ld[1].data > 0.1],ld[0].data[ld[1].data > 0.1]/np.max(ld[0].data))
    return np.nanmax([fit.slope,0])


def readAbundances(workDir, paramFile, elements):
    with open(workDir + paramFile) as f:
        lines = f.readlines()


    solarMOOG = [0,
           12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
            6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
            3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
            3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
            1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
            1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
           -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
            0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
            0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
           -5.00,-0.54,-5.00,-5.00,-5.00]



    dirList = os.listdir(workDir)
    abundances = np.array([])
    errorbars = np.array([])

    with open(workDir + paramFile) as f:
        lines = f.readlines()

    for i in elements:
        try:
            summary = np.genfromtxt(workDir + "%i" % i + "/abundanceSummary.txt", names = ["wav", "[X/H]", "v_broad", "chi2"], ndmin = 2)
            if np.max(summary["XH"]) - np.min(summary["XH"]) > 1.1:
                summary = summary[(summary["XH"] > 0.3 + np.min(summary["XH"])) & (summary["XH"] < np.max(summary["XH"])-0.3) ]

            if len(summary) > 15:
                summary = summary[summary["chi2"] < 3]
            
            abundances = np.append(abundances, np.median(summary["XH"]))
            errorbars = np.append(errorbars, np.std(summary["XH"])/np.sqrt(len(summary)))

        except:
            abundances = np.append(abundances, 0)
            errorbars = np.append(errorbars, 0)
    for i in range(len(elements)):
        if abundances[i] == 0:
            abundances[i] = solarMOOG[int(elements[i])]
    return elements, abundances, errorbars

def readInputMetallicity(paramFile):
    with open(paramFile) as f:
        lines = f.readlines()

    for line in lines:
        if "[M/H]" in line:
            return float(line[10:20]) 


def readInputAlphaEnhancement(paramFile):
    with open(paramFile) as f:
        lines = f.readlines()

    for line in lines:
        if "[alpha/M]" in line:
            return float(line[10:20]) 



solarMOOG = [0,
       12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
        6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
        3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
        3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
        1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
        1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
       -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
        0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
        0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
       -5.00,-0.54,-5.00,-5.00,-5.00]


