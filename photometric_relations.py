import numpy as np
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from extinction import ccm89
import logging, sys
import shapely.geometry as geom
from scipy import interpolate
from uncertainties import ufloat


def search_catalog_name(starname, vizier_object, catalog_list, catalog_name):
    result = vizier_object.query_object(starname, catalog = [catalog_name], radius = 10.*u.arcsec)
    if str(result) != 'Empty TableList':
        catalog_list.append(result)
    else:
        result = vizier_object.query_object(starname, catalog = [catalog_name], radius = 20.*u.arcsec)
        if str(result) != 'Empty TableList':
            catalog_list.append(result)

    return catalog_list


def search_catalog_coords(RA, DEC, vizier_object, catalog_list, catalog_name):
    result = vizier_object.query_region(SkyCoord(ra = RA, dec = DEC, unit = (u.deg, u.deg), frame = 'icrs'),\
                                        width = '10s', catalog = catalog_name)
    if str(result) != 'Empty TableList':
        catalog_list.append(result)
    else:
        result = vizier_object.query_region(SkyCoord(ra = RA, dec = DEC, unit = (u.deg, u.deg), frame = 'icrs'),\
                                            width = '20s', catalog = catalog_name)
        if str(result) != 'Empty TableList':
            catalog_list.append(result)
    return catalog_list



def retrieve_catalogs(starname, vizier_object, use_coords = False, RA = None, DEC = None):
    result = []

    if use_coords:
        result = search_catalog_coords(RA, DEC, vizier_object, result, "II/246/out")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "I/311/hip2")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "I/259/tyc2")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "J/MNRAS/373/13/table1")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "J/ApJS/168/128/survey")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "II/215/catalog")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "V/130/gcs3")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "J/MNRAS/403/1949/ubvri")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "I/337/tgas")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "II/237/colors")
        result = search_catalog_coords(RA, DEC, vizier_object, result, "II/336/apass9")

    else:
        starname = starname.replace('A', '').replace('-', ' ').replace('_', ' ')

        result = search_catalog_name(starname, vizier_object, result, "II/246/out")
        result = search_catalog_name(starname, vizier_object, result, "I/311/hip2")
        result = search_catalog_name(starname, vizier_object, result, "I/259/tyc2")
        result = search_catalog_name(starname, vizier_object, result, "J/MNRAS/373/13/table1")
        result = search_catalog_name(starname, vizier_object, result, "J/ApJS/168/128/survey")
        result = search_catalog_name(starname, vizier_object, result, "II/215/catalog")
        result = search_catalog_name(starname, vizier_object, result, "V/130/gcs3")
        result = search_catalog_name(starname, vizier_object, result, "J/MNRAS/403/1949/ubvri")
        result = search_catalog_name(starname, vizier_object, result, "I/337/tgas")
        result = search_catalog_name(starname, vizier_object, result, "II/237/colors")
        result = search_catalog_name(starname, vizier_object, result, "II/336/apass9")


    return result



def vizier_params(starname, use_coords = False, RA = None, DEC = None):
    photometry = {}
    photometry['name'] = starname

    try:
        e_hpmag = 0.01
        v=Vizier(columns=["**"])
        if use_coords:
            photometry['RA'] = RA
            photometry['DEC'] = DEC

        result = retrieve_catalogs(starname, v, use_coords, RA, DEC)

        #if str(result) != 'Empty TableList':
        if len(result) != 0:
            #name_cats = result.keys()
            name_cats = [r.keys()[0] for r in result]
            #print name_cats

            if len(result[0][0].columns) < 5:
                v = Vizier(columns=["*"])
                result = retrieve_catalogs(starname, v, use_coords, RA, DEC)

            # 2MASS
            if "II/246/out" in name_cats:
                i = name_cats.index("II/246/out")

                jmag = result[i][0]['Jmag'][0]
                kmag = result[i][0]['Kmag'][0]
                hmag = result[i][0]['Hmag'][0]
                e_jmag = result[i][0]['e_Jmag'][0]
                e_kmag = result[i][0]['e_Kmag'][0]
                e_hmag = result[i][0]['e_Hmag'][0]

                f_errors = [0.01, 0.01, 0.01]

                for j,e in enumerate([e_jmag, e_kmag, e_hmag]):
                    if np.isnan(float(e)) == False:
                        f_errors[j] = e
                if np.isnan(float(jmag)) == False:
                    photometry['J'] = [jmag, f_errors[0], "II/246/out", 1.25*10000., 0.0]
                if np.isnan(float(kmag)) == False:
                    photometry['K'] = [kmag, f_errors[1], "II/246/out", 2.5*10000., 0.0]
                if np.isnan(float(hmag)) == False:
                    photometry['H'] = [hmag, f_errors[2], "II/246/out", 1.75*10000., 0.0]

                photometry['RA'] = [result[i][0]['RAJ2000'][0], "II/246/out"]
                photometry['DEC'] = [result[i][0]['DEJ2000'][0], "II/246/out"]

            # GAIA
            if "I/337/tgas" in name_cats:
                i = name_cats.index("I/337/tgas")
                plx = result[i][0]['Plx'][0]
                e_plx = result[i][0]['e_Plx'][0]
                photometry['parallax'] = [plx*u.mas, e_plx*u.mas, "I/337/tgas"]
                photometry['RA'] = [result[i][0]['RA_ICRS'][0], "I/337/tgas"]
                photometry['DEC'] = [result[i][0]['DE_ICRS'][0], "I/337/tgas"]

            else:

                # HIPPARCOS
                if "I/311/hip2" in name_cats:
                    i = name_cats.index("I/311/hip2")
                    plx = result[i][0]['Plx'][0]
                    e_plx = result[i][0]['e_Plx'][0]
                    photometry['parallax'] = [plx*u.mas, e_plx*u.mas, "I/311/hip2"]
                    photometry['RA'] = [result[i][0]['RArad'][0], "I/311/hip2"]
                    photometry['DEC'] = [result[i][0]['DErad'][0], "I/311/hip2"]
                    e_hpmag = result[i][0]['e_Hpmag'][0]

                    if e_hpmag == 0.0:
                        e_hpmag = 0.01

            # The Tycho-2 Catalogue (Hog+ 2000)
            if "I/259/tyc2" in name_cats:
                i = name_cats.index("I/259/tyc2")
                bmag = result[i][0]['BTmag'][0]
                e_bmag = result[i][0]['e_BTmag'][0]
                vmag = result[i][0]['VTmag'][0]
                e_vmag = result[i][0]['e_VTmag'][0]
                if e_bmag == 0.0: e_bmag = 0.01
                if e_vmag == 0.0: e_vmag = 0.01
                if np.isnan(float(bmag)) == False:
                    photometry['Bt'] = [bmag, e_bmag, "I/259/tyc2", 450.*10., 0.0]
                if np.isnan(float(vmag)) == False:
                    photometry['Vt'] = [vmag, e_vmag, "I/259/tyc2", 550.*10., 0.0]
                photometry['RA'] = [result[i][0]['RA_ICRS_'][0], "I/259/tyc2"]
                photometry['DEC'] = [result[i][0]['DE_ICRS_'][0], "I/259/tyc2"]


            # Koen et al. 2010
            if "J/MNRAS/403/1949/ubvri" in name_cats:
                i = name_cats.index("J/MNRAS/403/1949/ubvri")
                vmag = result[i][0]['Vmag'][0]
                BV = result[i][0]['B-V'][0]
                UB = result[i][0]['U-B'][0]
                VR = result[i][0]['V-Rc'][0]
                VI = result[i][0]['V-Ic'][0]
                if e_hpmag == 0.0: e_hpmag = 0.01
                if np.isnan(float(vmag)) == False:
                    photometry['V'] = [vmag, e_hpmag, "J/MNRAS/403/1949/ubvri", 550.*10., 0.0]
                if np.isnan(float(BV + vmag)) == False:
                    photometry['B'] = [BV + vmag, e_hpmag, "J/MNRAS/403/1949/ubvri", 450.*10., 0.0]
                if np.isnan(float(vmag - VR)) == False:
                    photometry['R'] = [vmag - VR, e_hpmag, "J/MNRAS/403/1949/ubvri", 675.*10., 0.0]
                if np.isnan(float(vmag - VI)) == False:
                    photometry['I'] = [vmag - VI, e_hpmag, "J/MNRAS/403/1949/ubvri", 875.*10., 0.0]

            else:
                # Casagrande et al. 2006
                if "J/MNRAS/373/13/table1" in name_cats:
                    i = name_cats.index("J/MNRAS/373/13/table1")
                    vmag = result[i][0]['Vmag'][0]
                    BV = result[i][0]['B-V'][0]
                    VR = result[i][0]['V-Rc'][0]
                    RI = result[i][0]['R-Ic'][0]
                    J = result[i][0]['Jmag'][0]
                    e_J = result[i][0]['e_Jmag'][0]
                    H = result[i][0]['Hmag'][0]
                    e_H = result[i][0]['e_Hmag'][0]
                    K = result[i][0]['Ksmag'][0]
                    e_K = result[i][0]['e_Ksmag'][0]
                    plx = result[i][0]['Plx'][0]
                    e_plx = result[i][0]['e_Plx'][0]
                    if e_hpmag == 0.0: e_hpmag = 0.01
                    if e_J == 0.0: e_J = 0.01
                    if e_H == 0.0: e_H = 0.01
                    if e_K == 0.0: e_K = 0.01

                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = [vmag, e_hpmag, "J/MNRAS/373/13/table1", 550.*10., 0.0]
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = [BV + vmag, e_hpmag, "J/MNRAS/373/13/table1", 450.*10., 0.0]
                    if np.isnan(float(vmag - VR)) == False:
                        photometry['R'] = [vmag - VR, e_hpmag, "J/MNRAS/373/13/table1", 675.*10., 0.0]
                    if np.isnan(float(vmag - VR - RI)) == False:
                        photometry['I'] = [vmag - VR - RI, e_hpmag, "J/MNRAS/373/13/table1", 875.*10., 0.0]
                    if np.isnan(float(J)) == False:
                        photometry['J'] = [J, e_J, "J/MNRAS/373/13/table1", 1.25*10000., 0.0]
                    if np.isnan(float(H)) == False:
                        photometry['H'] = [H, e_H, "J/MNRAS/373/13/table1", 1.75*10000., 0.0]
                    if np.isnan(float(K)) == False:
                        photometry['K'] = [K, e_K, "J/MNRAS/373/13/table1", 2.5*10000., 0.0]
                    if (np.isnan(float(plx)) == False) and (np.isnan(float(e_plx)) == False):
                        photometry['parallax'] = [plx*u.mas, e_plx*u.mas, "J/MNRAS/373/13/table1"]

                # Beers et al. 2007
                elif "J/ApJS/168/128/survey" in name_cats:
                    i = name_cats.index("J/ApJS/168/128/survey")
                    vmag = result[i][0]['Vmag'][0]
                    e_vmag = result[i][0]['e_Vmag'][0]
                    BV = result[i][0]['B-V'][0]
                    e_BV = result[i][0]['e_B-V'][0]
                    VR = result[i][0]['V-Rc'][0]
                    e_VR = result[i][0]['e_V-Rc'][0]
                    VI = result[i][0]['V-Ic'][0]
                    e_VI = result[i][0]['e_V-Ic'][0]
                    J = result[i][0]['Jmag'][0]
                    e_J = result[i][0]['e_Jmag'][0]
                    H = result[i][0]['Hmag'][0]
                    e_H = result[i][0]['e_Hmag'][0]
                    K = result[i][0]['Ksmag'][0]
                    e_K = result[i][0]['e_Ksmag'][0]
                    if e_vmag == 0.0: e_vmag = 0.01
                    if e_J == 0.0: e_J = 0.01
                    if e_H == 0.0: e_H = 0.01
                    if e_K == 0.0: e_k = 0.01
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = [vmag, e_vmag, "J/ApJS/168/128/survey", 550.*10., 0.0]
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = [BV + vmag, np.sqrt(e_vmag**2. + e_BV**2.), "J/ApJS/168/128/survey", 450.*10., 0.0]
                    if np.isnan(float(vmag - VR)) == False:
                        photometry['R'] = [vmag - VR, np.sqrt(e_vmag**2. + e_VR**2.), "J/ApJS/168/128/survey", 675.*10., 0.0]
                    if np.isnan(float(vmag - VI)) == False:
                        photometry['I'] = [vmag - VI, np.sqrt(e_vmag**2. + e_VI**2.), "J/ApJS/168/128/survey", 875.*10., 0.0]
                    if np.isnan(float(J)) == False:
                        photometry['J'] = [J, e_J, "J/ApJS/168/128/survey", 1.25*10000., 0.0]
                    if np.isnan(float(H)) == False:
                        photometry['H'] = [H, e_H, "J/ApJS/168/128/survey", 1.75*10000., 0.0]
                    if np.isnan(float(K)) == False:
                        photometry['K'] = [K, e_K, "J/ApJS/168/128/survey", 2.5*10000., 0.0]


            # HAUCK
            if "II/215/catalog" in name_cats:
                i = name_cats.index("II/215/catalog")
                vmag = result[i][0]['Vmag'][0]
                e_vmag = result[i][0]['e_Vmag'][0]
                by = result[i][0]['b-y'][0]
                e_by = result[i][0]['e_b-y'][0]
                m1 = result[i][0]['m1'][0]
                e_m1 = result[i][0]['e_m1'][0]
                c1 = result[i][0]['c1'][0]
                e_c1 = result[i][0]['e_c1'][0]
                if e_vmag == 0.0 or np.isnan(float(e_vmag)): e_vmag = 0.01
                if e_by == 0.0 or np.isnan(float(e_by)): e_by = 0.01
                if e_m1 == 0.0 or np.isnan(float(e_m1)): e_m1 = 0.01
                if e_c1 == 0.0 or np.isnan(float(e_c1)): e_c1 = 0.01
                if ('V' in photometry.keys()) == False:
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = [vmag, e_vmag, "II/215/catalog", 550.*10., 0.0]
                if np.isnan(float(by)) == False and np.isnan(float(e_by)) == False:
                    photometry['b'] = [by, e_by, "II/215/catalog", 450.*10., 0.0]
                    photometry['y'] = [0., e_by, "II/215/catalog", 550.*10., 0.0]
                if np.isnan(float(m1)) == False and np.isnan(float(e_m1)) == False:
                    photometry['m1'] = [m1, e_m1, "II/215/catalog"]
                if np.isnan(float(c1)) == False and np.isnan(float(e_c1)) == False:
                    photometry['c1'] = [c1, e_c1, "II/215/catalog"]


            else:
                # GENEVA
                if "V/130/gcs3" in name_cats:
                    i = name_cats.index("V/130/gcs3")
                    vmag = result[i][0]['Vmag'][0]
                    by = result[i][0]['b-y'][0]
                    m1 = result[i][0]['m1'][0]
                    c1 = result[i][0]['c1'][0]
                    if e_hpmag == 0.0: e_hpmag = 0.01
                    if ('V' in photometry.keys()) == False:
                        if np.isnan(float(vmag)) == False:
                            photometry['V'] = [vmag, e_hpmag, "V/130/gcs3", 550.*10., 0.0]
                    if np.isnan(float(by)) == False:
                        photometry['b'] = [by, e_hpmag, "V/130/gcs3", 450.*10., 0.0]
                        photometry['y'] = [0., e_hpmag, "V/130/gcs3", 550.*10., 0.0]
                    if np.isnan(float(m1)) == False:
                        photometry['m1'] = [m1, e_hpmag, "V/130/gcs3"]
                    if np.isnan(float(c1)) == False:
                        photometry['c1'] = [c1, e_hpmag, "V/130/gcs3"]

            photo_keys = photometry.keys()

            if (('R' in photo_keys) == False and ('I' in photo_keys) == False and \
                    ('V' in photo_keys) == False and ('B' in photo_keys) == False) or \
                    (('J' in photo_keys) == False and ('H' in photo_keys) == False and \
                    ('K' in photo_keys) == False) or \
                    (('J' in photo_keys) and (photometry['J'][1] > 1.0)) or \
                    (('H' in photo_keys) and (photometry['H'][1] > 1.0)) or \
                    (('K' in photo_keys) and (photometry['K'][1] > 1.0)) or \
                    (('R' in photo_keys) == False and ('I' in photo_keys) == False and \
                    ('B' in photo_keys) == False):


                # Ducati 2002
                if "II/237/colors" in name_cats:
                    i = name_cats.index("II/237/colors")
                    if e_hpmag == 0.0:
                        e_hpmag = 0.01
                    vmag = result[i][0]['Vmag'][0]
                    UV = result[i][0]['U-V'][0]
                    BV = result[i][0]['B-V'][0]
                    RV = result[i][0]['R-V'][0]
                    IV = result[i][0]['I-V'][0]
                    JV = result[i][0]['J-V'][0]
                    e_J = result[i][0]['Jsig'][0]
                    HV = result[i][0]['H-V'][0]
                    e_H = result[i][0]['Hsig'][0]
                    KV = result[i][0]['K-V'][0]
                    e_K = result[i][0]['Ksig'][0]
                    if np.isnan(float(UV + vmag)) == False:
                        photometry['U'] = [UV + vmag, e_hpmag, "II/237/colors", 350.*10., 0.0]
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = [vmag, e_hpmag, "II/237/colors", 550.*10., 0.0]
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = [BV + vmag, e_hpmag, "II/237/colors", 450.*10., 0.0]
                    if np.isnan(float(RV + vmag)) == False:
                        photometry['R'] = [RV + vmag, e_hpmag, "II/237/colors", 675.*10., 0.0]
                    if np.isnan(float(IV + vmag)) == False:
                        photometry['I'] = [IV + vmag, e_hpmag, "II/237/colors", 875.*10., 0.0]

                    if (('J' in photo_keys) == False and ('H' in photo_keys) == False and \
                       ('K' in photo_keys) == False) or \
                       (('J' in photo_keys) and (photometry['J'][1] > 1.0)) or \
                       (('H' in photo_keys) and (photometry['H'][1] > 1.0)) or \
                       (('K' in photo_keys) and (photometry['K'][1] > 1.0)):
                        if np.isnan(float(JV + vmag)) == False:
                            if np.isnan(float(e_J)) == False:
                                photometry['J'] = [JV + vmag, e_J*1e-2, "II/237/colors", 1.25*10000., 0.0]
                            else:
                                photometry['J'] = [JV + vmag, e_hpmag, "II/237/colors", 1.25*10000., 0.0]
                        if np.isnan(float(HV + vmag)) == False:
                            if np.isnan(float(e_H)) == False:
                                photometry['H'] = [HV + vmag, e_H*1e-2, "II/237/colors", 1.25*10000., 0.0]
                            else:
                                photometry['H'] = [HV + vmag, e_hpmag, "II/237/colors", 1.25*10000., 0.0]
                        if np.isnan(float(KV + vmag)) == False:
                            if np.isnan(float(e_K)) == False:
                                photometry['K'] = [KV + vmag, e_K*1e-2, "II/237/colors", 1.75*10000., 0.0]
                            else:
                                photometry['K'] = [KV + vmag, e_hpmag, "II/237/colors", 2.5*10000., 0.0]

                elif "II/336/apass9" in name_cats:
                    i = name_cats.index("II/336/apass9")
                    if e_hpmag == 0.0:
                        e_hpmag = 0.01
                    vmag = result[i][0]['Vmag'][0]
                    e_vmag = result[i][0]['e_Vmag'][0]
                    bmag = result[i][0]['Bmag'][0]
                    e_bmag = result[i][0]['e_Bmag'][0]
                    if e_vmag == 0.0:
                        e_vmag = 0.01
                    if e_bmag == 0.0:
                        e_bmag = 0.01
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = [vmag, e_vmag, "II/336/apass9", 550.*10., 0.0]
                    if np.isnan(float(bmag)) == False:
                        photometry['B'] = [bmag, e_bmag, "II/336/apass9", 450.*10., 0.0]

        del result

        #check that are the magnitudes are valid
        photo_keys = photometry.keys()
        non_magnitudes = ('parallax', 'name', 'RA', 'DEC')
        for mag in photo_keys:
            if mag not in non_magnitudes:
                if np.isnan(float(photometry[mag][0])) or np.isnan(float(photometry[mag][1])):
                    p = dict(photometry)
                    del p[mag]
                    photometry = p
            elif mag == 'parallax':
                if np.isnan(float(photometry[mag][0].value)) or np.isnan(float(photometry[mag][1].value)):
                    p = dict(photometry)
                    del p[mag]
                    photometry = p

    except Exception as e:
        #print 'Problem'
        #print "Unexpected error:", sys.exc_info()[0]
        #print e
        pass

    photometry = correct_extinction(photometry)

    return photometry


#******************************************************************************
#******************************************************************************

def correct_extinction(photometry):
    k = photometry.keys()
    photo = photometry
    mag = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'b', 'y', 'Bt', 'Vt']

    #print k
    #print photometry

    if ('RA' in k) and ('DEC' in k):
        ra = photometry['RA'][0]
        dec = photometry['DEC'][0]
        c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame = 'icrs')
        l = c.galactic.l.degree
        b = c.galactic.b.degree

        if 'parallax' in k:
            p = photometry['parallax'][0]
            d = p.to(u.parsec, equivalencies=u.parallax()).value
            Av = findext_arenou(l, b, distance = d)
        else:
            Av = findext_arenou(l, b)

        #print Av

        for m in mag:
            #print m
            if m in k:
                wave = photo[m][3]
                A_wave = ccm89(np.array([wave]), Av, 3.1, unit='aa')[0]
                #print A_wave, photometry[m][4]
                #print type(A_wave), type(photometry[m][4])
                photometry[m][4] = A_wave
                #print m, wave, A_wave

    #print photometry

    return photometry


'''
Function to retrieve the Av extinction by using the Arenou et al. 2002
interstellar maps.

From https://github.com/IvS-KULeuven/ComboCode/blob/master/cc/ivs/sed/extinctionmodels.py
'''


def findext_arenou(ll, bb, distance=None, redlaw='cardelli1989', Rv=3.1, norm='Av',**kwargs):

    # make sure that the values for b and l are within the correct range
    if bb < -90. or bb > 90:
        logging.error("galactic lattitude outside [-90,90] degrees")
    elif ll < 0. or ll > 360:
        logging.error("galactic longitude outside [0,360] degrees")
    elif distance < 0 and distance is not None:
        logging.error("distance is negative")

    # find the Arenou paramaters in the Appendix of Arenou et al. (1992)
    alpha, beta, gamma, rr0, saa = _getarenouparams(ll, bb)
    logging.info("Arenou params: alpha = %.2f, beta = %.2f, gamma = %.2f, r0 = %.2f and saa = %.2f" %(alpha, beta, gamma, rr0, saa))

    # compute the visual extinction from the Arenou paramaters using Equation 5
    # and 5bis
    if distance is None:
        av = alpha*rr0 + beta*rr0**2.
    else:
        distance = distance/1e3 # to kparsec
        if distance <= rr0:
            av = alpha*distance + beta*distance**2.
        else:
            av = alpha*rr0 + beta*rr0**2. + (distance-rr0)*gamma

    #-- Marshall is standard in Ak, but you can change this:
    #redwave, redflux = get_law(redlaw,Rv=Rv,norm=norm,photbands=['JOHNSON.V'],**kwargs)

    #return av/redflux[0]
    return av


def _getarenouparams(ll,bb):
    """
    Input: galactic coordinates
    Output: Arenou 1992 alpha, beta, gamma, rr0, saa

    @param ll: Galactic Longitude (in degrees)
    @type ll: float
    @param bb: Galactic Lattitude (in degrees)
    @type bb: float
    """
    if -90 < bb < -60:
        gamma = 0
        if   0   <= ll < 29:
            alpha = 2.22538 ; beta = -6.00212 ; rr0 = 0.052 ; saa = 13.
        elif 29  <= ll < 57:
            alpha = 3.35436 ; beta = -14.74567 ; rr0 = 0.035 ; saa = 40
        elif 57  <= ll < 85:
            alpha = 2.77637 ; beta = -9.62706 ; rr0 = 0.042 ; saa = 15
        elif 85  <= ll < 110:
            alpha = 4.44190 ; beta = -19.92097 ; rr0 = 0.025 ; saa = 36
        elif 110 <= ll < 150:
            alpha = 4.46685 ; beta = -26.07305 ; rr0 = 0.026 ; saa = 28
        elif 150 <= ll < 180:
            alpha = 7.63699 ; beta = -46.10856  ; rr0 = 0.014 ; saa = 38
        elif 180 <= ll < 210:
            alpha = 2.43412 ; beta = -8.69913 ; rr0 = 0.050 ; saa = 36
        elif 210 <= ll < 240:
            alpha = 3.34481 ; beta = -13.93228 ; rr0 = 0.035 ; saa = 38
        elif 240 <= ll < 270:
            alpha = 1.40733 ; beta = -13.43418  ; rr0 = 0.091 ; saa = 30
        elif 270 <= ll < 300:
            alpha = 1.64466 ; beta = -3.97380 ; rr0 = 0.074 ; saa = 28
        elif 300 <= ll < 330:
            alpha = 2.12696 ; beta = -6.05682 ; rr0 = 0.056 ; saa = 14
        elif 330 <= ll < 360:
            alpha = 2.34636 ; beta = -8.17784 ; rr0 = 0.052 ; saa = 16

    if -60 < bb < -45:
        gamma = 0
        if   0  <= ll < 30:
            alpha = 2.77060 ; beta = -9.52310 ; rr0 = 0.145 ; saa = 16
        elif 30 <= ll < 60:
            alpha = 1.96533 ; beta = -9.52310 ; rr0 = 0.174 ; saa = 06
        elif 60 <= ll < 110:
            alpha = 1.93622 ; beta = -13.31757 ; rr0 = 0.073 ; saa = 26
        elif 110 <= ll < 180:
            alpha = 1.05414 ; beta = -2.236540 ; rr0 = 0.223 ; saa = 74
        elif 180 <= ll < 210:
            alpha = 1.39990 ; beta = -1.35325 ; rr0 = 0.252 ; saa = 10
        elif 210 <= ll < 240:
            alpha = 2,73481 ; beta = -11.70266 ; rr0 = 0.117 ; saa = 8
        elif 240 <= ll < 270:
            alpha = 2.99784 ; beta = -11.64272 ; rr0 = 0.129 ; saa = 3
        elif 270 <= ll < 300:
            alpha = 3.23802 ; beta = -11.63810 ; rr0 = 0.139 ; saa = 7
        elif 300 <= ll < 330:
            alpha = 1.72679 ; beta = -6.05085 ; rr0 = 0.143 ; saa = 7
        elif 330 <= ll < 360:
            alpha = 1.88890 ; beta = -5.51861 ; rr0 = 0.171 ; saa = 14

    if -45 < bb < -30:
        gamma = 0
        if 0  <= ll < 30:
            alpha = 1.98973 ; beta = -4.86206  ; rr0 = 0.205 ; saa = 6
        elif 30  <= ll < 60:
            alpha = 1.49901 ; beta = -3.75837  ; rr0 = 0.199 ; saa = 16
        elif 60  <= ll < 90:
            alpha = 0.90091 ; beta = -1.30459  ; rr0 = 0.329 ; saa = 73
        elif 90  <= ll < 120:
            alpha = 1.94200 ; beta = -6.26833  ; rr0 = 0.155 ; saa = 18
        elif 120 <= ll < 160:
            alpha = -0.37804 ; beta = 10.77372  ; rr0 = 0.210 ; saa = 100
        elif 160 <= ll < 200:
            alpha = -0.15710 ; beta = 5.03190   ; rr0 = 0.294 ; saa = 24
        elif 200 <= ll < 235:
            alpha = 3.20162 ; beta = -10.59297 ; rr0 = 0.151 ; saa = 9
        elif 235 <= ll < 265:
            alpha = 1.95079 ; beta = 4.73280   ; rr0 = 0.206 ; saa = 21
        elif 265 <= ll < 300:
            alpha = 1.91200 ; beta = 4.97640   ; rr0 = 0.192 ; saa = 17
        elif 300 <= ll < 330:
            alpha = 2.50487 ; beta = -9.63106  ; rr0 = 0.145 ; saa = 28
        elif 330 <= ll < 360:
            alpha = 2.44394 ; beta = -9.17612  ; rr0 = 0.133 ; saa = 7

    if -30 < bb < -15:
        gamma = 0
        if 0 <= ll < 20:
            alpha = 2.82440 ; beta = -4.78217 ; rr0 = 0.295 ; saa = 32
        elif 20  <= ll < 40:
            alpha = 3.84362 ; beta = -8.04690 ; rr0 = 0.239 ; saa = 46
        elif 40  <= ll < 80:
            alpha = 0.60365 ; beta = 0.07893  ; rr0 = 0.523 ; saa = 22
        elif 80  <= ll < 100:
            alpha = 0.58307 ; beta = -0.21053 ; rr0 = 0.523 ; saa = 53
        elif 100 <= ll < 120:
            alpha = 2.03861 ; beta = -4.40843 ; rr0 = 0.231 ; saa = 60
        elif 120 <= ll < 140:
            alpha = 1.14271 ; beta = -1.35635 ; rr0 = 0.421 ; saa = 34
        elif 140 <= ll < 160:
            alpha = 0.79908 ; beta = 1.48074  ; rr0 = 0.513 ; saa = 61
        elif 160 <= ll < 180:
            alpha = 0.94260 ; beta = 8.16346  ; rr0 = 0.441 ; saa = 42
        elif 180 <= ll < 200:
            alpha = 1.66398 ; beta = 0.26775  ; rr0 = 0.523 ; saa = 42
        elif 200 <= ll < 220:
            alpha = 1.08760 ; beta = -1.02443 ; rr0 = 0.523 ; saa = 45
        elif 220 <= ll < 240:
            alpha = 1.20087 ; beta = -2.45407 ; rr0 = 0.245 ; saa = 6
        elif 240 <= ll < 260:
            alpha = 1.13147 ; beta = -1.87916 ; rr0 = 0.301 ; saa = 16
        elif 260 <= ll < 280:
            alpha = 0.97804 ; beta = -2.92838 ; rr0 = 0.338 ; saa = 21
        elif 290 <= ll < 300:
            alpha = 1.40086 ; beta = -1.12403 ; rr0 = 0.523 ; saa = 19
        elif 300 <= ll < 320:
            alpha = 2.06355 ; beta = -3.68278 ; rr0 = 0.280 ; saa = 42
        elif 320 <= ll < 340:
            alpha = 1.59260 ; beta = -2.18754 ; rr0 = 0.364 ; saa = 15
        elif 310 <= ll < 360:
            alpha = 1.45589 ; beta = -1.90598 ; rr0 = 0.382 ; saa = 21

    if -15 < bb < -5:
        gamma = 0
        if 0   <= ll < 10:
            alpha = 2.56438 ; beta = -2.31586 ; rr0 = 0.554 ; saa = 37
        elif 10  <= ll < 20:
            alpha = 3.24095 ; beta = -2.78217 ; rr0 = 0.582 ; saa = 38
        elif 20  <= ll < 30:
            alpha = 2.95627 ; beta = -2.57422 ; rr0 = 0.574 ; saa = 41 ; gamma = 0.08336
        elif 30  <= ll < 40:
            alpha = 1.85158 ; beta = -0.67716 ; rr0 = 1.152 ; saa = 4
        elif 40  <= ll < 50:
            alpha = 1.60770 ; beta = 0.35279  ; rr0 = 0.661 ; saa = 24
        elif 50  <= ll < 60:
            alpha = 0.69920 ; beta = -0.09146 ; rr0 = 0.952 ; saa = 2  ; gamma = 0.12839
        elif 60  <= ll < 80:
            alpha = 1.36189 ; beta = -1.05290 ; rr0 = 0.647 ; saa = 45 ; gamma = 0.16258
        elif 80  <= ll < 90:
            alpha = 0.33179 ; beta = 0.37629  ; rr0 = 1.152 ; saa = 62
        elif 90  <= ll < 100:
            alpha = 1.70303 ; beta = -0.75246 ; rr0 = 1.132 ; saa = 31
        elif 100 <= ll < 110:
            alpha = 1.97414 ; beta = -1.59784 ; rr0 = 0.618 ; saa = 35 ; gamma = 0.12847
        elif 110 <= ll < 120:
            alpha = 1.07407 ; beta = -0.40066 ; rr0 = 1.152 ; saa = 14 ; gamma = 0.17698
        elif 120 <= ll < 130:
            alpha = 1.69495 ; beta = -1.00071 ; rr0 = 0.847 ; saa = 28 ; gamma = 0.08567
        elif 130 <= ll < 140:
            alpha = 1.51449 ; beta = -0.08441 ; rr0 = 0.897 ; saa = 12
        elif 140 <= ll < 150:
            alpha = 1.87894 ; beta = -0.73314 ; rr0 = 1.152 ; saa = 23
        elif 150 <= ll < 160:
            alpha = 1.43670 ; beta = 0.67706  ; rr0 = 0.778 ; saa = 3  ; gamma = 0.42624
        elif 160 <= ll < 180:
            alpha = 6.84802 ; beta = -5.06864 ; rr0 = 0.676 ; saa = 50
        elif 180 <= ll < 190:
            alpha = 4.16321 ; beta = -5.80016 ; rr0 = 0.359 ; saa = 51 ; gamma = 0.60387
        elif 190 <= ll < 200:
            alpha = 0.78135 ; beta = -0.27826 ; rr0 = 1.152 ; saa = 4
        elif 200 <= ll < 210:
            alpha = 0.85535 ; beta = 0.20848  ; rr0 = 1.152 ; saa = 17
        elif 210 <= ll < 220:
            alpha = 0.52521 ; beta = 0.65726  ; rr0 = 1.152 ; saa = 7
        elif 220 <= ll < 230:
            alpha = 0.88376 ; beta = -0.44519 ; rr0 = 0.993 ; saa = 40 ; gamma = 0.06013
        elif 230 <= ll < 240:
            alpha = 0.42228 ; beta = 0.26304  ; rr0 = 0.803 ; saa = 26
        elif 240 <= ll < 250:
            alpha = 0.71318 ; beta = -0.67229 ; rr0 = 0.530 ; saa = 55 ; gamma = 0.20994
        elif 250 <= ll < 260:
            alpha = 0.99606 ; beta = -0.70103 ; rr0 = 0.710 ; saa = 48 ; gamma = 0.01323
        elif 260 <= ll < 270:
            alpha = 0.91519 ; beta = -0.39690 ; rr0 = 1.152 ; saa = 48 ; gamma = 0.01961
        elif 270 <= ll < 280:
            alpha = 0.85791 ; beta = -0.29115 ; rr0 = 1.152 ; saa = 19
        elif 280 <= ll < 290:
            alpha = 1.44226 ; beta = -1.09775 ; rr0 = 0.657 ; saa = 39
        elif 290 <= ll < 300:
            alpha = 2.55486 ; beta = -1.68293 ; rr0 = 0.759 ; saa = 31
        elif 300 <= ll < 310:
            alpha = 3.18047 ; beta = -2.69796 ; rr0 = 0.589 ; saa = 40
        elif 210 <= ll < 320:
            alpha = 2.11235 ; beta = -1.77506 ; rr0 = 0.595 ; saa = 29
        elif 320 <= ll < 330:
            alpha = 1.75884 ; beta = -1.38574 ; rr0 = 0.635 ; saa = 25
        elif 330 <= ll < 340:
            alpha = 1.97257 ; beta = -1.55545 ; rr0 = 0.634 ; saa = 34 ; gamma = 0.00043
        elif 340 <= ll < 350:
            alpha = 1.41497 ; beta = -1.05722 ; rr0 = 0.669 ; saa = 46 ; gamma = 0.03264
        elif 350 <= ll < 360:
            alpha = 1.17795 ; beta = -0.95012 ; rr0 = 0.620 ; saa = 46 ; gamma = 0.03339

    if -5 < bb < 5:
        gamma = 0
        if 0   <= ll < 10:
            alpha = 2.62556 ; beta = -1.11097 ; rr0 = 1.182 ; saa = 57 ; gamma = 0.00340
        elif 10  <= ll < 20:
            alpha = 3.14461 ; beta = -1.01140 ; rr0 = 1.555 ; saa = 42
        elif 20  <= ll < 30:
            alpha = 4.26624 ; beta = -1.61242 ; rr0 = 1.323 ; saa = 34
        elif 30  <= ll < 40:
            alpha = 2.54447 ; beta = -0.12771 ; rr0 = 1.300 ; saa = 30
        elif 40  <= ll < 50:
            alpha = 2.27030 ; beta = -0.68720 ; rr0 = 1.652 ; saa = 52 ; gamma = 0.04928
        elif 50  <= ll < 60:
            alpha = 1.34359 ; beta = -0.05416 ; rr0 = 2.000 ; saa = 32
        elif 60  <= ll < 70:
            alpha = 1.76327 ; beta = -0.26407 ; rr0 = 2.000 ; saa = 37
        elif 70  <= ll < 80:
            alpha = 2.20666 ; beta = -0.41651 ; rr0 = 2.000 ; saa = 36
        elif 80  <= ll < 90:
            alpha = 1.50130 ; beta = -0.09855 ; rr0 = 1.475 ; saa = 45
        elif 90  <= ll < 100:
            alpha = 2.43965 ; beta = -0.82128 ; rr0 = 1.485 ; saa = 36 ; gamma = 0.01959
        elif 100 <= ll < 110:
            alpha = 3.35775 ; beta = -1.16400 ; rr0 = 0.841 ; saa = 35 ; gamma = 0.00298
        elif 110 <= ll < 120:
            alpha = 2.60621 ; beta = -0.68687 ; rr0 = 1.897 ; saa = 36
        elif 120 <= ll < 130:
            alpha = 2.90112 ; beta = -0.97988 ; rr0 = 1.480 ; saa = 32
        elif 130 <= ll < 140:
            alpha = 2.55377 ; beta = -0.71214 ; rr0 = 1.793 ; saa = 38
        elif 140 <= ll < 150:
            alpha = 3.12598 ; beta = -1.21437 ; rr0 = 1.287 ; saa = 23 ; gamma = 0.15298
        elif 150 <= ll < 160:
            alpha = 3.66930 ; beta = -2.29731 ; rr0 = 0.799 ; saa = 40 ; gamma = 0.33473
        elif 160 <= ll < 170:
            alpha = 2.15465 ; beta = -0.70690 ; rr0 = 1.524 ; saa = 37 ; gamma = 0.14017
        elif 170 <= ll < 180:
            alpha = 1.82465 ; beta = -0.60223 ; rr0 = 1.515 ; saa = 29 ; gamma = 0.20730
        elif 180 <= ll < 190:
            alpha = 1.76269 ; beta = -0.35945 ; rr0 = 2.000 ; saa = 28 ; gamma = 0.08052
        elif 190 <= ll < 200:
            alpha = 1.06085 ; beta = -0.14211 ; rr0 = 2.000 ; saa = 48
        elif 200 <= ll < 210:
            alpha = 1.21333 ; beta = -0.23225 ; rr0 = 2.000 ; saa = 57
        elif 210 <= ll < 220:
            alpha = 0.58326 ; beta = -0.06097 ; rr0 = 2.000 ; saa = 30 ; gamma = 0.36962
        elif 220 <= ll < 230:
            alpha = 0.74200 ; beta = -0.19293 ; rr0 = 1.923 ; saa = 48 ; gamma = 0.07459
        elif 230 <= ll < 240:
            alpha = 0.67520 ; beta = -0.21041 ; rr0 = 1.604 ; saa = 49 ; gamma = 0.16602
        elif 240 <= ll < 250:
            alpha = 0.62609 ; beta = -0.25312 ; rr0 = 1.237 ; saa = 73 ; gamma = 0.14437
        elif 250 <= ll < 260:
            alpha = 0.61415 ; beta = -0.13788 ; rr0 = 2.000 ; saa = 42 ; gamma = 0.26859
        elif 260 <= ll < 270:
            alpha = 0.58108 ; beta = 0.01195  ; rr0 = 2.000 ; saa = 40 ; gamma = 0.07661
        elif 270 <= ll < 280:
            alpha = 0.68352 ; beta = -0.10743 ; rr0 = 2.000 ; saa = 50 ; gamma = 0.00849
        elif 280 <= ll < 290:
            alpha = 0.61747 ; beta = 0.02675  ; rr0 = 2,000 ; saa = 49
        elif 290 <= ll < 300:
            alpha = 0.06827 ; beta = -0.26290 ; rr0 = 2.000 ; saa = 44
        elif 300 <= ll < 310:
            alpha = 1.53631 ; beta = -0.36833 ; rr0 = 2.000 ; saa = 37 ; gamma = 0.02960
        elif 210 <= ll < 320:
            alpha = 1.94300 ; beta = -0.71445 ; rr0 = 1.360 ; saa = 36 ; gamma = 0.15643
        elif 320 <= ll < 330:
            alpha = 1.22185 ; beta = -0.40185 ; rr0 = 1.520 ; saa = 48 ; gamma = 0.07354
        elif 330 <= ll < 340:
            alpha = 1.79429 ; beta = -0.48657 ; rr0 = 1.844 ; saa = 43
        elif 340 <= ll < 350:
            alpha = 2.29545 ; beta = -0.84096 ; rr0 = 1.365 ; saa = 32
        elif 350 <= ll < 360:
            alpha = 2.07408 ; beta = -0.64745 ; rr0 = 1.602 ; saa = 36 ; gamma = 0.12750


    if 5 < bb < 15:
        gamma = 0
        if 0   <= ll < 10:
            alpha = 2.94213 ; beta = -2.09158 ; rr0 = 0.703 ; saa = 41 ; gamma = 0.05490
        elif 10  <= ll < 30:
            alpha = 3.04627 ; beta = 7.71159  ; rr0 = 0.355 ; saa = 45
        elif 30  <= ll < 40:
            alpha = 3.78033 ; beta = -3.91956 ; rr0 = 0.482 ; saa = 42
        elif 40  <= ll < 50:
            alpha = 2.18119 ; beta = -2.4050  ; rr0 = 0.453 ; saa = 27
        elif 50  <= ll < 60:
            alpha = 1.45372 ; beta = -0.49522 ; rr0 = 1.152 ; saa = 31
        elif 60  <= ll < 70:
            alpha = 1.05051 ; beta = -1.01704 ; rr0 = 0.516 ; saa = 2
        elif 70  <= ll < 80:
            alpha = 0.48416 ; beta = -0.27182 ; rr0 = 0.891 ; saa = 94 ; gamma = 0.08639
        elif 80  <= ll < 90:
            alpha = 0.61963 ; beta = 0.41697  ; rr0 = 1.152 ; saa = 35 ; gamma = 0.47171
        elif 90  <= ll < 100:
            alpha = 4.40348 ; beta = -2.95011 ; rr0 = 0.745 ; saa = 52
        elif 100 <= ll < 120:
            alpha = 2.50938 ; beta = -0.56541 ; rr0 = 1.152 ; saa = 27
        elif 120 <= ll < 130:
            alpha = 0.44180 ; beta = 1.58923  ; rr0 = 0.949 ; saa = 4
        elif 130 <= ll < 140:
            alpha = 3.96081 ; beta = -3.37374 ; rr0 = 0.587 ; saa = 40 ; gamma = 0.34109
        elif 140 <= ll < 160:
            alpha = 2.53335 ; beta = -0.40541 ; rr0 = 1.152 ; saa = 38
        elif 160 <= ll < 170:
            alpha = 2.03760 ; beta = -0.66136 ; rr0 = 1.152 ; saa = 23
        elif 170 <= ll < 200:
            alpha = 1.06946 ; beta = -0.87395 ; rr0 = 0.612 ; saa = 29 ; gamma = 0.29230
        elif 200 <= ll < 210:
            alpha = 0.86348 ; beta = -0.65870 ; rr0 = 0.655 ; saa = 79 ; gamma = 0.09089
        elif 210 <= ll < 230:
            alpha = 0.30117 ; beta = -0.16136 ; rr0 = 0.933 ; saa = 17 ; gamma = 0.07495
        elif 230 <= ll < 240:
            alpha = 0.75171 ; beta = -0.57143 ; rr0 = 0.658 ; saa = 12 ; gamma = 0.00534
        elif 240 <= ll < 250:
            alpha = 1.97427 ; beta = -2.02654 ; rr0 = 0.487 ; saa = 67
        elif 250 <= ll < 260:
            alpha = 1.25208 ; beta = -1.47763 ; rr0 = 0.424 ; saa = 19 ; gamma = 0.09089
        elif 260 <= ll < 270:
            alpha = 0.89448 ; beta = -0.43870 ; rr0 = 1.019 ; saa = 5
        elif 270 <= ll < 280:
            alpha = 0.81141 ; beta = -0.51001 ; rr0 = 0.795 ; saa = 27 ; gamma = 0.03505
        elif 280 <= ll < 290:
            alpha = 0.83781 ; beta = -0.44138 ; rr0 = 0.949 ; saa = 50 ; gamma = 0.02820
        elif 290 <= ll < 300:
            alpha = 1.10600 ; beta = -0.86263 ; rr0 = 0.641 ; saa = 28 ; gamma = 0.03402
        elif 300 <= ll < 310:
            alpha = 1.37040 ; beta = -1.02779 ; rr0 = 0.667 ; saa = 28 ; gamma = 0.05608
        elif 310 <= ll < 320:
            alpha = 1.77590 ; beta = -1.26951 ; rr0 = 0.699 ; saa = 37 ; gamma = 0.06972
        elif 320 <= ll < 330:
            alpha = 1.20865 ; beta = -0.70679 ; rr0 = 0.855 ; saa = 35 ; gamma = 0.02902
        elif 330 <= ll < 340:
            alpha = 2.28830 ; beta = -1.71890 ; rr0 = 0.666 ; saa = 42 ; gamma = 0.22887
        elif 340 <= ll < 350:
            alpha = 3.26278 ; beta = -0.94181 ; rr0 = 1.152 ; saa = 38
        elif 350 <= ll < 360:
            alpha = 2.58100 ; beta = -1.69237 ; rr0 = 0.763 ; saa = 53

    if 15 < bb < 30:
        gamma = 0
        if 0   <= ll < 20:
            alpha = 6.23279  ; beta = -10.30384 ; rr0 = 0.302 ; saa = 42
        elif 20  <= ll < 40:
            alpha = -4.47693 ; beta = -7.28366  ; rr0 = 0.307 ; saa = 29
        elif 40  <= ll < 60 :
            alpha =  1.22938 ; beta = -1.19030  ; rr0 = 0.516 ; saa = 5
        elif 60  <= ll < 80 :
            alpha =  0.84291 ; beta = -1.59338  ; rr0 = 0.265 ; saa = 4
        elif 80  <= ll < 100 :
            alpha =  0.23996 ; beta = 0.06304   ; rr0 = 0.523 ; saa = 32
        elif 100 <= ll < 140 :
            alpha =  0.40062 ; beta = -1.75628  ; rr0 = 0.114 ; saa = 16
        elif 140 <= ll < 180 :
            alpha =  0.56898 ; beta = -0.53331  ; rr0 = 0.523 ; saa = 41
        elif 180 <= ll < 200 :
            alpha = -0.95721 ; beta = 11.6917   ; rr0 = 0.240 ; saa = 2
        elif 200 <= ll < 220 :
            alpha = -0.19051 ; beta = 1.45670   ; rr0 = 0.376 ; saa = 1
        elif 220 <= ll < 240 :
            alpha =  2.31305 ; beta = -7.82531  ; rr0 = 0.148 ; saa = 95
        elif 240 <= ll < 260:
            alpha =  1.39169 ; beta = -1.72984  ; rr0 = 0.402 ; saa = 6
        elif 260 <= ll < 260:
            alpha =  1.59418 ; beta = -1.28296  ; rr0 = 0.523 ; saa = 36
        elif 280 <= ll < 300 :
            alpha =  1.57082 ; beta = -197295   ; rr0 = 0.398 ; saa = 10
        elif 300 <= ll < 320 :
            alpha =  1.95998 ; beta = -3.26159  ; rr0 = 0.300 ; saa = 11
        elif 320 <= ll < 340:
            alpha =  2.59567 ; beta = -4.84133  ; rr0 = 0.168 ; saa = 37
        elif 340 <= ll < 360:
            alpha =  5.30273 ; beta = -7.43033  ; rr0 = 0.357 ; saa = 37

    if 30 < bb < 45:
        gamma = 0
        if 0   <= ll < 20:
            alpha =  2.93960 ; beta  = -6.48019  ; rr0 = 0.227 ; saa = 77
        elif 20  <= ll < 50:
            alpha =  1.65864 ; beta  = -9.99317  ; rr0 = 0.083 ; saa = 99
        elif 50  <= ll < 80:
            alpha =  1.71831 ; beta  = -7.25286  ; rr0 = 0.118 ; saa = 28
        elif 80  <= ll < 110:
            alpha =  1.33617 ; beta  = -10.39799 ; rr0 = 0.064 ; saa = 99
        elif 110 <= ll < 160:
            alpha = -0.31330 ; beta  = 1.35622   ; rr0 = 0.329 ; saa = 24
        elif 160 <= ll < 190:
            alpha =  1.51984 ; beta  = -8.69502  ; rr0 = 0.087 ; saa = 99
        elif 190 <= ll < 220:
            alpha = -0.50758 ; beta  = 4.73320   ; rr0 = 0.250 ; saa = 78
        elif 220 <= ll < 250:
            alpha =  1.25864 ; beta  = -12.59627 ; rr0 = 0.050 ; saa = 70
        elif 250 <= ll < 280:
            alpha =  1.54243 ; beta  = -3.75065  ; rr0 = 0.205 ; saa = 10
        elif 280 <= ll < 320:
            alpha =  2.72258 ; beta  = -7.47806  ; rr0 = 0.182 ; saa = 5
        elif 320 <= ll < 340:
            alpha =  2.81435 ; beta  = -5.52139  ; rr0 = 0.255 ; saa = 10
        elif 340 <= ll < 360:
            alpha =  2.23818 ; beta  = 0.81772   ; rr0 = 0.329 ; saa = 19

    if 45 < bb < 60:
        gamma = 0
        if 0   <= ll < 60:
            alpha = 1.38587 ; beta  = -9.06536  ; rr0 = 0.076 ; saa = 3
        elif 60  <= ll < 90:
            alpha = 2.28570 ; beta  = -9.88812  ; rr0 = 0.116 ; saa = 3
        elif 90  <= ll < 110:
            alpha = 1.36385 ; beta  = -8.10127  ; rr0 = 0.084 ; saa = 4
        elif 110 <= ll < 170:
            alpha = 0.05943 ; beta  = -1.08126  ; rr0 = 0.027 ; saa = 50
        elif 170 <= ll < 200:
            alpha = 1.40171 ; beta  = -3.21783  ; rr0 = 0.218 ; saa = 99
        elif 200 <= ll < 230:
            alpha = 0.14718 ; beta  = 3.92670   ; rr0 = 0.252 ; saa = 14
        elif 230 <= ll < 290:
            alpha = 0.57124 ; beta  = -4.30242  ; rr0 = 0.066 ; saa = 10
        elif 290 <= ll < 330:
            alpha = 3.69891 ; beta  = -19.62204 ; rr0 = 0.094 ; saa = 5
        elif 330 <= ll < 360:
            alpha = 1.19563 ; beta  = -0.45043  ; rr0 = 0.252 ; saa = 9

    if 60 < bb < 90:
        gamma = 0
        if 0   <= ll < 30:
            alpha = 0.69443 ;  beta = -0.27600  ; rr0 = 0.153 ; saa = 99
        elif 30  <= ll < 60:
            alpha = 1.11811 ;  beta = 0.71179   ; rr0 = 0.085 ; saa = 73
        elif 60  <= ll < 90:
            alpha = 1.10427 ;  beta = -2.37654  ; rr0 = 0.123 ; saa = 99
        elif 90  <= ll < 120:
            alpha = -0.42211 ; beta = 5.24037   ; rr0 = 0.184 ; saa = 12
        elif 120 <= ll < 150:
            alpha = 0.87576 ;  beta = -4.38033  ; rr0 = 0.100 ; saa = 35
        elif 150 <= ll < 180:
            alpha = 1.27477 ;  beta = -4.98307  ; rr0 = 0.128 ; saa = 72
        elif 180 <= ll < 210:
            alpha = 1.19512 ;  beta = -6.58464  ; rr0 = 0.091 ; saa = 49
        elif 210 <= ll < 240:
            alpha = 0.97581 ;  beta = -4.89869  ; rr0 = 0.100 ; saa = 95
        elif 240 <= ll < 270:
            alpha = 0.54379 ;  beta = -0.84403  ; rr0 = 0.207 ; saa = 35
        elif 270 <= ll < 300:
            alpha = 0.85054 ;  beta = 13.01249  ; rr0 = 0.126 ; saa = 39
        elif 300 <= ll < 330:
            alpha = 0.74347 ;  beta = 1.39825   ; rr0 = 0.207 ; saa = 10
        elif 330 <= ll < 360:
            alpha = 0.77310 ;  beta = -4.45005  ; rr0 = 0.087 ; saa = 16

    return alpha, beta, gamma, rr0, saa



#******************************************************************************
#******************************************************************************

def stellar_class(photometry):

    # Intrinsic colors of dwarf and giant stars, for different spectral types.
    # From Bessell and Brett 1988

    JH_dwarfs = np.array([-0.05, 0.0, 0.02,0.06, 0.09, 0.13, 0.165, 0.23, 0.285, 0.305, 0.32, 0.33, 0.37, \
                 0.45, 0.5, 0.58, 0.61, 0.66, 0.695, 0.68, 0.665, 0.62, 0.60, 0.62, 0.66])
    HK_dwarfs = np.array([-0.035, 0.00, 0.005, 0.015, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.052, 0.055, \
                 0.06, 0.075, 0.09, 0.105, 0.11, 0.13, 0.165, 0.20, 0.21, 0.25, 0.275, 0.32, 0.37])

    JH_giants = np.array([0.37, 0.47, 0.50, 0.50, 0.54, 0.58, 0.63, 0.68, 0.73, 0.79, 0.83, 0.85, \
                 0.87, 0.90, 0.93, 0.95, 0.96, 0.96])
    HK_giants = np.array([0.065, 0.08, 0.085, 0.085, 0.095, 0.10, 0.115, 0.14, 0.15, 0.165, 0.19, 0.205, \
                 0.215, 0.235, 0.245, 0.285, 0.30, 0.31])

    # Cenvert the Bessell and Brett magnitudes to 2mass filters using the
    # relationships from Carpenter 2001

    JH_dwarfs_2mass = 0.990*JH_dwarfs - 0.049
    JH_giants_2mass = 0.990*JH_giants - 0.049
    HK_dwarfs_2mass = 1.00*HK_dwarfs + 0.034
    HK_giants_2mass = 1.00*HK_giants + 0.034

    line_dwarfs = geom.LineString(zip(HK_dwarfs_2mass, JH_dwarfs_2mass))
    line_giants = geom.LineString(zip(HK_giants_2mass, JH_giants_2mass))

    sp_class = 'dwarf'

    # Detect if JHK are present

    photo_keys = photometry.keys()
    if ('J' in photo_keys) and ('H' in photo_keys) and ('K' in photo_keys):
        J = photometry['J'][0] - photometry['J'][4]
        H = photometry['H'][0] - photometry['H'][4]
        K = photometry['K'][0] - photometry['K'][4]

        #print J,H,K

        JH = J-H
        HK = H-K

        # Compute distance from curves

        if HK > 0.11:
            point = geom.Point(HK, JH)
            d_dwarf = point.distance(line_dwarfs)
            d_giant = point.distance(line_giants)

            if d_giant < d_dwarf:
                sp_class = 'giant'

    return sp_class


def ini_logg(T, sp_type):

    if sp_type == 'dwarf':
        return 4.68*1e-8*T**2. - 8.33*1e-4*T + 7.547

    return max(1.5, -2.8*1e-7*T**2. + 3.79*1e-3*T - 9.335)


def ini_met(photometry):

    met = 0.0

    # From Martell & Laughlin 2002

    photo_keys = photometry.keys()
    if ('b' in photo_keys) and ('y' in photo_keys) and ('m1' in photo_keys) and ('c1' in photo_keys):
        b = photometry['b'][0] - photometry['b'][4]
        y = photometry['y'][0] - photometry['y'][4]
        by = b-y
        m1 = photometry['m1'][0]
        c1 = photometry['c1'][0]

        if (0.288 < by < 0.571) and (0.058 < m1 < 0.497) and (0.116 < c1 < 0.745):
            met = -10.424602 + 31.059003*by + 42.184476*m1 + 15.351995*c1 \
                  -11.239435*by**2. - 29.218135*m1**2. - 11.457610*c1**2. \
                  -138.92376*by*m1 - 52.033290*by*c1 + 11.259341*m1*c1 \
                  -46.087731*by**3. + 26.065099*m1**3. - 1.1017830*c1**3. \
                  +138.48588*by**2.*m1 + 39.012001*by**2.*c1 \
                  +23.225562*m1**2.*by - 69.146876*m1**2.*c1 \
                  +20.456093*c1**2.*by - 3.3302478*c1**2.*m1 \
                  +70.168761*by*m1*c1

            if (met > 0.5) or (met < -2.0):
                met = 0.0


    return met


#******************************************************************************
#******************************************************************************

def mamajek(photometry):
    # Define the spline representations

    tck_bv = (np.array([3000.,  3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                       3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                       3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                       5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                       5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                       6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                       7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                       8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                       14500.,  14500.,  14500.,  14500.]),\
              np.array([1.91      ,  1.89889695,  1.64625059,  1.69918967,  1.53801881,
                        1.55319352,  1.53106075,  1.48801053,  1.48838346,  1.47779774,
                        1.44527783,  1.4129107 ,  1.3861779 ,  1.30383319,  1.24749974,
                        1.12702398,  1.12023921,  0.98822216,  0.90145692,  0.84143165,
                        0.83585252,  0.74482766,  0.76955031,  0.70259105,  0.7107385 ,
                        0.68402323,  0.67252669,  0.657925  ,  0.65103668,  0.60442773,
                        0.5946758 ,  0.53880897,  0.53974386,  0.5054044 ,  0.48018914,
                        0.43561453,  0.42026151,  0.38487402,  0.36939307,  0.34017362,
                        0.28593417,  0.25906698,  0.24421401,  0.22166016,  0.17370659,
                        0.15529447,  0.14294398,  0.07675079,  0.07944651,  0.03617357,
                        -0.00639459, -0.04057817, -0.10613727, -0.10490602, -0.1217432 ,
                        -0.14      ,  0.        ,  0.        ,  0.        ,  0.]), 3)

    tck_ub = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                        3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                        3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                        5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                        5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                        6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                        7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                        8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                        14500.,  14500.,  14500.,  14500.]),\
              np.array([1.3       ,  1.20179034,  1.24160075,  1.23336   ,  1.17798138,
                        1.18406638,  1.17700818,  1.16693519,  1.17013625,  1.17072352,
                        1.17986289,  1.2021882 ,  1.21267162,  1.2295794 ,  1.2017036 ,
                        1.05346897,  1.03969546,  0.80481124,  0.6171521 ,  0.49146118,
                        0.4750023 ,  0.29442981,  0.33979346,  0.22523313,  0.23872882,
                        0.19140577,  0.17276234,  0.1465829 ,  0.13315183,  0.066202  ,
                        0.05692161, -0.0015156 ,  0.00900488, -0.01655514, -0.02194147,
                        -0.03008711, -0.03032406, -0.01427563, -0.00573432,  0.01545909,
                        0.06019658,  0.07873128,  0.08417486,  0.0893996 ,  0.09696149,
                        0.10135862,  0.09885859,  0.07296584,  0.06625457,  0.02453844,
                        0.00297067, -0.09506513, -0.32830524, -0.34885643, -0.4398635 ,
                        -0.504     ,  0.        ,  0.        ,  0.        ,  0.]), 3)
    tck_vrc = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                         3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                         3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                         5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                         5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                         6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                         7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                         8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                         14500.,  14500.,  14500.,  14500.]),\
               np.array([1.65600000e+00,   1.32770139e+00,   1.38299370e+00,
                         1.26214059e+00,   1.12637239e+00,   1.08196910e+00,
                         1.05922109e+00,   9.78020530e-01,   9.83105652e-01,
                         9.62037137e-01,   9.26381055e-01,   8.95803690e-01,
                         8.57189517e-01,   8.02404948e-01,   7.68251851e-01,
                         6.65413009e-01,   6.59434399e-01,   5.33274971e-01,
                         4.93438490e-01,   4.56602040e-01,   4.53096266e-01,
                         4.06592212e-01,   4.18801877e-01,   3.88421712e-01,
                         3.90584365e-01,   3.79183024e-01,   3.73163872e-01,
                         3.66587469e-01,   3.63795482e-01,   3.40100095e-01,
                         3.35535146e-01,   3.04888078e-01,   3.04710217e-01,
                         2.87939207e-01,   2.73458064e-01,   2.51675089e-01,
                         2.40855510e-01,   2.19520334e-01,   2.10105091e-01,
                         1.94043816e-01,   1.60928322e-01,   1.45505907e-01,
                         1.36444436e-01,   1.24121696e-01,   9.55016366e-02,
                         8.68505135e-02,   7.96215784e-02,   4.18827642e-02,
                         4.72293943e-02,   1.09376982e-02,  -5.68741127e-04,
                         -1.09590175e-02,  -4.90412980e-02,  -4.35706726e-02,
                         -5.70515816e-02,  -6.20000000e-02,   0.00000000e+00,
                         0.00000000e+00,   0.00000000e+00,   0.00000000e+00]), 3)
    tck_vic = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                         3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                         3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                         5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                         5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                         6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                         7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                         8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                         14500.,  14500.,  14500.,  14500.]),\
               np.array([3.66400000e+00,   3.06242198e+00,   3.14656116e+00,
                         2.88392932e+00,   2.55639210e+00,   2.42552295e+00,
                         2.36793462e+00,   2.09680131e+00,   2.10740788e+00,
                         2.02854987e+00,   1.89751198e+00,   1.78441302e+00,
                         1.65940613e+00,   1.52030218e+00,   1.43397845e+00,
                         1.23361981e+00,   1.22358199e+00,   1.00591487e+00,
                         9.39956011e-01,   8.77891597e-01,   8.69357489e-01,
                         7.93108035e-01,   8.12970737e-01,   7.57444973e-01,
                         7.63387713e-01,   7.41341381e-01,   7.31803241e-01,
                         7.19266831e-01,   7.14400179e-01,   6.72923293e-01,
                         6.62428793e-01,   6.06546956e-01,   6.08918038e-01,
                         5.74350574e-01,   5.49242862e-01,   5.04453730e-01,
                         4.85671214e-01,   4.44003309e-01,   4.27209183e-01,
                         3.91710468e-01,   3.29824669e-01,   2.98986092e-01,
                         2.80944030e-01,   2.55096915e-01,   2.00997866e-01,
                         1.80830665e-01,   1.67657117e-01,   9.19376557e-02,
                         1.01149638e-01,   2.84913731e-02,  -1.09392951e-03,
                         -2.01745535e-02,  -1.11433319e-01,  -9.72864498e-02,
                         -1.28847171e-01,  -1.45000000e-01,   0.00000000e+00,
                         0.00000000e+00,   0.00000000e+00,   0.00000000e+00]), 3)
    tck_vks = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                         3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                         3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                         5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                         5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                         6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                         7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                         8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                         14500.,  14500.,  14500.,  14500.]),\
               np.array([6.50000000e+00,   5.62931496e+00,   5.75367026e+00,
                         5.34076308e+00,   4.79071352e+00,   4.61391909e+00,
                         4.51720164e+00,   4.13159547e+00,   4.14709596e+00,
                         4.03097681e+00,   3.85578132e+00,   3.70546709e+00,
                         3.54424972e+00,   3.31949749e+00,   3.16721742e+00,
                         2.81656851e+00,   2.79358299e+00,   2.39537651e+00,
                         2.17786419e+00,   2.02134802e+00,   1.99884413e+00,
                         1.78559326e+00,   1.84401087e+00,   1.68769762e+00,
                         1.70675113e+00,   1.64431540e+00,   1.61765528e+00,
                         1.58260890e+00,   1.56619519e+00,   1.45865442e+00,
                         1.43726664e+00,   1.30978554e+00,   1.31255319e+00,
                         1.23323500e+00,   1.17676105e+00,   1.07373032e+00,
                         1.03727999e+00,   9.50782384e-01,   9.14314901e-01,
                         8.42414011e-01,   7.12893077e-01,   6.47709841e-01,
                         6.12165798e-01,   5.57520113e-01,   4.37245386e-01,
                         3.91233825e-01,   3.60535900e-01,   1.93847001e-01,
                         2.05967449e-01,   6.93649668e-02,   3.93919827e-02,
                         -5.51783950e-03,  -2.63670696e-01,  -2.23184229e-01,
                         -3.14197676e-01,  -3.58000000e-01,   0.00000000e+00,
                         0.00000000e+00,   0.00000000e+00,   0.00000000e+00]), 3)
    tck_jh = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                        3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                        3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                        5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                        5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                        6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                        7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                        8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                        14500.,  14500.,  14500.,  14500.]),\
              np.array([0.588     ,  0.57925679,  0.55974591,  0.55662775,  0.55739209,
                        0.58028578,  0.58376255,  0.60634672,  0.60546662,  0.6120368 ,
                        0.61989396,  0.62481195,  0.63594655,  0.61500676,  0.60674412,
                        0.55966082,  0.55230917,  0.49293409,  0.43609132,  0.40202698,
                        0.39806407,  0.34706382,  0.35960846,  0.32332711,  0.3280911 ,
                        0.31165243,  0.30654383,  0.29678003,  0.29401661,  0.26762904,
                        0.26081307,  0.22934639,  0.23111651,  0.20988225,  0.19754575,
                        0.17130757,  0.16340732,  0.1446157 ,  0.13839745,  0.12210648,
                        0.09356932,  0.08068036,  0.07127585,  0.06038954,  0.03801209,
                        0.02849096,  0.02386206, -0.00845365, -0.00723668, -0.02947868,
                        -0.0309659 , -0.04121424, -0.06064657, -0.06648605, -0.07545834,
                        -0.081     ,  0.        ,  0.        ,  0.        ,  0.]), 3)
    tck_hk = (np.array([3000.,   3000.,   3000.,   3000.,   3100.,   3200.,   3250.,
                        3410.,   3500.,   3550.,   3650.,   3700.,   3850.,   3880.,
                        3970.,   4050.,   4230.,   4410.,   4620.,   4840.,   5040.,
                        5170.,   5280.,   5340.,   5490.,   5530.,   5590.,   5660.,
                        5680.,   5720.,   5770.,   5880.,   5920.,   6040.,   6170.,
                        6240.,   6340.,   6510.,   6640.,   6720.,   6810.,   7030.,
                        7200.,   7440.,   7500.,   7800.,   8000.,   8080.,   8270.,
                        8550.,   8840.,   9200.,   9700.,  10400.,  10700.,  12500.,
                        14500.,  14500.,  14500.,  14500.]),\
              np.array([ 0.329     ,  0.30083187,  0.3061432 ,  0.28706635,  0.257594  ,
                         0.2529637 ,  0.2486753 ,  0.22822261,  0.22937567,  0.22361754,
                         0.20932367,  0.19499336,  0.1796351 ,  0.16107941,  0.14881621,
                         0.12776045,  0.12821934,  0.10896355,  0.10001716,  0.09330098,
                         0.09300045,  0.08391132,  0.08457407,  0.07918705,  0.07958247,
                         0.07690231,  0.0745582 ,  0.07391399,  0.07305143,  0.06821834,
                         0.06662477,  0.06210192,  0.06123406,  0.06016717,  0.05572764,
                         0.05437476,  0.0521449 ,  0.05093807,  0.04959402,  0.04716756,
                         0.04441532,  0.04422156,  0.04038873,  0.04118725,  0.03772925,
                         0.03817168,  0.03710862,  0.03284367,  0.03506379,  0.02983915,
                         0.02818454,  0.02373756,  0.00656967,  0.00374458, -0.00280642,
                         -0.007     ,  0.        ,  0.        ,  0.        ,  0.]), 3)
    tck_btvt = (np.array([3850.,   3850.,   3850.,   3850.,   3970.,   4050.,   4230.,
                          4410.,   4620.,   4840.,   5040.,   5170.,   5280.,   5340.,
                          5490.,   5530.,   5590.,   5660.,   5680.,   5720.,   5770.,
                          5880.,   5920.,   6040.,   6170.,   6240.,   6340.,   6510.,
                          6640.,   6720.,   6810.,   7030.,   7200.,   7440.,   7500.,
                          7800.,   8000.,   8080.,   8270.,   8550.,   8840.,   9200.,
                          9700.,  10400.,  10700.,  12500.,  14500.,  14500.,  14500.,
                          14500.]),\
                np.array([1.65      ,  1.62529982,  1.61107882,  1.49837621,  1.4387566 ,
                          1.36306363,  1.32576584,  1.16207044,  1.06271882,  0.97207866,
                          0.97009969,  0.85213091,  0.88670237,  0.79478911,  0.8051525 ,
                          0.76802269,  0.75579062,  0.73434453,  0.72579215,  0.66145907,
                          0.64843825,  0.567861  ,  0.56983969,  0.52834313,  0.4990927 ,
                          0.4526401 ,  0.43614294,  0.40085542,  0.38573171,  0.3529186 ,
                          0.31238673,  0.28288031,  0.26840011,  0.24412267,  0.19710763,
                          0.18078778,  0.17080152,  0.09332344,  0.09761852,  0.04697616,
                          0.00787136, -0.02889995, -0.10493612, -0.10123972, -0.12268253,
                          -0.142     ,  0.        ,  0.        ,  0.        ,  0.]), 3)

    T = np.array([])
    err_mag = np.array([])
    photo_keys = photometry.keys()
    mult_zeros = []
    colors = ['B-V', 'U-B', 'V-R', 'V-I', 'V-K', 'J-H', 'H-K', 'Bt-Vt']
    tcks = [tck_bv, tck_ub, tck_vrc, tck_vic, tck_vks, tck_jh, tck_hk, tck_btvt]

    for i in range(8):
        c1, c2 = colors[i].split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            C = C1-C2
            e_C = np.sqrt(photometry[c1][1]**2. + photometry[c2][1]**2.)

            tck_new = (tcks[i][0], tcks[i][1]-C, tcks[i][2])
            zeros = interpolate.sproot(tck_new)

            if len(zeros) == 1:
                T = np.append(T, zeros[0])
                err_mag = np.append(err_mag, e_C)
            if len(zeros) > 1:
                mult_zeros.append((i, colors[i], zeros, err_C))

    if len(T) > 0:

        T_average = np.average(T, weights = 1./err_mag)

        if len(mult_zeros) > 0:
            for i in range(len(mult_zeros)):
                d = np.abs(mult_zeros[i][2]-T_average)

                T = np.append(T, mult_zeros[i][2][np.argmin(d)])
                err_mag = np.append(err_mag, mult_zeros[i][3])

        T_average = np.average(T, weights = 1./err_mag)
	err_T_average = np.sqrt(np.average((T - T_average)**2., weights = 1./err_mag))
	if err_T_average == 0.0:
            err_T_average = 100.

    else:
        T_average = 0.0
	err_T_average = 0.0

    return T_average, err_T_average



def gonzalez_hernandez(photometry, met):
    colors = ['B-V', 'V-R', 'V-J', 'V-H', 'V-K', 'J-K']
    coefficients = [( 0.4967,  0.7260, -0.1563,  0.0255, -0.0585, -0.0061),\
                    ( 0.4530,  1.4347, -0.5883, -0.0156, -0.0096, -0.0039),\
                    ( 0.4629,  0.4124, -0.0417, -0.0012,  0.0094,  0.0013),\
                    ( 0.5321,  0.2649, -0.0146, -0.0069,  0.0211,  0.0009),\
                    ( 0.5293,  0.2489, -0.0119, -0.0042,  0.0135,  0.0010),\
                    ( 0.6517,  0.6312,  0.0168, -0.0381,  0.0256,  0.0013)]
    T = np.array([])
    err_mag = np.array([])
    photo_keys = photometry.keys()
    T_final = 0.0
    err_T_final = 999.0
    for i in range(len(colors)):
        c1, c2 = colors[i].split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            C = C1-C2
            e_C = np.sqrt(photometry[c1][1]**2. + photometry[c2][1]**2.)

            color = ufloat(C, e_C)

            b = coefficients[i]
            theta_eff = b[0] + b[1]*color + b[2]*color**2 + b[3]*color*met + b[4]*met + b[5]*met**2
            if theta_eff.n != 0.0 and theta_eff.s < err_T_final:
                T = 5040./theta_eff
                T_final = T.n
                err_T_final = T.s
                #T = np.append(T, 5040./theta_eff)
                #err_mag = np.append(err_mag, e_C)

    #if len(T) > 0:

    #    T_average = np.average(T, weights = 1./err_mag)
    #else:
    #    T_average = 0.0

    #return T_average
    return T_final, err_T_final


#******************************************************************************
#******************************************************************************



def coefs_mann(type_color, met = True):
    if met:
        if type_color == 'V-J':
            return [2.515, -1.054, 0.2965, -0.04150, 0.002245, 0.05262]
        elif type_color == 'V-I':
            return [1.901, -0.6564, 0.1471, -0.01274, 0.0, 0.04697]
        else:
            return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    else:
        if type_color == 'V-J':
            return [2.769, -1.421, 0.4284, -0.06133, 0.003310, 0.1333, 0.05416]
        elif type_color == 'V-I':
            return [1.568, -0.4381, 0.07749, -0.005610, 0.0, 0.2441, -0.09257]
        else:
            return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


def mann(photometry, met):
    photo_keys = photometry.keys()
    if met != None:
        if ('V' in photo_keys) and ('J' in photo_keys):
            i = coefs_mann('V-J', met = True)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['J'][0] - photometry['J'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['J'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        elif ('V' in photo_keys) and ('I' in photo_keys):
            i = coefs_mann('V-I', met = True)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['I'][0] - photometry['I'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['I'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        else:
            i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            c = 0.0
            color = ufloat(c, 999.0)

        T = i[0] + i[1]*color + i[2]*color**2 + i[3]*color**3 + i[4]*color**4 + i[5]*met

    else:
        if ('V' in photo_keys) and ('J' in photo_keys):
            i = coefs_mann('V-J', met = False)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['J'][0] - photometry['J'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['J'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        elif ('V' in photo_keys) and ('I' in photo_keys):
            i = coefs_mann('V-I', met = False)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['I'][0] - photometry['I'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['I'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        else:
            i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            c = 0.0
            color = ufloat(c, 999.0)

        #print i, c

        if ('J' in photo_keys) and ('H' in photo_keys):
            C1 = photometry['J'][0] - photometry['J'][4]
            C2 = photometry['H'][0] - photometry['H'][4]
            c2 = C1-C2
            err_c2 = np.sqrt(photometry['J'][1]**2. + photometry['H'][1]**2.)
            color2 = ufloat(c2, err_c2)

            T = i[0] + i[1]*color + i[2]*color**2 + i[3]*color**2 + i[4]*color**3 +\
                i[5]*color2 + i[6]*color2**2

        else:
            T = 0.0

    T = T*3500.
    return T.n, T.s


#******************************************************************************
#******************************************************************************



def coefs_casagrande(type_color, color):
    if type_color == 'B-V' and (0.18 <= color.n <= 1.29):
        return [0.5665, 0.4809, -0.0060, -0.0613, -0.0042, -0.0055]

    elif type_color == 'V-R' and (0.24 <= color.n <= 0.80):
        return [0.4386, 1.4614, -0.7014, -0.0807, 0.0142, -0.0015]

    elif type_color == 'R-I' and (0.23 <= color.n <= 0.68):
        return [0.3296, 1.9716, -1.0225, -0.0298, 0.0329, 0.0035]

    elif type_color == 'V-I' and (0.46 <= color.n <= 1.47):
        return [0.4033, 0.8171, -0.1987, -0.0409, 0.0319, 0.0012]

    elif type_color == 'V-J' and (0.61 <= color.n <= 2.44):
        return [0.4669, 0.3849, -0.0350, -0.0140, 0.0225, 0.0011]

    elif type_color == 'V-H' and (0.67 <= color.n <= 3.01):
        return [0.5251, 0.2553, -0.0119, -0.0187, 0.0410, 0.0025]

    elif type_color == 'V-K' and (0.78 <= color.n <= 3.15):
        return [0.5057, 0.2600, -0.0146, -0.0131, 0.0288, 0.0016]

    elif type_color == 'J-K' and (0.07 <= color.n <= 0.80):
        return [0.6393, 0.6104, 0.0920, -0.0330, 0.0291, 0.0020]

    elif type_color == 'Bt-Vt' and (0.19 <= color.n <= 1.49):
        return [0.5839, 0.4000, -0.0067, -0.0282, -0.0346, -0.0087]

    elif type_color == 'Vt-J' and (0.77 <= color.n <= 2.56):
        return [0.4525, 0.3797, -0.0357, -0.0082, 0.0123, -0.0009]

    elif type_color == 'Vt-H' and (0.77 <= color.n <= 3.16):
        return [0.5286, 0.2354, -0.0073, -0.0182, 0.0401, 0.0021]

    elif type_color == 'Vt-K' and (0.99 <= color.n <= 3.29):
        return [0.4892, 0.2634, -0.0165, -0.0121, 0.0249, -0.0001]

    elif type_color == 'b-y' and (0.18 <= color.n <= 0.72):
        return [0.5796, 0.4812, 0.5747, -0.0633, 0.0042, -0.0055]

    else:
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]



def casagrande(photometry, met):
    possible_colors = ['Vt-K', 'V-K', 'Vt-H', 'V-H', 'V-J', 'Vt-J', 'V-I', \
                       'b-y', 'V-R', 'B-V', 'Bt-Vt', 'R-I', 'J-K']
    photo_keys = photometry.keys()
    T_array = 0.0
    T_array2 = []
    color_p = 'no_color'
    T_final = 0.0
    err_T_final = 100000.
    color_p_final = 'no_color'
    for c in possible_colors:
        #print c
        c1, c2 = c.split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            C1 = ufloat(photometry[c1][0] - photometry[c1][4], photometry[c1][1])
            C2 = ufloat(photometry[c2][0] - photometry[c2][4], photometry[c2][1])
            color = C1-C2
            i = coefs_casagrande(c, color)
            t = i[0] + i[1]*color + i[2]*color**2. + i[3]*color*met\
                + i[4]*met + i[5]*met**2.
            #print i, c, color, t, (t < 5040./3500.), (t > 5040./10000.)
            if t.n!=0.0 and (t.n < 5040./3500.) and (t.n > 5040./10000.):
                T_array = 5040./t
                T_array2.append(T_array)
                #print T_array, c
                color_p = c
                if T_array.s < err_T_final:
                    T_final = T_array.n
                    err_T_final = T_array.s
                    color_p_final = c
                del t
                #break


    del possible_colors, photo_keys
    #print T_array, color_p
    #print np.median(T_array2)
    #print T_final, err_T_final, color_p_final

    #return T_array.n, color_p
    #return np.median(T_array2).n, color_p
    return T_final, err_T_final, color_p_final


#******************************************************************************
#******************************************************************************


def use_relation(photometry):

    photo_keys = photometry.keys()
    limit_colors = {'B-V': 1.340, 'Bt-Vt': 1.543, 'V-R': 0.826, 'V-I': 1.580,\
                    'V-K': 3.418, 'J-H': 0.557, 'H-K': 0.169}

    use_mann = False
    for c in limit_colors.keys():
        c1,c2 = c.split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            val = C1-C2
            if val >= limit_colors[c]:
                use_mann = True
                break

    if use_mann == True:
        return 'mann'
    else:
        return 'casagrande'


#******************************************************************************
#******************************************************************************


def check_relation(photometry, xmetal, exception):

    relation = use_relation(photometry)
    T_c = 0.0

    if relation == 'casagrande':
        T_c, err_T_c, color_c = casagrande(photometry, xmetal)

    if (T_c == 0.0) or (relation == 'mann'):
        if exception == 1:
            T_c, err_T_c = mann(photometry, xmetal)
        else:
            T_c, err_T_c = mann(photometry, met = None)

        #print T_c

        if T_c == 0.0 or T_c > 4000.:
            T_c, err_T_c, color_c = casagrande(photometry, xmetal)
            relation = 'casagrande'

        else:
            color_c = 'any'

    return T_c, err_T_c, color_c, relation


#******************************************************************************
#******************************************************************************


def correct_casagrande(Tc, color, inst):
    # harps, feros, hires, uves
    corrections = {'B-V': (47.3, 54.7, 141.8, 53.8),\
                   'V-R': (31.3, 35.5, 301.3, 0.4),\
                   'R-I': (22.1, 35.8, 209.2, 65.6),\
                   'V-I': (24.4, 12.5, 276.8, 19.2),\
                   'V-J': (1.0, -19.7, 72.3, 18.1),\
                   'V-H': (18.2, -15.8, 103.3, 48.1),\
                   'V-K': (31.3, 19.8, 91.4, 72.9),\
                   'J-K': (87.4, 102.5, 196.8, 195.6),\
                   'Bt-Vt': (28.6, 26.9, 64.3, 26.7),\
                   'Vt-J': (7.6, -8.1, 45.6, 8.3),\
                   'Vt-H': (28.7, 2.0, 72.4, 40.6),\
                   'Vt-K': (26.5, 22.2, 74.9, 62.5),\
                   'b-y': (22.9, 9.7, 75.0, 59.3)}

    x = 0
    if inst == 'harps':
        x = corrections[color][0]
    elif inst == 'feros':
        x = corrections[color][1]
    elif inst == 'hires':
        x = corrections[color][2]
    elif inst == 'uves':
        x = corrections[color][3]
    return Tc + x
