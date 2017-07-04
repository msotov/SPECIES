import numpy as np
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord


def vizier_params(starname, use_coords = False, RA = None, DEC = None):
    photometry = {}
    photometry['name'] = starname

    try:
        e_hpmag = 0.01
        v=Vizier(columns=["**"])
        if use_coords:
            result = v.query_region(coord.SkyCoord(ra = RA, dec = DEC,\
                                        unit = (u.deg, u.deg),\
                                        frame = 'icrs'),\
                                        width = '10s',\
                                        catalog = ["II/246/out",\
                                                   "I/311/hip2",\
                                                   "I/259/tyc2",\
                                                   "J/MNRAS/373/13/table1",\
                                                   "J/ApJS/168/128/survey",\
                                                   "II/215/catalog",\
                                                   "V/130/gcs3",\
                                                   "J/MNRAS/403/1949/ubvri",\
                                                   "I/337/tgas"])

            if str(result) == 'Empty TableList':
                result = v.query_region(coord.SkyCoord(ra = RA, dec = DEC,\
                                            unit = (u.deg, u.deg),\
                                            frame = 'icrs'),\
                                            width = '20s',\
                                            catalog = ["II/246/out",\
                                                       "I/311/hip2",\
                                                       "I/259/tyc2",\
                                                       "J/MNRAS/373/13/table1",\
                                                       "J/ApJS/168/128/survey",\
                                                       "II/215/catalog",\
                                                       "V/130/gcs3",\
                                                       "J/MNRAS/403/1949/ubvri",\
                                                       "I/337/tgas"])


        else:

            result = v.query_object(starname.replace('A', '').replace('-', ' ').replace('_', ' '), \
                                        catalog = ["II/246/out",\
                                                   "I/311/hip2",\
                                                   "I/259/tyc2",\
                                                   "J/MNRAS/373/13/table1",\
                                                   "J/ApJS/168/128/survey",\
                                                   "II/215/catalog",\
                                                   "V/130/gcs3",\
                                                   "J/MNRAS/403/1949/ubvri",\
                                                   "I/337/tgas"],
                                        radius = 10.*u.arcsec)

            if str(result) == 'Empty TableList':
                result = v.query_object(starname.replace('A', '').replace('-', ' ').replace('_', ' '), \
                                            catalog = ["II/246/out",\
                                                       "I/311/hip2",\
                                                       "I/259/tyc2",\
                                                       "J/MNRAS/373/13/table1",\
                                                       "J/ApJS/168/128/survey",\
                                                       "II/215/catalog",\
                                                       "V/130/gcs3",\
                                                       "J/MNRAS/403/1949/ubvri",\
                                                       "I/337/tgas"],
                                            radius = 20.*u.arcsec)

        if str(result) != 'Empty TableList':
            name_cats = result.keys()

            # 2MASS
            if "II/246/out" in name_cats:
                i = name_cats.index("II/246/out")

                jmag = result[i]['Jmag'][0]
                kmag = result[i]['Kmag'][0]
                hmag = result[i]['Hmag'][0]
                e_jmag = result[i]['e_Jmag'][0]
                e_kmag = result[i]['e_Kmag'][0]
                e_hmag = result[i]['e_Hmag'][0]

                f_errors = [0.01, 0.01, 0.01]

                for j,e in enumerate([e_jmag, e_kmag, e_hmag]):
                    if np.isnan(float(e)) == False:
                        f_errors[j] = e
                if np.isnan(float(jmag)) == False:
                    photometry['J'] = (jmag, f_errors[0])
                if np.isnan(float(kmag)) == False:
                    photometry['K'] = (kmag, f_errors[1])
                if np.isnan(float(hmag)) == False:
                    photometry['H'] = (hmag, f_errors[2])

            # GAIA
            if "I/337/tgas" in name_cats:
                i = name_cats.index("I/337/tgas")
                plx = result[i]['Plx'][0]
                e_plx = result[i]['e_Plx'][0]
                photometry['parallax'] = (plx, e_plx)

            else:

                # HIPPARCOS
                if "I/311/hip2" in name_cats:
                    i = name_cats.index("I/311/hip2")
                    plx = result[i]['Plx'][0]
                    e_plx = result[i]['e_Plx'][0]
                    photometry['parallax'] = (plx, e_plx)
                    e_hpmag = result[i]['e_Hpmag'][0]

                    if e_hpmag == 0.0:
                        e_hpmag = 0.01

            # The Tycho-2 Catalogue (Hog+ 2000)
            if "I/259/tyc2" in name_cats:
                i = name_cats.index("I/259/tyc2")
                bmag = result[i]['BTmag'][0]
                e_bmag = result[i]['e_BTmag'][0]
                vmag = result[i]['VTmag'][0]
                e_vmag = result[i]['e_VTmag'][0]
                if e_bmag == 0.0: e_bmag = 0.01
                if e_vmag == 0.0: e_vmag = 0.01
                if np.isnan(float(bmag)) == False:
                    photometry['Bt'] = (bmag, e_bmag)
                if np.isnan(float(vmag)) == False:
                    photometry['Vt'] = (vmag, e_vmag)

            # Koen et al. 2010
            if "J/MNRAS/403/1949/ubvri" in name_cats:
                i = name_cats.index("J/MNRAS/403/1949/ubvri")
                vmag = result[i]['Vmag'][0]
                BV = result[i]['B-V'][0]
                UB = result[i]['U-B'][0]
                VR = result[i]['V-Rc'][0]
                VI = result[i]['V-Ic'][0]
                if e_hpmag == 0.0: e_hpmag = 0.01
                if np.isnan(float(vmag)) == False:
                    photometry['V'] = (vmag, e_hpmag)
                if np.isnan(float(BV + vmag)) == False:
                    photometry['B'] = (BV + vmag, e_hpmag)
                if np.isnan(float(vmag - VR)) == False:
                    photometry['R'] = (vmag - VR, e_hpmag)
                if np.isnan(float(vmag - VI)) == False:
                    photometry['I'] = (vmag - VI, e_hpmag)

            else:
                # Casagrande et al. 2006
                if "J/MNRAS/373/13/table1" in name_cats:
                    i = name_cats.index("J/MNRAS/373/13/table1")
                    vmag = result[i]['Vmag'][0]
                    BV = result[i]['B-V'][0]
                    VR = result[i]['V-Rc'][0]
                    RI = result[i]['R-Ic'][0]
                    J = result[i]['Jmag'][0]
                    e_J = result[i]['e_Jmag'][0]
                    H = result[i]['Hmag'][0]
                    e_H = result[i]['e_Hmag'][0]
                    K = result[i]['Ksmag'][0]
                    e_K = result[i]['e_Ksmag'][0]
                    if e_hpmag == 0.0: e_hpmag = 0.01
                    if e_J == 0.0: e_J = 0.01
                    if e_H == 0.0: e_H = 0.01
                    if e_K == 0.0: e_K = 0.01

                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = (vmag, e_hpmag)
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = (BV + vmag, e_hpmag)
                    if np.isnan(float(vmag - VR)) == False:
                        photometry['R'] = (vmag - VR, e_hpmag)
                    if np.isnan(float(vmag - VR - RI)) == False:
                        photometry['I'] = (vmag - VR - RI, e_hpmag)
                    if np.isnan(float(J)) == False:
                        photometry['J'] = (J, e_J)
                    if np.isnan(float(H)) == False:
                        photometry['H'] = (H, e_H)
                    if np.isnan(float(K)) == False:
                        photometry['K'] = (K, e_K)

                # Beers et al. 2007
                elif "J/ApJS/168/128/survey" in name_cats:
                    i = name_cats.index("J/ApJS/168/128/survey")
                    vmag = result[i]['Vmag'][0]
                    e_vmag = result[i]['e_Vmag'][0]
                    BV = result[i]['B-V'][0]
                    e_BV = result[i]['e_B-V'][0]
                    VR = result[i]['V-Rc'][0]
                    e_VR = result[i]['e_V-Rc'][0]
                    VI = result[i]['V-Ic'][0]
                    e_VI = result[i]['e_V-Ic'][0]
                    J = result[i]['Jmag'][0]
                    e_J = result[i]['e_Jmag'][0]
                    H = result[i]['Hmag'][0]
                    e_H = result[i]['e_Hmag'][0]
                    K = result[i]['Ksmag'][0]
                    e_K = result[i]['e_Ksmag'][0]
                    if e_vmag == 0.0: e_vmag = 0.01
                    if e_J == 0.0: e_J = 0.01
                    if e_H == 0.0: e_H = 0.01
                    if e_K == 0.0: e_k = 0.01
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = (vmag, e_vmag)
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = (BV + vmag, np.sqrt(e_vmag**2. + e_BV**2.))
                    if np.isnan(float(vmag - VR)) == False:
                        photometry['R'] = (vmag - VR, np.sqrt(e_vmag**2. + e_VR**2.))
                    if np.isnan(float(vmag - VI)) == False:
                        photometry['I'] = (vmag - VI, np.sqrt(e_vmag**2. + e_VI**2.))
                    if np.isnan(float(J)) == False:
                        photometry['J'] = (J, e_J)
                    if np.isnan(float(H)) == False:
                        photometry['H'] = (H, e_H)
                    if np.isnan(float(K)) == False:
                        photometry['K'] = (K, e_K)


            # HAUCK
            if "II/215/catalog" in name_cats:
                i = name_cats.index("II/215/catalog")
                vmag = result[i]['Vmag'][0]
                e_vmag = result[i]['e_Vmag'][0]
                by = result[i]['b-y'][0]
                e_by = result[i]['e_b-y'][0]
                if e_vmag == 0.0: e_vmag = 0.01
                if e_by == 0.0: e_by = 0.01
                if ('V' in photometry.keys()) == False:
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = (vmag, e_vmag)
                if np.isnan(float(by)) == False:
                    photometry['b'] = (by, e_by)
                    photometry['y'] = (0., e_by)


            else:
                # GENEVA
                if "V/130/gcs3" in name_cats:
                    i = name_cats.index("V/130/gcs3")
                    vmag = result[i]['Vmag'][0]
                    by = result[i]['b-y'][0]
                    if e_hpmag == 0.0: e_hpmag = 0.01
                    if ('V' in photometry.keys()) == False:
                        if np.isnan(float(vmag)) == False:
                            photometry['V'] = (vmag, e_hpmag)
                    if np.isnan(float(by)) == False:
                        photometry['b'] = (by, e_hpmag)
                        photometry['y'] = (0., e_hpmag)

        del result

        if ('R' in photometry.keys()) == False and ('I' in photometry.keys()) == False and \
                ('V' in photometry.keys()) == False and ('B' in photometry.keys()) == False:

            if use_coords:
                results_UBV = v.query_region(coord.SkyCoord(ra = RA, dec = DEC,\
                                             unit = (u.deg, u.deg),\
                                             frame = 'icrs'),\
                                             width = '20s',\
                                             catalog = ["II/237/colors"])

                if str(results_UBV) != 'Empty TableList':
                    if e_hpmag == 0.0:
                        e_hpmag = 0.01
                    vmag = results_UBV[0]['Vmag'][0]
                    UV = results_UBV[0]['U-V'][0]
                    BV = results_UBV[0]['B-V'][0]
                    RV = results_UBV[0]['R-V'][0]
                    IV = results_UBV[0]['I-V'][0]
                    if np.isnan(float(UV + vmag)) == False:
                        photometry['U'] = (UV + vmag, e_hpmag)
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = (vmag, e_hpmag)
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = (BV + vmag, e_hpmag)
                    if np.isnan(float(RV + vmag)) == False:
                        photometry['R'] = (RV + vmag, e_hpmag)
                    if np.isnan(float(IV + vmag)) == False:
                        photometry['I'] = (IV + vmag, e_hpmag)



                else:
                    results_UBV = v.query_region(coord.SkyCoord(ra = RA, dec = DEC,\
                                                 unit = (u.deg, u.deg),\
                                                 frame = 'icrs'),\
                                                 width = '20s',\
                                                 catalog = ["II/336/apass9"])

                    if str(results_UBV) != 'Empty TableList':
                        if e_hpmag == 0.0:
                            e_hpmag = 0.01
                        vmag = results_UBV[0]['Vmag'][0]
                        e_vmag = results_UBV[0]['e_Vmag'][0]
                        bmag = results_UBV[0]['Bmag'][0]
                        e_bmag = results_UBV[0]['e_Bmag'][0]
                        if e_vmag == 0.0:
                            e_vmag = 0.01
                        if e_bmag == 0.0:
                            e_bmag = 0.01
                        if np.isnan(float(vmag)) == False:
                            photometry['V'] = (vmag, e_vmag)
                        if np.isnan(float(bmag)) == False:
                            photometry['B'] = (bmag, e_bmag)


            # Ducati, 2002
            else:
                results_UBV = v.query_object(starname.replace('A', '').replace('-', ' ').replace('_', ' '),
                                                catalog=["II/237/colors"],\
                                                radius = 20.*u.arcsec)
                if str(results_UBV) != 'Empty TableList':
                    if e_hpmag == 0.0:
                        e_hpmag = 0.01
                    vmag = results_UBV[0]['Vmag'][0]
                    UV = results_UBV[0]['U-V'][0]
                    BV = results_UBV[0]['B-V'][0]
                    RV = results_UBV[0]['R-V'][0]
                    IV = results_UBV[0]['I-V'][0]
                    if np.isnan(float(UV + vmag)) == False:
                        photometry['U'] = (UV + vmag, e_hpmag)
                    if np.isnan(float(vmag)) == False:
                        photometry['V'] = (vmag, e_hpmag)
                    if np.isnan(float(BV + vmag)) == False:
                        photometry['B'] = (BV + vmag, e_hpmag)
                    if np.isnan(float(RV + vmag)) == False:
                        photometry['R'] = (RV + vmag, e_hpmag)
                    if np.isnan(float(IV + vmag)) == False:
                        photometry['I'] = (IV + vmag, e_hpmag)


                else:
                    results_UBV = v.query_object(starname.replace('A', '').replace('-', ' ').replace('_', ' '),
                                                    catalog=["II/336/apass9"],\
                                                    radius = 20.*u.arcsec)
                    if str(results_UBV) != 'Empty TableList':
                        if e_hpmag == 0.0:
                            e_hpmag = 0.01
                        vmag = results_UBV[0]['Vmag'][0]
                        e_vmag = results_UBV[0]['e_Vmag'][0]
                        bmag = results_UBV[0]['Bmag'][0]
                        e_bmag = results_UBV[0]['e_Bmag'][0]
                        if e_vmag == 0.0:
                            e_vmag = 0.01
                        if e_bmag == 0.0:
                            e_bmag = 0.01
                        if np.isnan(float(vmag)) == False:
                            photometry['V'] = (vmag, e_vmag)
                        if np.isnan(float(bmag)) == False:
                            photometry['B'] = (bmag, e_bmag)


            del results_UBV

        #check that are the magnitudes are valid
        photo_keys = photometry.keys()
        for mag in photo_keys:
            try:
                if np.isnan(float(photometry[mag][0])) or np.isnan(float(photometry[mag][1])):
                    p = dict(photometry)
                    del p[mag]
                    photometry = p
            except ValueError:
                pass

    except:
        pass

    return photometry


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

#******************************************************************************
#******************************************************************************


def mann(photometry, met):
    photo_keys = photometry.keys()
    if met != None:
        if ('V' in photo_keys) and ('J' in photo_keys):
            i = coefs_mann('V-J', met = True)
            c = photometry['V'][0] - photometry['J'][0]
        elif ('V' in photo_keys) and ('I' in photo_keys):
            i = coefs_mann('V-I', met = True)
            c = photometry['V'][0] - photometry['I'][0]
        else:
            i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            c = 0.0

        T = i[0] + i[1]*c + i[2]*c**2 + i[3]*c**3 + i[4]*c**4 + i[5]*met

    else:
        if ('V' in photo_keys) and ('J' in photo_keys):
            i = coefs_mann('V-J', met = False)
            c = photometry['V'][0] - photometry['J'][0]
        elif ('V' in photo_keys) and ('I' in photo_keys):
            i = coefs_mann('V-I', met = False)
            c = photometry['V'][0] - photometry['I'][0]
        else:
            i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            c = 0.0

        if ('J' in photo_keys) and ('H' in photo_keys):
            c2 = photometry['J'][0] - photometry['H'][0]

            T = i[0] + i[1]*c + i[2]*c**2 + i[3]*c**2 + i[4]*c**3 +\
                i[5]*c2 + i[6]*c2**2

        else:

            T = 0.0

    return T*3500.


#******************************************************************************
#******************************************************************************



def coefs_casagrande(type_color, color):
    if type_color == 'B-V' and (0.18 <= color <= 1.29):
        return [0.5665, 0.4809, -0.0060, -0.0613, -0.0042, -0.0055]

    elif type_color == 'V-R' and (0.24 <= color <= 0.80):
        return [0.4386, 1.4614, -0.7014, -0.0807, 0.0142, -0.0015]

    elif type_color == 'R-I' and (0.23 <= color <= 0.68):
        return [0.3296, 1.9716, -1.0225, -0.0298, 0.0329, 0.0035]

    elif type_color == 'V-I' and (0.46 <= color <= 1.47):
        return [0.4033, 0.8171, -0.1987, -0.0409, 0.0319, 0.0012]

    elif type_color == 'V-J' and (0.61 <= color <= 2.44):
        return [0.4669, 0.3849, -0.0350, -0.0140, 0.0225, 0.0011]

    elif type_color == 'V-H' and (0.67 <= color <= 3.01):
        return [0.5251, 0.2553, -0.0119, -0.0187, 0.0410, 0.0025]

    elif type_color == 'V-K' and (0.78 <= color <= 3.15):
        return [0.5057, 0.2600, -0.0146, -0.0131, 0.0288, 0.0016]

    elif type_color == 'J-K' and (0.07 <= color <= 0.80):
        return [0.6393, 0.6104, 0.0920, -0.0330, 0.0291, 0.0020]

    elif type_color == 'Bt-Vt' and (0.19 <= color <= 1.49):
        return [0.5839, 0.4000, -0.0067, -0.0282, -0.0346, -0.0087]

    elif type_color == 'Vt-J' and (0.77 <= color <= 2.56):
        return [0.4525, 0.3797, -0.0357, -0.0082, 0.0123, -0.0009]

    elif type_color == 'Vt-H' and (0.77 <= color <= 3.16):
        return [0.5286, 0.2354, -0.0073, -0.0182, 0.0401, 0.0021]

    elif type_color == 'Vt-K' and (0.99 <= color <= 3.29):
        return [0.4892, 0.2634, -0.0165, -0.0121, 0.0249, -0.0001]

    elif type_color == 'b-y' and (0.18 <= color <= 0.72):
        return [0.5796, 0.4812, 0.5747, -0.0633, 0.0042, -0.0055]

    else:
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


#******************************************************************************
#******************************************************************************


def casagrande(photometry, met):
    #possible_colors = ['V-J', 'Vt-J', 'b-y', 'V-H', 'Vt-H', 'Bt-Vt', 'V-K',\
    #                   'Vt-K', 'J-K', 'B-V', 'V-R', 'R-I', 'V-I']
    possible_colors = ['Vt-K', 'V-K', 'Vt-H', 'V-H', 'V-J', 'Vt-J', 'V-I', \
                       'b-y', 'V-R', 'B-V', 'Bt-Vt', 'R-I', 'J-K']
    photo_keys = photometry.keys()
    T_array = 0.0
    color_p = 'no_color'
    for c in possible_colors:
        c1, c2 = c.split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            color = photometry[c1][0] - photometry[c2][0]
            i = coefs_casagrande(c, color)
            t = i[0] + i[1]*color + i[2]*color**2. + i[3]*color*met\
                + i[4]*met + i[5]*met**2.
            if t!=0.0 and (t < 5040./3500.) and (t > 5040./10000.):
                T_array = 5040./t
                #print T_array, c
                color_p = c
                del t
                break


    del possible_colors, photo_keys

    return T_array, color_p


#******************************************************************************
#******************************************************************************


def use_relation(photometry):

    '''

           B-V    Bt-Vt  U-B     V-Rc   V-Ic    V-Ks    J-H     H-K     Ks-W1
    K7V    1.340  1.543  1.222   0.826  1.580   3.418   0.624   0.169   0.062
    K8V    1.380  1.592  1.213   0.857  1.664   3.544   0.631   0.180   0.069
    K9V    1.420  1.634  1.197   0.902  1.808   3.737   0.624   0.198   0.082
    M0V    1.431  1.650  1.190   0.913  1.848   3.790   0.622   0.203   0.086
    M1V    1.480  ...    1.171   0.968  2.051   4.065   0.610   0.225   0.103
    M1.5V  1.486  ...    1.170   0.978  2.089   4.12    0.607   0.228   0.105
    M2V    1.500  ...    1.170   1.001  2.173   4.24    0.600   0.234   0.110
    M2.5V  1.522  ...    1.175   1.041  2.306   4.43    0.589   0.244   0.117
    M3V    1.544  ...    1.181   1.079  2.420   4.60    0.579   0.252   0.122
    M3.5V  1.602  ...    1.200   1.178  2.680   5.00    0.558   0.269   0.132
    M4V    1.661  ...    1.222   1.241  2.831   5.25    0.557   0.282   0.139
    M4.5V  1.72   ...    1.23    1.345  3.073   5.64    0.564   0.301   ...
    M5V    1.874  ...    1.24    1.446  3.277   5.94    0.580   0.311   ...
    M5.5V  1.91   ...    1.3     1.656  3.664   6.50    0.588   0.329   ...
    M6V    2.00   ...    1.3     1.950  4.100   7.30    0.605   0.352   ...
    M6.5V  2.06   ...    ...     2.003  4.284   7.60    0.609   0.364   ...
    M7V    2.06   ...    ...     2.180  4.520   8.05    0.613   0.386   ...
    M7.5V  2.17   ...    ...     2.160  4.560   8.45    0.650   0.422   ...
    M8V    2.20   ...    ...     2.150  4.600   8.73    0.677   0.447   ...
    M9V    ...    ...    ...     1.890  4.370   8.85    0.749   0.481   ...
    '''


    photo_keys = photometry.keys()
    limit_colors = {'B-V': 1.340, 'Bt-Vt': 1.543, 'V-R': 0.826, 'V-I': 1.580,\
                    'V-K': 3.418, 'J-H': 0.557, 'H-K': 0.169}

    use_mann = False
    for c in limit_colors.keys():
        c1,c2 = c.split('-')
        if (c1 in photo_keys) and (c2 in photo_keys):
            val = photometry[c1][0] - photometry[c2][0]
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

    if relation == 'casagrande':
        T_c, color_c = casagrande(photometry, xmetal)

    else:
        if exception == 1:
            T_c = mann(photometry, xmetal)
        else:
            T_c = mann(photometry, met = None)

        if T_c == 0.0 or T_c > 4000.:
            T_c, color_c = casagrande(photometry, xmetal)
            relation = 'casagrande'

        else:
            color_c = 'any'

    return T_c, color_c, relation


#******************************************************************************
#******************************************************************************


def correct_casagrande(Tc, color, inst):
    x = 0
    if color == 'B-V':
        if inst == 'harps':
            x = 47.3
        elif inst == 'feros':
            x = 54.7
        elif inst == 'hires':
            x = 141.8
        elif inst == 'uves':
            x = 53.8

    elif color == 'V-R':
        if inst == 'harps':
            x = 31.3
        elif inst == 'feros':
            x = 35.5
        elif inst == 'hires':
            x = 301.3
        elif inst == 'uves':
            x = 0.4

    elif color == 'R-I':
        if inst == 'harps':
            x = 22.1
        elif inst == 'feros':
            x = 35.8
        elif inst == 'hires':
            x = 209.2
        elif inst == 'uves':
            x = 65.6

    elif color == 'V-I':
        if inst == 'harps':
            x = 24.4
        elif inst == 'feros':
            x = 12.5
        elif inst == 'hires':
            x = 276.8
        elif inst == 'uves':
            x = 19.2

    elif color == 'V-J':
        if inst == 'harps':
            x = 1.0
        elif inst == 'feros':
            x = -19.7
        elif inst == 'hires':
            x = 72.3
        elif inst == 'uves':
            x = 18.1

    elif color == 'V-H':
        if inst == 'harps':
            x = 18.2
        elif inst == 'feros':
            x = -15.8
        elif inst == 'hires':
            x = 103.3
        elif inst == 'uves':
            x = 48.1

    elif color == 'V-K':
        if inst == 'harps':
            x = 31.3
        elif inst == 'feros':
            x = 19.8
        elif inst == 'hires':
            x = 91.4
        elif inst == 'uves':
            x = 72.9

    elif color == 'J-K':
        if inst == 'harps':
            x = 87.4
        elif inst == 'feros':
            x = 102.5
        elif inst == 'hires':
            x = 196.8
        elif inst == 'uves':
            x = 195.6

    elif color == 'Vt-Bt':
        if inst == 'harps':
            x = 28.6
        elif inst == 'feros':
            x = 26.9
        elif inst == 'hires':
            x = 64.3
        elif inst == 'uves':
            x = 26.7

    elif color == 'Vt-J':
        if inst == 'harps':
            x = 7.6
        elif inst == 'feros':
            x = -8.1
        elif inst == 'hires':
            x = 45.6
        elif inst == 'uves':
            x = 8.3

    elif color == 'Vt-H':
        if inst == 'harps':
            x = 28.7
        elif inst == 'feros':
            x = 2.0
        elif inst == 'hires':
            x = 72.4
        elif inst == 'uves':
            x = 40.6

    elif color == 'Vt-K':
        if inst == 'harps':
            x = 26.5
        elif inst == 'feros':
            x = 22.2
        elif inst == 'hires':
            x = 74.9
        elif inst == 'uves':
            x = 62.5

    else:
        if inst == 'harps':
            x = 22.9
        elif inst == 'feros':
            x = 9.7
        elif inst =='hires':
            x = 75.0
        elif inst == 'uves':
            x = 59.3

    return Tc + x
