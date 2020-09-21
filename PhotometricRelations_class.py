from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import zip
from builtins import range
from dataclasses import dataclass, field
from typing import List
import warnings
import numpy as np
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
import shapely.geometry as geom
from scipy import interpolate
from uncertainties import ufloat, unumpy
import mwdust

warnings.simplefilter("ignore")

@dataclass
class Quantity:
    name: str
    value: float = np.nan
    error: float = np.nan
    source: str = 'unknown'

    def setval(self, v, e, s, condition=True):
        if hasattr(v, 'value'):
            if ~np.isnan(float(v.value)) and ~np.isnan(float(e.value)) and condition:
                self.value = v
                self.error = e
                self.source = s
        else:
            if ~np.isnan(float(v)) and ~np.isnan(float(e)) and condition:
                self.value = v
                self.error = e
                self.source = s

    @property
    def isvalid(self):
        return self.source != 'unknown'

    def totuple(self):
        return [self.value, self.error, self.source]

@dataclass
class Magnitude(Quantity):
    pfilter: str = 'unknown'
    Av: float = 0.0

    def setmag(self, v, e, s, f, condition=True):
        self.setval(v, max(e, 0.01), s, condition=condition)
        if ~np.isnan(self.value):
            self.pfilter = f

    def totuplemag(self):
        return self.totuple()+[self.pfilter, self.Av]

    def setAv(self, v):
        if ~np.isnan(self.value) and ~np.isnan(float(v)) and (v >= 0.0):
            self.Av = v


_magnitudes = ['V', 'B', 'Bt', 'Vt', 'J', 'H', 'K', 'b', 'y', 'm1', 'c1',\
               'G', 'W1', 'W2', 'W3', 'R', 'I', 'BP', 'RP', 'E_BPRP']

def _make_mag_list():
    return [Magnitude(s) for s in _magnitudes]
def _make_coord_list():
    return [Quantity('RA'), Quantity('DEC')]
def _make_proper_motion_list():
    return [Quantity('pmRA'), Quantity('pmDEC')]

@dataclass
class Star:
    name: str
    magnitudes: List[Magnitude] = field(default_factory=_make_mag_list)
    coordinates: List[Quantity] = field(default_factory=_make_coord_list)
    propermotion: List[Quantity] = field(default_factory=_make_proper_motion_list)
    parallax: Quantity = Quantity('Parallax')
    galactic_coords: List[float] = field(default_factory=list)

    def getmag(self, mag):
        return self.magnitudes[_magnitudes.index(mag)]

    def getmags(self, mags):
        return [self.getmag(m) for m in mags]

    @property
    def validmags(self):
        return [self.getmag(m).name for m in _magnitudes if self.getmag(m).isvalid]

    def change_coords(self, ra, dec, s):
        [RA, DEC] = self.coordinates
        RA.setval(ra, 0.0, s)
        DEC.setval(dec, 0.0, s)

    def change_propermotion(self, pmra, pmdec, err_ra, err_dec, s):
        if ~np.isnan(pmra) and ~np.isnan(pmdec) and ~np.isnan(err_ra) and ~np.isnan(err_dec):
            pmRA, pmDEC = self.propermotion
            pmRA.setval(pmra, err_ra, s)
            pmDEC.setval(pmdec, err_dec, s)

    def todict(self):
        d = {'name': self.name}
        if all([c.isvalid for c in self.coordinates]):
            d['RA'] = [self.coordinates[0].value, self.coordinates[0].source]
            d['DEC'] = [self.coordinates[1].value, self.coordinates[1].source]
        if all([c.isvalid for c in self.propermotion]):
            d['pmRA'] = self.propermotion[0].totuple()
            d['pmDEC'] = self.propermotion[1].totuple()
        if self.parallax.isvalid:
            d['parallax'] = self.parallax.totuple()
        valid = self.validmags
        for v in valid:
            d[v] = self.getmag(v).totuplemag()
        if self.galactic_coords:
            d['galactic_coords'] = self.galactic_coords
        return d



def search_catalog_name(starname, vizier_object, catlist, catname):
    def _look_r(starname, catname, catlist, rlist, restrict=False):
        for r in rlist:
            result = vizier_object.query_object(starname, catalog=[catname], radius=r*u.arcsec)
            if str(result) != 'Empty TableList':
                if restrict:
                    p = result[0]['Plx']
                    if any(p >= 2.0):
                        catlist.append(result)
                        break
                else:
                    catlist.append(result)
                    break
        return catlist

    if catname == "I/345/gaia2":
        catlist = _look_r(starname, catname, catlist, [15., 25., 50.], restrict=True)
    elif catname == "I/284/out":
        catlist = _look_r(starname, catname, catlist, [10., 20.])
    elif catname == "I/311/hip2":
        catlist = _look_r(starname, catname, catlist, [5., 10., 20., 30.])
    else:
        catlist = _look_r(starname, catname, catlist, [5., 10., 20., 30., 50.])
    return catlist


def search_catalog_coords(RA, DEC, vizier_object, catlist, catname):
    def _look_r(C, catname, catlist, rlist, restrict=False):
        for r in rlist:
            result = vizier_object.query_region(SkyCoord(ra=C[0], dec=C[1], unit=(u.deg, u.deg),
                                                         frame='icrs'),\
                                                width=r, catalog=[catname])
            if str(result) != 'Empty TableList':
                if restrict:
                    p = result[0]['Plx']
                    if any(p >= 2.0):
                        catlist.append(result)
                        break
                else:
                    catlist.append(result)
                    break
        return catlist

    if catname == "I/345/gaia2":
        catlist = _look_r([RA, DEC], catname, catlist, ['15s', '25s', '50s'], restrict=True)
    elif catname == "I/284/out":
        catlist = _look_r([RA, DEC], catname, catlist, ['10s', '20s'])
    elif catname == "I/311/hip2":
        catlist = _look_r([RA, DEC], catname, catlist, ['5s', '10s', '20s', '30s'])
    else:
        catlist = _look_r([RA, DEC], catname, catlist, ['5s', '10s', '20s', '30s', '50s'])
    return catlist


def retrieve_catalogs(starname, vizier_object, use_coords=False, RA=None, DEC=None):
    result = []
    catalogues = ["II/246/out", "I/311/hip2", "I/259/tyc2", "J/MNRAS/373/13/table1",
                  "J/ApJS/168/128/survey", "II/215/catalog", "V/130/gcs3",
                  "J/MNRAS/403/1949/ubvri", "I/337/tgas", "II/237/colors",
                  "I/345/gaia2", "II/336/apass9", "I/284/out", "II/328/allwise"]
    if use_coords:
        for c in catalogues:
            result = search_catalog_coords(RA, DEC, vizier_object, result, c)
    else:
        starname = starname.replace('A', '').replace('-', ' ').replace('_', ' ')
        for c in catalogues:
            result = search_catalog_name(starname, vizier_object, result, c)
    return result


def vizier_params(starname, use_coords=False, RA=None, DEC=None):
    star = Star(starname)
    V, B, Bt, Vt, J, H, K, b, y, m1, c1, G, W1, W2, W3, R, I, BP, RP, E_BPRP = star.getmags(_magnitudes)
    if use_coords:
        star.change_coords(RA, DEC, "from user or fits header")

    try:
        e_hpmag = 0.01
        v = Vizier(columns=["**"])
        result = retrieve_catalogs(starname, v, use_coords, RA, DEC)

        if result:
            if len(result[0][0].columns) < 5:
                v = Vizier(columns=["*"])
                result = retrieve_catalogs(starname, v, use_coords, RA, DEC)

            name_cats = [list(r.keys())[0] for r in result]

            # 2MASS
            if "II/246/out" in name_cats:
                r = result[name_cats.index("II/246/out")][0]

                if len(r['Jmag']) > 1:
                    nA = np.zeros(len(r['Qflg']))
                    for k, q in enumerate(r['Qflg']):
                        nA[k] = q.count('A')
                    maxA = np.argmax(nA)
                else:
                    maxA = 0
                    
                if r['Qflg'][maxA].count('U') == 0:
                    J.setmag(r['Jmag'][maxA], r['e_Jmag'][maxA], "II/246/out", '2MASS J',\
                             condition=r['Qflg'][maxA][0] != 'U')
                    H.setmag(r['Hmag'][maxA], r['e_Hmag'][maxA], "II/246/out", '2MASS H',\
                             condition=r['Qflg'][maxA][1] != 'U')
                    K.setmag(r['Kmag'][maxA], r['e_Kmag'][maxA], "II/246/out", '2MASS Ks',\
                             condition=r['Qflg'][maxA][2] != 'U')

                    star.change_coords(r['RAJ2000'][maxA], r['DEJ2000'][maxA], "II/246/out")

            # GAIA DR2
            if "I/345/gaia2" in name_cats:
                r = result[name_cats.index("I/345/gaia2")][0]
                i2 = np.where(r['Teff'] > 0.0)[0]
                if len(r['Plx'][i2]) >= 1:
                    if len(r['Plx'][i2]) > 1:
                        i0 = np.where(r['Plx'][i2] > 2.0)[0]
                        if len(i0) > 1:
                            iav0 = np.where(np.abs(r['RV'][i2][i0]) > 0.0)[0]
                            if len(iav0) > 1:
                                itemp = np.argmax(r['Teff'][i2][i0[iav0]])
                                i0 = int(i0[iav0[itemp]])
                            elif len(iav0) == 1:
                                i0 = int(i0[iav0])
                            else:
                                i0 = int(i0[0])
                        elif len(i0) == 1:
                            i0 = int(i0[0])
                        else:
                            i0 = 0
                    else:
                        i0 = 0
                    if r['Plx'][i2][i0] > 2.0:
                        star.parallax.setval(r['Plx'][i2][i0]*u.mas, r['e_Plx'][i2][i0]*u.mas,
                                             "I/345/gaia2")
                        star.change_coords(r['RA_ICRS'][i2][i0], r['DE_ICRS'][i2][i0],
                                           "I/345/gaia2")
                        G.setmag(r['Gmag'][i2][i0], r['e_Gmag'][i2][i0], "I/345/gaia2", "GAIA G")
                        G.setAv(r['AG'][i2][i0])
                        star.change_propermotion(r['pmRA'][i2][i0], r['pmDE'][i2][i0],
                                                 r['e_pmRA'][i2][i0], r['e_pmDE'][i2][i0],
                                                 "I/345/gaia2")
                        BP.setmag(r['BPmag'][i2][i0], r['e_BPmag'][i2][i0], "I/345/gaia2", "GAIA BP")
                        RP.setmag(r['RPmag'][i2][i0], r['e_RPmag'][i2][i0], "I/345/gaia2", "GAIA RP")
                        E_BPRP.setmag(r['E_BP-RP_'][i2][i0], 0.1, "I/345/gaia2", "GAIA E(BP-RP)")


                elif r['Plx'][0] > 2.0:
                    star.parallax.setval(r['Plx'][0]*u.mas, r['e_Plx'][0]*u.mas, "I/345/gaia2",
                                         condition=r['Plx'][0] > 0.0)
                    star.change_coords(r['RA_ICRS'][0], r['DE_ICRS'][0], "I/345/gaia2")
                    G.setmag(r['Gmag'][0], r['e_Gmag'][0], "I/345/gaia2", "GAIA G")
                    G.setAv(r['AG'][0])
                    star.change_propermotion(r['pmRA'][0], r['pmDE'][0],
                                             r['e_pmRA'][0], r['e_pmDE'][0], "I/345/gaia2")
                    BP.setmag(r['BPmag'][0], r['e_BPmag'][0], "I/345/gaia2", "GAIA BP")
                    RP.setmag(r['RPmag'][0], r['e_RPmag'][0], "I/345/gaia2", "GAIA RP")
                    E_BPRP.setmag(r['E_BP-RP_'][0], 0.1, "I/345/gaia2", "GAIA E(BP-RP)")


                # Correct for the systematic from Stassun & Torres 2018
                if star.parallax.isvalid:
                    plx_corr = ufloat(star.parallax.value.value,
                                      star.parallax.error.value)\
                               + ufloat(82, 33)*1e-3
                    star.parallax.setval(plx_corr.n*u.mas, plx_corr.s*u.mas, "I/345/gaia2")

            # GAIA
            elif "I/337/tgas" in name_cats:
                r = result[name_cats.index("I/337/tgas")][0]
                star.parallax.setval(r['Plx'][0]*u.mas, r['e_Plx'][0]*u.mas, "I/337/tgas",
                                     condition=r['Plx'][0] > 0.0)
                star.change_coords(r['RA_ICRS'][0], r['DE_ICRS'][0], "I/337/tgas")

            # HIPPARCOS
            if ("I/311/hip2" in name_cats) and (not star.parallax.isvalid):
                r = result[name_cats.index("I/311/hip2")][0]
                star.parallax.setval(r['Plx'][0]*u.mas, r['e_Plx'][0]*u.mas, "I/311/hip2",
                                     condition=r['Plx'][0] > 0.0)
                star.change_coords(r['RArad'][0], r['DErad'][0], "I/311/hip2")
                if all([c.isvalid is False for c in star.propermotion]):
                    star.change_propermotion(r['pmRA'][0], r['pmDE'][0],
                                             r['e_pmRA'][0], r['e_pmDE'][0], "I/311/hip2")
                e_hpmag = r['e_Hpmag'][0] if r['e_Hpmag'][0] != 0.0 else 0.01


            # USNO-B1.0 (Monet+ 2003) for proper motions
            if "I/284/out" in name_cats:
                r = result[name_cats.index("I/284/out")][0]
                if all([c.isvalid is False for c in star.propermotion]):
                    star.change_propermotion(r['pmRA'][0], r['pmDE'][0],
                                             r['e_pmRA'][0], r['e_pmDE'][0], "I/284/out")

            # WISE
            if "II/328/allwise" in name_cats:
                r = result[name_cats.index("II/328/allwise")][0]
                if len(r['qph']) > 1:
                    counts = np.zeros(len(r['qph']))
                    flag_vals = {'A':6, 'B':5, 'C':4, 'U':3, 'X':2, 'Z':1}
                    for q, q_val in enumerate(r['qph']):
                        counts[q] = sum([flag_vals[l]*q_val.count(l) for l in flag_vals])
                    i0 = np.argmax(counts)
                    del counts
                else:
                    i0 = 0
                if all([~np.isnan(float(r['e_%smag' % k][i0])) for k in ['W1', 'W2', 'W3']]):
                    W1.setmag(r['W1mag'][i0], r['e_W1mag'][i0], "II/328/allwise", 'WISE-1')
                    W2.setmag(r['W2mag'][i0], r['e_W2mag'][i0], "II/328/allwise", 'WISE-2')
                    W3.setmag(r['W3mag'][i0], r['e_W3mag'][i0], "II/328/allwise", 'WISE-3')


            # The Tycho-2 Catalogue (Hog+ 2000)
            if "I/259/tyc2" in name_cats:
                r = result[name_cats.index("I/259/tyc2")][0]
                Bt.setmag(r['BTmag'][0], r['e_BTmag'][0], "I/259/tyc2", 'HIPPARCOS BT')
                Vt.setmag(r['VTmag'][0], r['e_VTmag'][0], "I/259/tyc2", 'HIPPARCOS VT')
                star.change_coords(r['RA_ICRS_'][0], r['DE_ICRS_'][0], "I/259/tyc2")


            # Koen et al. 2010
            if "J/MNRAS/403/1949/ubvri" in name_cats:
                r = result[name_cats.index("J/MNRAS/403/1949/ubvri")][0]
                V.setmag(r['Vmag'][0], 0.01, "J/MNRAS/403/1949/ubvri", 'Landolt V')
                B.setmag(r['B-V'][0] + r['Vmag'][0], 0.01, "J/MNRAS/403/1949/ubvri", 'Landolt B')
                R.setmag(r['Vmag'][0] - r['V-Rc'][0], 0.01, "J/MNRAS/403/1949/ubvri", 'Landolt R')
                I.setmag(r['Vmag'][0] - r['V-Ic'][0], 0.01, "J/MNRAS/403/1949/ubvri", 'Landolt I')

            else:
                # Casagrande et al. 2006
                if "J/MNRAS/373/13/table1" in name_cats:
                    r = result[name_cats.index("J/MNRAS/373/13/table1")][0]
                    V.setmag(r['Vmag'][0], 0.01, "J/MNRAS/373/13/table1", 'Landolt V')
                    B.setmag(r['B-V'][0] + r['Vmag'][0], 0.01, "J/MNRAS/373/13/table1", 'Landolt B')
                    R.setmag(r['Vmag'][0] - r['V-Rc'][0], 0.01,
                             "J/MNRAS/373/13/table1", 'Landolt R')
                    I.setmag(r['Vmag'][0] - r['V-Rc'][0] - r['R-Ic'][0], 0.01,
                             "J/MNRAS/373/13/table1", 'Landolt I')
                    J.setmag(r['Jmag'][0], r['e_Jmag'][0], "J/MNRAS/373/13/table1", '2MASS J')
                    H.setmag(r['Hmag'][0], r['e_Hmag'][0], "J/MNRAS/373/13/table1", '2MASS H')
                    K.setmag(r['Ksmag'][0], r['e_Ksmag'][0], "J/MNRAS/373/13/table1", '2MASS Ks')
                    if not star.parallax.isvalid:
                        star.parallax.setval(r['Plx'][0]*u.mas, r['e_Plx'][0]*u.mas,
                                             "J/MNRAS/373/13/table1", condition=r['Plx'][0] > 0.0)

                # Beers et al. 2007
                elif "J/ApJS/168/128/survey" in name_cats:
                    r = result[name_cats.index("J/ApJS/168/128/survey")][0]
                    V.setmag(r['Vmag'][0], r['e_Vmag'][0], "J/ApJS/168/128/survey", 'Landolt V')
                    B.setmag(r['B-V'][0] + r['Vmag'][0],
                             np.sqrt(r['e_Vmag'][0]**2 + r['e_B-V'][0]**2),
                             "J/ApJS/168/128/survey", 'Landolt B')
                    R.setmag(r['Vmag'][0] - r['V-Rc'][0],
                             np.sqrt(r['e_Vmag'][0]**2. + r['e_V-Rc'][0]**2.),
                             "J/ApJS/168/128/survey", 'Landolt R')
                    I.setmag(r['Vmag'][0] - r['V-Ic'][0],
                             np.sqrt(r['e_Vmag'][0]**2. + r['e_V-Ic'][0]**2.),
                             "J/ApJS/168/128/survey", 'Landolt I')
                    J.setmag(r['Jmag'][0], r['e_Jmag'][0], "J/ApJS/168/128/survey", '2MASS J')
                    H.setmag(r['Hmag'][0], r['e_Hmag'][0], "J/ApJS/168/128/survey", '2MASS H')
                    K.setmag(r['Ksmag'][0], r['e_Ksmag'][0], "J/ApJS/168/128/survey", '2MASS Ks')


            # HAUCK
            if "II/215/catalog" in name_cats:
                r = result[name_cats.index("II/215/catalog")][0]
                if not V.isvalid:
                    V.setmag(r['Vmag'][0], r['e_Vmag'][0], "II/215/catalog", 'Landolt V')
                b.setmag(r['b-y'][0], r['e_b-y'][0], "II/215/catalog", 'Stromgren b')
                y.setmag(0., r['e_b-y'][0], "II/215/catalog", 'Stromgren y')
                m1.setmag(r['m1'][0], r['e_m1'][0], "II/215/catalog", 'f')
                c1.setmag(r['c1'][0], r['e_c1'][0], "II/215/catalog", 'f')

            else:
                # GENEVA
                if "V/130/gcs3" in name_cats:
                    r = result[name_cats.index("V/130/gcs3")][0]
                    if not V.isvalid:
                        V.setmag(r['Vmag'][0], 0.01, "V/130/gcs3", 'Landolt V')
                    b.setmag(r['b-y'][0], 0.01, "V/130/gcs3", 'Stromgren b')
                    y.setmag(0., e_hpmag, "V/130/gcs3", 'Stromgren y')
                    m1.setmag(r['m1'][0], 0.01, "V/130/gcs3", 'f')
                    c1.setmag(r['c1'][0], 0.01, "V/130/gcs3", 'f')

            if all([not M.isvalid for M in [R, I, V, B]]) or\
               any([M.error > 0.2 for M in [J, H, K]]) or\
               all([not M.isvalid for M in [R, I, B]]):

                # Ducati 2002
                if "II/237/colors" in name_cats:
                    r = result[name_cats.index("II/237/colors")][0]
                    V.setmag(r['Vmag'][0], e_hpmag, "II/237/colors", 'Landolt V')
                    B.setmag(r['B-V'][0] + r['Vmag'][0], 0.01, "II/237/colors", 'Landolt B')
                    R.setmag(r['R-V'][0] + r['Vmag'][0], 0.01, "II/237/colors", 'Landolt R')
                    I.setmag(r['I-V'][0] + r['Vmag'][0], 0.01, "II/237/colors", 'Landolt I')

                    if any([M.error > 0.2 for M in [J, H, K]]):
                        J.setmag(r['J-V'][0] + r['Vmag'][0], r['Jsig'][0]*1e-2,
                                 "II/237/colors", '2MASS J')
                        H.setmag(r['H-V'][0] + r['Vmag'][0], r['Hsig'][0]*1e-2,
                                 "II/237/colors", '2MASS H')
                        K.setmag(r['K-V'][0] + r['Vmag'][0], r['Ksig'][0]*1e-2,
                                 "II/237/colors", '2MASS Ks')

                elif "II/336/apass9" in name_cats:
                    r = result[name_cats.index("II/336/apass9")][0]
                    if Vt.isvalid:
                        V.setmag(r['Vmag'][0], r['e_Vmag'][0], "II/336/apass9", 'Landolt V',
                                 condition=np.abs(Vt.value - r['Vmag'][0]) < 0.25)
                    else:
                        V.setmag(r['Vmag'][0], r['e_Vmag'][0], "II/336/apass9", 'Landolt V')
                    if Bt.isvalid:
                        B.setmag(r['Bmag'][0], r['e_Bmag'][0], "II/336/apass9", 'Landolt B',
                                 condition=np.abs(Bt.value - r['e_Bmag'][0]) < 0.25)
                    else:
                        B.setmag(r['Bmag'][0], r['e_Bmag'][0], "II/336/apass9", 'Landolt B')

        del result

    except Exception as e:
        print('error', e)

    # Get galactic coordinates if possible
    if all([c.isvalid for c in star.coordinates]):
        ra, dec = star.coordinates
        l, b = get_galactic_coords(ra, dec)
        star.galactic_coords = [l, b]

    photometry = star.todict()
    photometry = correct_extinction(photometry)
    del star

    return photometry


#******************************************************************************
#******************************************************************************

def get_galactic_coords(ra, dec):
    c = SkyCoord(ra=ra.value*u.deg, dec=dec.value*u.deg, frame='icrs')
    l = c.galactic.l.degree
    b = c.galactic.b.degree
    del c
    return l, b

def get_Av(photometry):
    A_V = 0.0
    if 'galactic_coords' in photometry:
        l, b = photometry['galactic_coords']

        if 'parallax' in photometry:
            p = photometry['parallax'][0]
            d = p.to(u.pc, equivalencies=u.parallax()).value/1000. # in Kpc
            ext_map = mwdust.Combined15()
        else:
            d = 1.
            ext_map = mwdust.SFD()

        A_V = 2.742*ext_map(l, b, d)
        del ext_map

    return A_V

def correct_extinction(photometry):
    def _correct_wave(l, b, d, wave, ext_map):
        if wave in ('HIPPARCOS BT', 'HIPPARCOS VT'):
            ext_map._filter = 'Landolt B'
            A_B = ext_map(l, b, d)
            ext_map._filter = 'Landolt V'
            A_V = ext_map(l, b, d)

            # Transform A_wave to the Hipparcos system
            if wave == 'HIPPARCOS BT':
                A_wave = 1.090/0.850*(A_B-A_V) + A_V
            else:
                A_wave = A_V + 0.090/0.850*(A_B-A_V)

        else:
            ext_map._filter = wave
            A_wave = ext_map(l, b, d)

        return max(A_wave[0], 0.0)

    mag = ['B', 'V', 'R', 'I', 'J', 'H', 'K', 'b', 'y', 'Bt', 'Vt', 'W1', 'W2']

    if 'galactic_coords' in photometry:
        l, b = photometry['galactic_coords']
        d = 1.0
        if 'parallax' in photometry:
            p = photometry['parallax'][0]
            d = p.to(u.pc, equivalencies=u.parallax()).value/1000. # in Kpc
            ext_map = mwdust.Combined15()
        else:
            ext_map = mwdust.SFD()

        for m in mag:
            if m in photometry:
                wave = photometry[m][3]
                A_wave = _correct_wave(l, b, d, wave, ext_map)
                if not np.isnan(A_wave) and\
                    (np.isnan(photometry[m][4]) or photometry[m][4] == 0.0):
                    photometry[m][4] = A_wave
        if 'G' in photometry:
            ext_map._filter=None
            E_BV = ext_map(l, b, d)[0]
            if not np.isnan(E_BV):
                photometry['G'][4] = max(2.35*E_BV, 0.0)

        if all([m in photometry for m in ['BP', 'RP', 'E_BPRP', 'G']]):
            ext_map._filter=None
            E_BV = ext_map(l, b, d)[0]
            A0 = 3.1*E_BV
            G_BVRPm = photometry['BP'][0]-photometry['RP'][0]
            E_BPRP = photometry['E_BPRP'][0]
            G_BVRP0 = G_BVRPm-E_BPRP

            c_BP = [1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
            c_RP = [0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]
            c_G = [0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]

            def Ax(c1, c2, c3, c4, c5, c6, c7):
                return c1 + c2*G_BVRP0 + c3*G_BVRP0**2 + c4*G_BVRP0**3 + c5*A0 + c6*A0**2 + c7*G_BVRP0*A0

            A_BP = Ax(*c_BP)*A0
            A_RP = Ax(*c_RP)*A0
            A_G = Ax(*c_G)*A0
            photometry['BP'][4] = A_BP
            photometry['RP'][4] = A_RP
            #print(photometry['G'][4], A_G, E_BV*2.35)
            photometry['G'][4] = A_G

        del ext_map

    return photometry


#******************************************************************************
#******************************************************************************

def stellar_class(photometry):

    # Intrinsic colors of dwarf and giant stars, for different spectral types.
    # From Bessell and Brett 1988

    JH_dwarfs = np.array([-0.05, 0.0, 0.02, 0.06, 0.09, 0.13, 0.165, 0.23, 0.285, 0.305, 0.32,\
                          0.33, 0.37, 0.45, 0.5, 0.58, 0.61, 0.66, 0.695, 0.68, 0.665, 0.62,\
                          0.60, 0.62, 0.66])
    HK_dwarfs = np.array([-0.035, 0.00, 0.005, 0.015, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,\
                          0.052, 0.055, 0.06, 0.075, 0.09, 0.105, 0.11, 0.13, 0.165, 0.20, 0.21,\
                          0.25, 0.275, 0.32, 0.37])

    JH_giants = np.array([0.37, 0.47, 0.50, 0.50, 0.54, 0.58, 0.63, 0.68, 0.73, 0.79, 0.83, 0.85,\
                          0.87, 0.90, 0.93, 0.95, 0.96, 0.96])
    HK_giants = np.array([0.065, 0.08, 0.085, 0.085, 0.095, 0.10, 0.115, 0.14, 0.15, 0.165, 0.19,\
                          0.205, 0.215, 0.235, 0.245, 0.285, 0.30, 0.31])

    # Cenvert the Bessell and Brett magnitudes to 2mass filters using the
    # relationships from Carpenter 2001

    i_dwarfs = np.where(HK_dwarfs > 0.14)[0]
    i_giants = np.where(HK_giants > 0.14)[0]

    JH_dwarfs_2mass = 0.990*JH_dwarfs[i_dwarfs] - 0.049
    JH_giants_2mass = 0.990*JH_giants[i_giants] - 0.049
    HK_dwarfs_2mass = 1.00*HK_dwarfs[i_dwarfs] + 0.034
    HK_giants_2mass = 1.00*HK_giants[i_giants] + 0.034

    line_dwarfs = geom.LineString(list(zip(HK_dwarfs_2mass, JH_dwarfs_2mass)))
    line_giants = geom.LineString(list(zip(HK_giants_2mass, JH_giants_2mass)))

    sp_class = 'dwarf'

    # Detect if JHK are present
    if ('J' in photometry) and ('H' in photometry) and ('K' in photometry):
        J = photometry['J'][0] - photometry['J'][4]
        H = photometry['H'][0] - photometry['H'][4]
        K = photometry['K'][0] - photometry['K'][4]
        JH = J-H
        HK = H-K

        # Compute distance from curves
        if HK > 0.14:
            point = geom.Point(HK, JH)
            d_dwarf = point.distance(line_dwarfs)
            d_giant = point.distance(line_giants)
            del point

            if d_giant < d_dwarf:
                sp_class = 'giant'

    del line_dwarfs, line_giants
    return sp_class


def stellar_class_pm(photometry):
    # Use the proper motion of the star for classification, following
    # Collier  Cameron et al. 2007

    sp_class = 'dwarf'

    if all([k in photometry for k in ['J', 'H', 'pmRA', 'pmDEC']]):
        try:
            pmRA = photometry['pmRA'][0]/1000.
            pmDEC = photometry['pmDEC'][0]/1000.
            mu = np.sqrt(pmRA**2. + pmDEC**2.)
            J = photometry['J'][0] - photometry['J'][4]
            H = photometry['H'][0] - photometry['H'][4]

            RPM = J + 5.*np.log10(mu)

            if (-15. <= RPM <= 10.) and (0.0 <= (J-H) <= 1.0):
                ycurve = -141.25*(J-H)**4. + 473.18*(J-H)**3.\
                         -583.60*(J-H)**2. + 313.42*(J-H) - 58.

                if RPM < ycurve:
                    sp_class = 'giant'
            else:
                return stellar_class(photometry)
        except:
            return stellar_class(photometry)
    else:
        return stellar_class(photometry)

    return sp_class


def ini_logg(T, sp_type):
    if sp_type == 'dwarf':
        return 4.68*1e-8*T**2. - 8.33*1e-4*T + 7.547
    return max(1.0, -5.82*1e-7*T**2. + 6.73*1e-3*T - 10.65)


def ini_met(photometry):
    met = 0.0

    # From Martell & Laughlin 2002
    if 'b' in photometry and 'y' in photometry and 'm1' in photometry and 'c1' in photometry:
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

def correct_mamajek(Tc, color, inst):
    # harps, feros, hires, uves
    corrections = {'B-V': (-3.67, 54.8, -15.79, -6.12),\
                   'V-K': (15.23, 11.03, 43.25, 42.25),\
                   'J-H': (15.47, 32.62, 78.61, 59.06),\
                   'Bt-Vt': (64.61, 84.23, 108.15, 71.98)}

    x = 0
    if inst in ['harps', 'feros', 'hires', 'uves']:
        x = corrections[color][['harps', 'feros', 'hires', 'uves'].index(inst)]
    return Tc + x

def mamajek(photometry, inst):
    # Define the spline representations

    tck_bv = (np.array([3000., 3000., 3000., 3000., 3100., 3200., 3250.,
                        3410., 3500., 3550., 3650., 3700., 3850., 3880.,
                        3970., 4050., 4230., 4410., 4620., 4840., 5040.,
                        5170., 5280., 5340., 5490., 5530., 5590., 5660.,
                        5680., 5720., 5770., 5880., 5920., 6040., 6170.,
                        6240., 6340., 6510., 6640., 6720., 6810., 7030.,
                        7200., 7440., 7500., 7800., 8000., 8080., 8270.,
                        8550., 8840., 9200., 9700., 10400., 10700., 12500.,
                        14500., 14500., 14500., 14500.]),\
              np.array([1.91, 1.89889695, 1.64625059, 1.69918967, 1.53801881,
                        1.55319352, 1.53106075, 1.48801053, 1.48838346, 1.47779774,
                        1.44527783, 1.4129107, 1.3861779, 1.30383319, 1.24749974,
                        1.12702398, 1.12023921, 0.98822216, 0.90145692, 0.84143165,
                        0.83585252, 0.74482766, 0.76955031, 0.70259105, 0.7107385,
                        0.68402323, 0.67252669, 0.657925, 0.65103668, 0.60442773,
                        0.5946758, 0.53880897, 0.53974386, 0.5054044, 0.48018914,
                        0.43561453, 0.42026151, 0.38487402, 0.36939307, 0.34017362,
                        0.28593417, 0.25906698, 0.24421401, 0.22166016, 0.17370659,
                        0.15529447, 0.14294398, 0.07675079, 0.07944651, 0.03617357,
                        -0.00639459, -0.04057817, -0.10613727, -0.10490602, -0.1217432,
                        -0.14, 0., 0., 0., 0.]), 3)
    tck_vks = (np.array([3000., 3000., 3000., 3000., 3100., 3200., 3250.,
                         3410., 3500., 3550., 3650., 3700., 3850., 3880.,
                         3970., 4050., 4230., 4410., 4620., 4840., 5040.,
                         5170., 5280., 5340., 5490., 5530., 5590., 5660.,
                         5680., 5720., 5770., 5880., 5920., 6040., 6170.,
                         6240., 6340., 6510., 6640., 6720., 6810., 7030.,
                         7200., 7440., 7500., 7800., 8000., 8080., 8270.,
                         8550., 8840., 9200., 9700., 10400., 10700., 12500.,
                         14500., 14500., 14500., 14500.]),\
               np.array([6.50000000e+00, 5.62931496e+00, 5.75367026e+00,
                         5.34076308e+00, 4.79071352e+00, 4.61391909e+00,
                         4.51720164e+00, 4.13159547e+00, 4.14709596e+00,
                         4.03097681e+00, 3.85578132e+00, 3.70546709e+00,
                         3.54424972e+00, 3.31949749e+00, 3.16721742e+00,
                         2.81656851e+00, 2.79358299e+00, 2.39537651e+00,
                         2.17786419e+00, 2.02134802e+00, 1.99884413e+00,
                         1.78559326e+00, 1.84401087e+00, 1.68769762e+00,
                         1.70675113e+00, 1.64431540e+00, 1.61765528e+00,
                         1.58260890e+00, 1.56619519e+00, 1.45865442e+00,
                         1.43726664e+00, 1.30978554e+00, 1.31255319e+00,
                         1.23323500e+00, 1.17676105e+00, 1.07373032e+00,
                         1.03727999e+00, 9.50782384e-01, 9.14314901e-01,
                         8.42414011e-01, 7.12893077e-01, 6.47709841e-01,
                         6.12165798e-01, 5.57520113e-01, 4.37245386e-01,
                         3.91233825e-01, 3.60535900e-01, 1.93847001e-01,
                         2.05967449e-01, 6.93649668e-02, 3.93919827e-02,
                         -5.51783950e-03, -2.63670696e-01, -2.23184229e-01,
                         -3.14197676e-01, -3.58000000e-01, 0.00000000e+00,
                         0.00000000e+00, 0.00000000e+00, 0.00000000e+00]), 3)
    tck_jh = (np.array([3000., 3000., 3000., 3000., 3100., 3200., 3250.,
                        3410., 3500., 3550., 3650., 3700., 3850., 3880.,
                        3970., 4050., 4230., 4410., 4620., 4840., 5040.,
                        5170., 5280., 5340., 5490., 5530., 5590., 5660.,
                        5680., 5720., 5770., 5880., 5920., 6040., 6170.,
                        6240., 6340., 6510., 6640., 6720., 6810., 7030.,
                        7200., 7440., 7500., 7800., 8000., 8080., 8270.,
                        8550., 8840., 9200., 9700., 10400., 10700., 12500.,
                        14500., 14500., 14500., 14500.]),\
              np.array([0.588, 0.57925679, 0.55974591, 0.55662775, 0.55739209,
                        0.58028578, 0.58376255, 0.60634672, 0.60546662, 0.6120368,
                        0.61989396, 0.62481195, 0.63594655, 0.61500676, 0.60674412,
                        0.55966082, 0.55230917, 0.49293409, 0.43609132, 0.40202698,
                        0.39806407, 0.34706382, 0.35960846, 0.32332711, 0.3280911,
                        0.31165243, 0.30654383, 0.29678003, 0.29401661, 0.26762904,
                        0.26081307, 0.22934639, 0.23111651, 0.20988225, 0.19754575,
                        0.17130757, 0.16340732, 0.1446157, 0.13839745, 0.12210648,
                        0.09356932, 0.08068036, 0.07127585, 0.06038954, 0.03801209,
                        0.02849096, 0.02386206, -0.00845365, -0.00723668, -0.02947868,
                        -0.0309659, -0.04121424, -0.06064657, -0.06648605, -0.07545834,
                        -0.081, 0., 0., 0., 0.]), 3)
    tck_btvt = (np.array([3850., 3850., 3850., 3850., 3970., 4050., 4230.,
                          4410., 4620., 4840., 5040., 5170., 5280., 5340.,
                          5490., 5530., 5590., 5660., 5680., 5720., 5770.,
                          5880., 5920., 6040., 6170., 6240., 6340., 6510.,
                          6640., 6720., 6810., 7030., 7200., 7440., 7500.,
                          7800., 8000., 8080., 8270., 8550., 8840., 9200.,
                          9700., 10400., 10700., 12500., 14500., 14500., 14500.,
                          14500.]),\
                np.array([1.65, 1.62529982, 1.61107882, 1.49837621, 1.4387566,
                          1.36306363, 1.32576584, 1.16207044, 1.06271882, 0.97207866,
                          0.97009969, 0.85213091, 0.88670237, 0.79478911, 0.8051525,
                          0.76802269, 0.75579062, 0.73434453, 0.72579215, 0.66145907,
                          0.64843825, 0.567861, 0.56983969, 0.52834313, 0.4990927,
                          0.4526401, 0.43614294, 0.40085542, 0.38573171, 0.3529186,
                          0.31238673, 0.28288031, 0.26840011, 0.24412267, 0.19710763,
                          0.18078778, 0.17080152, 0.09332344, 0.09761852, 0.04697616,
                          0.00787136, -0.02889995, -0.10493612, -0.10123972, -0.12268253,
                          -0.142, 0., 0., 0., 0.]), 3)

    T = np.array([])
    err_mag = np.array([])
    mult_zeros = []
    colors = ['B-V', 'V-K', 'J-H', 'Bt-Vt']
    tcks = [tck_bv, tck_vks, tck_jh, tck_btvt]

    for i in range(4):
        c1, c2 = colors[i].split('-')
        if (c1 in photometry) and (c2 in photometry):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            C = C1-C2
            e_C = np.sqrt(photometry[c1][1]**2. + photometry[c2][1]**2.)

            tck_new = (tcks[i][0], tcks[i][1]-C, tcks[i][2])
            zeros = interpolate.sproot(tck_new)

            correction = correct_mamajek(0.0, colors[i], inst)
            zeros = zeros + correction

            if len(zeros) == 1:
                T = np.append(T, zeros[0])
                err_mag = np.append(err_mag, e_C)
            if len(zeros) > 1:
                mult_zeros.append((i, colors[i], zeros, e_C))

    if T.size > 0:
        T_average = np.average(T, weights=1./err_mag)

        if mult_zeros:
            for m in mult_zeros:
                d = np.abs(m[2]-T_average)
                T = np.append(T, m[2][np.argmin(d)])
                err_mag = np.append(err_mag, m[3])

        T_average = np.average(T, weights=1./err_mag)
        err_T_average = np.sqrt(np.average((T - T_average)**2., weights=1./err_mag))
        if err_T_average == 0.0:
            err_T_average = 100.

    else:
        T_average = 0.0
        err_T_average = 0.0

    return T_average, err_T_average


def alonso1999(photometry, met):
    colors = np.array(['V-K', 'V-K', 'J-H', 'J-K'])
    coefficients = [[0.5558, 0.2105, 0.001981, -0.009965, 0.01325, -0.002726],\
                    [0.3770, 0.3660, -0.03170, -0.003074, -0.002765, -0.002973],\
                    [0.5977, 1.015, -0.1020, -0.01029, 0.03006, 0.01013],\
                    [0.5816, 0.9134, -0.1443, 0.0000, 0.0000, 0.0000]]
    e_Teff = np.array([40., 25., 170., 125.])
    color_ranges = {'V-K':[{'FeH': [-0.5, 0.2], 'color': [0.20, 2.50], 'r': 0},
                           {'FeH': [-1.5, -0.5], 'color': [1.00, 2.50], 'r': 0},
                           {'FeH': [-2.5, -1.5], 'color': [1.20, 2.50], 'r': 0},
                           {'FeH': [-3.0, -2.5], 'color': [1.70, 2.50], 'r': 0},
                           {'FeH': [-0.5, 0.2], 'color': [2.00, 4.90], 'r': 1},
                           {'FeH': [-1.5, -0.5], 'color': [2.00, 4.60], 'r': 1},
                           {'FeH': [-2.5, -1.5], 'color': [2.00, 3.40], 'r': 1},
                           {'FeH': [-3.0, -2.5], 'color': [2.00, 2.80], 'r': 1}],
                    'J-H':[{'FeH': [-0.5, 0.2], 'color': [0.00, 0.90]},
                           {'FeH': [-1.5, -0.5], 'color': [0.20, 0.80]},
                           {'FeH': [-2.5, -1.5], 'color': [0.30, 0.70]},
                           {'FeH': [-3.0, -2.5], 'color': [0.35, 0.65]}],
                    'J-K':[{'FeH': [-0.5, 0.2], 'color': [0.00, 1.10]},
                           {'FeH': [-1.5, -0.5], 'color': [0.20, 1.00]},
                           {'FeH': [-2.5, -1.5], 'color': [0.30, 0.90]},
                           {'FeH': [-3.0, -2.5], 'color': [0.40, 0.80]}]}
    #corrections = np.array([0.0, 0.0, 0.0, 0.0])
    a_corr = np.array([2.70107158e-01, 8.86192858e-02, 2.57890150e-01])
    b_corr = np.array([3.58390280e+03, 4.48537163e+03, 3.72508078e+03])

    list_colors = ['V-K', 'J-H', 'J-K']
    T_array = np.ones(len(list_colors))*np.nan
    err_T_array = np.ones(len(list_colors))*np.nan
    T_final = 0.0
    err_T_final = 0.0

    for i, lc in enumerate(list_colors):
        c1, c2 = lc.split('-')
        if (c1 in photometry) and (c2 in photometry):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            C = C1-C2
            icolor = np.where(colors == lc)[0]
            coefs = None
            err_T_int = None

            # If there is only one relation for the color
            if len(icolor) == 1:
                # Check metallicity and color ranges
                in_ranges = 0
                for k in color_ranges[lc]:
                    range_m = k['FeH']
                    range_c = k['color']
                    if (range_m[0] <= met <= range_m[1]) and (range_c[0] <= C <= range_c[1]):
                        in_ranges = 1
                if in_ranges > 0:
                    coefs = coefficients[int(icolor)]
                    err_T_int = e_Teff[int(icolor)]

            # There are two equations for the color, depending on the color value
            else:
                in_ranges = 0
                r = -99
                for k in color_ranges[lc]:
                    range_m = k['FeH']
                    range_c = k['color']
                    range_r = k['r']
                    if (range_m[0] <= met <= range_m[1]) and (range_c[0] <= C <= range_c[1]):
                        in_ranges += 1
                        r = range_r
                if in_ranges == 1:
                    coefs = coefficients[icolor[r]]
                    err_T_int = e_Teff[icolor[r]]
                elif in_ranges > 1:
                    imin = np.argmin(e_Teff[icolor])
                    coefs = coefficients[icolor[imin]]
                    err_T_int = e_Teff[icolor[imin]]

            if coefs is not None and err_T_int is not None:
                theta = coefs[0] + coefs[1]*C + coefs[2]*C**2 + coefs[3]*C*met\
                        + coefs[4]*met + coefs[5]*met**2

                if theta != 0.0:
                    Teff = 5040./theta + ufloat(0.0, err_T_int)
                    T_array[i] = Teff.n
                    err_T_array[i] = Teff.s

    ii = np.where(~np.isnan(T_array))[0]
    if ii.size > 0:
        Tcorr = [np.polyval([a_corr[ii][i_], b_corr[ii][i_]], T_array[ii][i_])\
                 for i_ in range(len(ii))]
        Tmean = unumpy.uarray(Tcorr, err_T_array[ii])
        #T_final = np.average(T_array[ii], weights=1./err_T_array[ii])
        T_final = np.average(Tcorr, weights=1./err_T_array[ii])
        #T_final = np.median(T_array[ii])
        err_T_final = np.median(Tmean).s

        #Tmean = unumpy.uarray(T_array[ii]+corrections[ii], err_T_array[ii])
        #T_final = np.median(Tmean).n
        #err_T_final = np.median(Tmean).s

    return T_final, err_T_final


#******************************************************************************
#******************************************************************************

def correct_mann(Tc, color, inst, met):
    # harps, feros, hires, uves
    corrections = {'V-J': ((41.3, 26.3, 89.7, 53.4), (-87.8, -73.3, -87.8, -48.2)),\
                   'V-I': ((0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0))}

    m = 0 if met else 1
    x = 0
    if inst in ['harps', 'feros', 'hires', 'uves']:
        x = corrections[m][color][['harps', 'feros', 'hires', 'uves'].index(inst)]
    return Tc + x


def coefs_mann(type_color, met=True):
    if met:
        if type_color == 'V-J':
            return [2.515, -1.054, 0.2965, -0.04150, 0.002245, 0.05262]
        if type_color == 'V-I':
            return [1.901, -0.6564, 0.1471, -0.01274, 0.0, 0.04697]
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if type_color == 'V-J':
        return [2.769, -1.421, 0.4284, -0.06133, 0.003310, 0.1333, 0.05416]
    if type_color == 'V-I':
        return [1.568, -0.4381, 0.07749, -0.005610, 0.0, 0.2441, -0.09257]
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


def mann(photometry, met):
    if met is not None:
        if ('V' in photometry) and ('J' in photometry):
            i = coefs_mann('V-J', met=True)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['J'][0] - photometry['J'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['J'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        elif ('V' in photometry) and ('I' in photometry):
            i = coefs_mann('V-I', met=True)
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
        if ('V' in photometry) and ('J' in photometry):
            i = coefs_mann('V-J', met=False)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['J'][0] - photometry['J'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['J'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        elif ('V' in photometry) and ('I' in photometry):
            i = coefs_mann('V-I', met=False)
            C1 = photometry['V'][0] - photometry['V'][4]
            C2 = photometry['I'][0] - photometry['I'][4]
            c = C1-C2
            err_c = np.sqrt(photometry['I'][1]**2. + photometry['V'][1]**2.)
            color = ufloat(c, err_c)
        else:
            i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            c = 0.0
            color = ufloat(c, 999.0)


        if ('J' in photometry) and ('H' in photometry):
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
    coefs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if type_color == 'B-V' and (0.18 <= color.n <= 1.29):
        coefs = [0.5665, 0.4809, -0.0060, -0.0613, -0.0042, -0.0055]
    elif type_color == 'V-J' and (0.61 <= color.n <= 2.44):
        coefs = [0.4669, 0.3849, -0.0350, -0.0140, 0.0225, 0.0011]
    elif type_color == 'V-H' and (0.67 <= color.n <= 3.01):
        coefs = [0.5251, 0.2553, -0.0119, -0.0187, 0.0410, 0.0025]
    elif type_color == 'V-K' and (0.78 <= color.n <= 3.15):
        coefs = [0.5057, 0.2600, -0.0146, -0.0131, 0.0288, 0.0016]
    elif type_color == 'J-K' and (0.07 <= color.n <= 0.80):
        coefs = [0.6393, 0.6104, 0.0920, -0.0330, 0.0291, 0.0020]
    elif type_color == 'Bt-Vt' and (0.19 <= color.n <= 1.49):
        coefs = [0.5839, 0.4000, -0.0067, -0.0282, -0.0346, -0.0087]
    elif type_color == 'Vt-J' and (0.77 <= color.n <= 2.56):
        coefs = [0.4525, 0.3797, -0.0357, -0.0082, 0.0123, -0.0009]
    elif type_color == 'Vt-H' and (0.77 <= color.n <= 3.16):
        coefs = [0.5286, 0.2354, -0.0073, -0.0182, 0.0401, 0.0021]
    elif type_color == 'Vt-K' and (0.99 <= color.n <= 3.29):
        coefs = [0.4892, 0.2634, -0.0165, -0.0121, 0.0249, -0.0001]
    elif type_color == 'b-y' and (0.18 <= color.n <= 0.72):
        coefs = [0.5796, 0.4812, 0.5747, -0.0633, 0.0042, -0.0055]
    return coefs

def correct_casagrande(Tc, color, inst):
    # harps, feros, hires, uves
    corrections = {'Vt-K': [16.66, 25.67, 34.48, 31.58],
                   'V-K': [23.26, 26.81, 47.89, 52.3],
                   'Vt-H': [11.48, 18.96, 35.03, 45.29],
                   'V-H': [21.02, 23.92, 55.81, 53.8],
                   'V-J': [-7.26, -17.38, -10.37, 25.0],
                   'Vt-J': [-18.77, -4.76, -16.9, 13.75],
                   'b-y': [14.09, 24.5, 44.49, 29.27],
                   'B-V': [36.54, 65.18, 37.37, 19.43],
                   'Bt-Vt': [19.56, 19.39, 23.88, 11.53],
                   'J-K': [85.5, 124.89, 191.43, 120.83]}

    x = 0
    if inst in ['harps', 'feros', 'hires', 'uves']:
        x = corrections[color][['harps', 'feros', 'hires', 'uves'].index(inst)]
    return Tc + x

def casagrande(photometry, met, inst):
    possible_colors = ['Vt-K', 'V-K', 'Vt-H', 'V-H', 'V-J', 'Vt-J',\
                       'b-y', 'B-V', 'Bt-Vt', 'J-K']
    T_final = 0.0
    err_T_final = 100000.
    color_p_final = []
    T_array = []
    for c in possible_colors:
        c1, c2 = c.split('-')
        if (c1 in photometry) and (c2 in photometry):
            C1 = ufloat(photometry[c1][0] - photometry[c1][4], photometry[c1][1])
            C2 = ufloat(photometry[c2][0] - photometry[c2][4], photometry[c2][1])
            color = C1-C2
            i = coefs_casagrande(c, color)
            t = i[0] + i[1]*color + i[2]*color**2. + i[3]*color*met + i[4]*met + i[5]*met**2.
            if t.n != 0.0 and (t.n < 5040./3500.) and (t.n > 5040./10000.):
                T_array_v1 = (5040./t).n
                err_T_array_v1 = (5040./t).s
                T_corrected = correct_casagrande(T_array_v1, c, inst)
                T_array.append(ufloat(T_corrected, err_T_array_v1))
                color_p_final.append('%s: %.1f +- %.1f' % (c, T_corrected, err_T_array_v1))
                del t

    T_array = np.array(T_array)
    if T_array.size > 0:
        T_mean = np.mean(T_array)
        T_final = T_mean.n
        err_T_final = T_mean.s
        color_p_final = ", ".join(color_p_final)# color_p_final[:-2]

    del possible_colors, photometry
    return T_final, err_T_final, color_p_final


#******************************************************************************
#******************************************************************************


def use_relation(photometry):
    limit_colors = {'B-V': 1.340, 'Bt-Vt': 1.543, 'V-R': 0.826, 'V-I': 1.580,\
                    'V-K': 3.418, 'J-H': 0.557, 'H-K': 0.169}

    use_mann = False
    for c in list(limit_colors.keys()):
        c1, c2 = c.split('-')
        if (c1 in photometry) and (c2 in photometry):
            C1 = photometry[c1][0] - photometry[c1][4]
            C2 = photometry[c2][0] - photometry[c2][4]
            val = C1-C2
            if val >= limit_colors[c]:
                use_mann = True
                break

    if use_mann:
        return 'mann'
    return 'casagrande'


#******************************************************************************
#******************************************************************************


def check_relation(photometry, xmetal, exception, inst):
    relation = use_relation(photometry)
    T_c = 0.0

    if relation == 'casagrande':
        T_c, err_T_c, color_c = casagrande(photometry, xmetal, inst)

    if (T_c == 0.0) or (relation == 'mann'):
        if exception == 1:
            T_c, err_T_c = mann(photometry, xmetal)
        else:
            T_c, err_T_c = mann(photometry, met=None)

        if T_c == 0.0 or T_c > 4000.:
            T_c, err_T_c, color_c = casagrande(photometry, xmetal, inst)
            relation = 'casagrande'
        else:
            color_c = 'any'

    return T_c, err_T_c, color_c, relation
