'''
Last modified in 09/05/2017
'''

import numpy as np
import glob, os, re, sys
from scipy.interpolate import UnivariateSpline


def split_file(file_m, ametal, path):
    i = 0
    f = open(file_m)
    for line in f:
        line = line.strip()
        m1 = re.search(r'TEFF   (\d*\.)  GRAVITY (\d*\.\d*) LTE.*', line)
        m2 = re.search(r'EFF   (\d*\.)  GRAVITY (\d*\.\d*) LTE.*', line)
        if m1 or m2:
            if i > 0:
                new_file.close()
            if m1:
                Teff = float(m1.group(1))
                logg = float(m1.group(2))
            else:
                Teff = float(m2.group(1))
                logg = float(m2.group(2))

            new_file = open(os.path.join(path, ametal + '_' + str(int(Teff)) + '_' + str(logg) + '.dat'), 'w')
            i += 1
        try:
            new_file.writelines(line + '\n')
        except UnboundLocalError:
            pass

    try:
        new_file.close()
    except UnboundLocalError:
        pass
    f.close()


def read_file(f, MT, N):
    mod = open(f)
    tabp = np.zeros((MT, N))
    i = 0
    for line in mod:
        m = re.search('\A[0-9]\.[0-9]*',line)
        if m:
            columnas = line.split()
            for j in range(N):
                try:
                    tabp[i][j] = float(columnas[j])
                except IndexError:
                    break
            i += 1

            del columnas


    mod.close()
    return tabp



def interpol(starname, xteff, xlogg, xmetal, micro):
    T = np.array([3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,\
                  5750.,6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,\
                  8000.,8250.,8250.,8500.,8750.,9000.,9250.,9500.,9750.,\
                  10000.,10250.,10500.,10750.,11000.,11250.,11500.,11750.,\
                  12000.,12250.,12500.,12750.,13000.,14000.,15000.])
    G = np.linspace(0.,5.,11)
    XM = np.array([-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3,\
                   -0.2,-0.1,0.,0.1,0.2,0.3,0.5,1.0])

    M = 72
    N = 9

    xlogg = min(xlogg, 5.0)
    xlogg = max(xlogg, 0.0)
    micro = min(micro, 10.0)
    micro = max(micro, 0.0)
    xteff = max(xteff, 3500)
    xteff = min(xteff, 15000)
    xmetal = max(xmetal, -3.0)
    xmetal = min(xmetal, 1.0)

    tabA = np.zeros((M,N))
    tabB = np.zeros((M,N))
    tabC = np.zeros((M,N))
    tabD = np.zeros((M,N))

    tab1 = np.zeros((M,N))
    tab2 = np.zeros((M,N))
    tab3 = np.zeros((M,N))
    tab4 = np.zeros((M,N))

    tab = np.zeros((M,N))

    TabY = np.zeros(4)

    #----------------------------------------------
    it = np.where(T <= xteff)[0]
    it = int(it[-1])

    ig = np.where(G <= xlogg)[0]
    ig = int(ig[-1])

    imt = np.where(XM < xmetal)[0]
    imt = int(imt[-1])
    #----------------------------------------------

    if xlogg == 5.:
        ig = ig - 1
        xfrac = 1.
    else:
        xfrac = (xlogg-G[ig])/(G[ig+1]-G[ig])

    if xmetal == 1.:
        imt = imt - 1
        yfrac = 1.
    else:
        yfrac = (xmetal-XM[imt])/(XM[imt+1]-XM[imt])

    #----------------------------------------------

    if xteff < 3500.:
        print 'No Kurucz model for less than 3500 K.'
    elif (8751. <= xteff <= 8999.) == True:
        print 'No interpolation possible between 8750 and 9000 K.'
    elif xteff >= 9500.:
        print 'Interpolation possible but not produced by the program.'
    else:
        TabX = np.zeros(4)

        if (3750. <= xteff < 8500) == True:
            TabX[0] = T[it - 1]
            TabX[1] = T[it]
            TabX[2] = T[it + 1]
            TabX[3] = T[it + 2]
        elif xteff < 3750.:
            TabX[0] = 3500.
            TabX[1] = 3750.
            TabX[2] = 4000.
            TabX[3] = 4250.
            T[it - 1] = TabX[0]
            T[it] = TabX[1]
            T[it + 1] = TabX[2]
            T[it + 2] = TabX[3]
        elif (8500. <= xteff <= 8750.) == True:
            TabX[0] = 8000.
            TabX[1] = 8250.
            TabX[2] = 8500.
            TabX[3] = 8750.
            T[it - 1] = TabX[0]
            T[it] = TabX[1]
            T[it + 1] = TabX[2]
            T[it + 2] = TabX[3]
        elif (9000. <= xteff < 9250.) == True:
            TabX[0] = 9000.
            TabX[1] = 9250.
            TabX[2] = 9500.
            TabX[3] = 9750.
            T[it - 1] = TabX[0]
            T[it] = TabX[1]
            T[it + 1] = TabX[2]
            T[it + 2] = TabX[3]
        elif xteff >= 9250.:
            TabX[0] = T[it - 1]
            TabX[1] = T[it]
            TabX[2] = T[it + 1]
            TabX[3] = T[it + 2]

        nbreak = 0

        if xteff <= 8750.:
            MT = 72
            N = 9
        else:
            MT = 64
            N = 7

        for k in range(4):
            atemp1 = str(int(T[it - 1]))
            atemp2 = str(int(T[it]))
            atemp3 = str(int(T[it + 1]))
            atemp4 = str(int(T[it + 2]))

            if k == 1:
                ig = ig + 1
            elif k == 2:
                imt = imt + 1
                ig = ig - 1
            elif k == 3:
                ig = ig + 1


            if XM[imt] < 0.:
                signo = 'm'
            else:
                signo = 'p'
            avar = str(int(abs(XM[imt]*10)))
            if abs(XM[imt]*10) <= 5.:
                ametal = signo + '0' + avar
            else:
                ametal = signo + avar



            path = './atm_models/atlas9/grids/grid' + ametal
            path_odfnew = path + 'odfnew'

            if os.path.isdir(path_odfnew):
                if len(glob.glob(path_odfnew + '/' + ametal + '*.dat')) == 0:
                    files_m = glob.glob(path_odfnew + '/a' + ametal + '*odfnew.dat')[0]
                    split_file(files_m, ametal, path_odfnew)
                temps = glob.glob(os.path.join(path_odfnew, ametal + '_*_' + str(G[ig]) + '.dat'))
                try:
                    max_temp = max(temps)
                except ValueError:
                    e = sys.exc_info()[0]
                    print 'Error is %s' % e, ametal, G[ig]
                    raise

                if (int(atemp1) > max_temp) or (int(atemp2) > max_temp) or (int(atemp3) > max_temp) or (int(atemp4) > max_temp):
                    indice = int(np.where(temps == max_temp)[0])
                    atemp1 = str(int(temps[indice - 3]))
                    atemp2 = str(int(temps[indice - 2]))
                    atemp3 = str(int(temps[indice - 1]))
                    atemp4 = str(int(temps[indice]))
                try:
                    tabA = read_file(os.path.join(path_odfnew,'%s_%s_%.1f.dat' % (ametal, atemp1, G[ig])), MT, N)
                    tabB = read_file(os.path.join(path_odfnew,'%s_%s_%.1f.dat' % (ametal, atemp2, G[ig])), MT, N)
                    tabC = read_file(os.path.join(path_odfnew,'%s_%s_%.1f.dat' % (ametal, atemp3, G[ig])), MT, N)
                    tabD = read_file(os.path.join(path_odfnew,'%s_%s_%.1f.dat' % (ametal, atemp4, G[ig])), MT, N)
                except:
                    nbreak += 1
                    break

            else:
                if len(glob.glob(path + '/' + ametal + '*.dat')) == 0:
                    files_m = glob.glob(path + '/a' + ametal + '*.dat')[0]
                    split_file(files_m, ametal, path)
                temps = glob.glob(os.path.join(path, ametal + '_*_' + str(G[ig]) + '.dat'))
                try:
                    max_temp = max(temps)
                except ValueError:
                    e = sys.exc_info()[0]
                    print 'Error is %s' % e, ametal, G[ig]
                    raise

                if (int(atemp1) > max_temp) or (int(atemp2) > max_temp) or (int(atemp3) > max_temp) or (int(atemp4) > max_temp):
                    indice = int(np.where(temps == max_temp)[0])
                    atemp1 = str(int(temps[indice - 3]))
                    atemp2 = str(int(temps[indice - 2]))
                    atemp3 = str(int(temps[indice - 1]))
                    atemp4 = str(int(temps[indice]))
                try:
                    tabA = read_file(os.path.join(path,'%s_%s_%.1f.dat' % (ametal, atemp1, G[ig])), MT, N)
                    tabB = read_file(os.path.join(path,'%s_%s_%.1f.dat' % (ametal, atemp2, G[ig])), MT, N)
                    tabC = read_file(os.path.join(path,'%s_%s_%.1f.dat' % (ametal, atemp3, G[ig])), MT, N)
                    tabD = read_file(os.path.join(path,'%s_%s_%.1f.dat' % (ametal, atemp4, G[ig])), MT, N)
                except:
                    nbreak += 1
                    break



            for i in range(MT):
                for j in range(N):
                    TabY[0] = tabA[i][j]
                    TabY[1] = tabB[i][j]
                    TabY[2] = tabC[i][j]
                    TabY[3] = tabD[i][j]
                    tck = UnivariateSpline(TabX, TabY, s = 1, k = 3)
                    Yfy = UnivariateSpline.__call__(tck, xteff)

                    if k == 0: tab1[i][j] = Yfy
                    if k == 1: tab2[i][j] = Yfy
                    if k == 2: tab3[i][j] = Yfy
                    if k == 3: tab4[i][j] = Yfy

                    del tck, Yfy


            del temps, atemp1, atemp2, atemp3, atemp4

        if nbreak == 0:

            output = open('./atm_models/' + starname + '.atm', 'w')

            output.writelines('KURUCZ\n')
            output.writelines('          Teff= %i          log g= %4.2f\n' % (xteff, xlogg))
            output.writelines('NTAU        72\n')

            for i in range(MT):
                for j in range(N - 2):
                    tab[i][j] = (tab1[i][j]**(1.-xfrac)*tab2[i][j]**xfrac)**(1.-yfrac)
                    tab[i][j] = tab[i][j]*(tab3[i][j]**(1.-xfrac)*tab4[i][j]**xfrac)**yfrac
                output.writelines(' %14.8E  %7.1f %9.3E %9.3E %9.3E %9.3E %9.3E\n' % \
                                 (tab[i][0], tab[i][1], tab[i][2], tab[i][3],\
                                 tab[i][4], tab[i][5], tab[i][6]))


            logefe = 7.50 + xmetal


            output.writelines('    %5.3fe+05\n' % micro)
            output.writelines('NATOMS     1  %5.2f\n' % xmetal)
            output.writelines('      26.0   %4.2f\n' % logefe)
            output.writelines('NMOL      19\n')
            output.writelines('      606.0    106.0    607.0    608.0    107.0    108.0    112.0    707.0\n')
            output.writelines('      708.0    808.0     12.1  60808.0  10108.0    101.0      6.1      7.1\n')
            output.writelines('        8.1    822.0     22.1')

            output.close()

        del TabX, nbreak

    del tabA, tabB, tabC, tabD, tab1, tab2, tab3, tab4, tab, TabY, T, G, XM
