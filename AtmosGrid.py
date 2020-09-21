#AtmGrid.py

import _pickle as pickle
import re, glob, os
import numpy as np

def split_file(file_m, ametal, path):
    i = 0
    with open(file_m) as f:
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

                new_file = open(os.path.join(path, '%s_%s_%s.dat' %
                                             (ametal, str(int(Teff)), str(logg))), 'w')
                i += 1
            try:
                new_file.writelines(line + '\n')
            except UnboundLocalError:
                pass

    try:
        new_file.close()
    except UnboundLocalError:
        pass


def read_file(f, MT, N):
    tabp = np.ones((MT, N))*np.nan
    try:
        mod = open(f)
    except FileNotFoundError:
        return tabp
    i = 0
    for line in mod:
        m = re.search('\A[0-9]\.[0-9]*', line)
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

def create_grid_file(file_grid):
    print('\t\tCreating the atmospheric model grid. This will be done just once.')
    T = np.array([3500., 3750., 4000., 4250., 4500., 4750., 5000., 5250., 5500.,\
                  5750., 6000., 6250., 6500., 6750., 7000., 7250., 7500., 7750.,\
                  8000., 8250., 8500., 8750., 9000., 9250., 9500., 9750.,\
                  10000., 10250., 10500., 10750., 11000., 11250., 11500., 11750.,\
                  12000., 12250., 12500., 12750., 13000., 14000., 15000.])
    G = np.linspace(0., 5., 11)
    XM = np.array([-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3,\
                   -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.5, 1.0])

    Tgg = []
    Ggg = []
    Mgg = []

    ii, jj, kk = len(T), len(G), len(XM)

    for i in range(ii):
        for j in range(jj):
            for k in range(kk):
                Tgg.append(T[i])
                Ggg.append(G[j])
                Mgg.append(XM[k])

    MT = 72
    N = 9

    col0 = []
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []

    for t,g,feh in zip(Tgg, Ggg, Mgg):
        signo='p'
        if feh < 0:
            signo = 'm'
        else:
            signp = 'p'
        
        avar = str(int(abs(feh*10)))
        if abs(feh*10) <= 5.:
            ametal = signo + '0' + avar
        else:
            ametal = signo + avar

        path = './atm_models/atlas9/grids/grid' + ametal
        path_odfnew = path + 'odfnew'
    
        if os.path.isdir(path_odfnew):
            path_file = path_odfnew
            if len(glob.glob(path_file + '/' + ametal + '*.dat')) == 0:
                files_m = glob.glob(path_file + '/a' + ametal + '*odfnew.dat')[0]
                split_file(files_m, ametal, path_file)

        else:
            path_file = path
            if len(glob.glob(path_file + '/' + ametal + '*.dat')) == 0:
                files_m = glob.glob(path_file + '/a' + ametal + '*.dat')[0]
                split_file(files_m, ametal, path_file)
    
        tab = read_file(os.path.join(path_file, '%s_%s_%.1f.dat' %
                                                (ametal, str(int(t)), g)), MT, N)
        col0.append(tab.T[0])
        col1.append(tab.T[1])
        col2.append(tab.T[2])
        col3.append(tab.T[3])
        col4.append(tab.T[4])
        col5.append(tab.T[5])
        col6.append(tab.T[6])
        del tab

    Tgg = np.array(Tgg)
    Mgg = np.array(Mgg)
    Ggg = np.array(Ggg)
    col0 = np.array(col0)
    col1 = np.array(col1)
    col2 = np.array(col2)
    col3 = np.array(col3)
    col4 = np.array(col4)
    col5 = np.array(col5)
    col6 = np.array(col6)

    dic = {'tgrid': Tgg, 'ggrid': Ggg, 'mgrid': Mgg,
           'col0': col0, 'col1': col1, 'col2': col2, 'col3': col3,\
           'col4': col4, 'col5': col5, 'col6': col6}

    with open(file_grid, 'wb') as grid:
        pickle.dump(dic, grid)

    del T, G, XM, Tgg, Mgg, Ggg, col0, col1, col2, col3, col4, col5, col6, dic


grid = None

# Check if grid file exists
file_grid = './atm_models/atlas9/ATLAS9_grid.pickle'
if not os.path.isfile(file_grid):
    create_grid_file(file_grid)

with open(file_grid, 'rb') as f:
    grid = pickle.load(f)

