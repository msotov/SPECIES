'''
Last modified in 17/05/2016
'''

import numpy as np
import sys

cmd = sys.argv
file_star = cmd[-1]

linelist = np.genfromtxt('../Spectra/linelist.dat', dtype = None, skip_header = 2,\
                         names = ('line', 'excit', 'loggf', 'num', 'ion'))
line = linelist['line']
excit = linelist['excit']
loggf = linelist['loggf']
num = linelist['num']
ion = linelist['ion']

file_ew = np.genfromtxt('../EW/' + file_star, dtype = None, usecols = (0, 4, 5),\
                        names = ('line', 'ew', 'ew_err'))
line_ew_b = file_ew['line']
ew_b = file_ew['ew']
ew_err_b = file_ew['ew_err']

#Take only the lines that 10. <= EW <=150 and the error in the EW is lower than the EW
ilines = np.where((ew_b >= 10.) & (ew_b <= 150.))[0]
line_ew_o = line_ew_b[ilines]
ew_o = ew_b[ilines]
ew_err_o = ew_err_b[ilines]

#line_ew_o = np.array([line_ew_b[i] for i in range(len(ew_b)) if (10. <= ew_b[i] <= 150.) == True])
#ew_o = np.array([ew_b[i] for i in range(len(ew_b)) if (10. <= ew_b[i] <= 150.) == True])
#ew_err_o = np.array([ew_err_b[i] for i in range(len(ew_b)) if (10. <= ew_b[i] <= 150.) == True])

div = ew_err_o/ew_o
indices = np.where(div <=1.)[0]
line_ew = line_ew_o[indices]
ew = ew_o[indices]
ew_err = ew_err_o[indices]

output = open('lines.' + file_star, 'w')
output.writelines(' %s\n' % file_star)

for i in range(len(line_ew)):
    index = np.where(line == line_ew[i])[0]
    if len(index) > 0:
        index = int(index)
        output.writelines('  %7.2f    %4.1f       %5.2f     %6.3f                       %7.4f\n' %\
                          (line[index], ion[index], excit[index], loggf[index], ew[i]))

output.close()

del linelist, line, excit, loggf, num, ion, file_ew, line_ew_b, ew_b, ew_err_b,\
    ilines, line_ew_o, ew_o, ew_err_o, indices, line_ew, ew, ew_err
