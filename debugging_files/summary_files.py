import glob, re
import numpy as np

files = glob.glob('*.log')
list_stars1 = np.genfromtxt('../list_cheps.dat', dtype = None, delimiter = ', ')
list_stars2 = np.genfromtxt('../list_s.dat', dtype = None, delimiter = ', ')

list_stars = np.append(list_stars1, list_stars2)

parent_p = '19730'

parent_process = []
star_name = []
num_star = []
process_id = []
finished = []

for s,f in enumerate(files):
    h = open(f)
    finished_file = False
    for i,linea in enumerate(h):
        if i == 0:
            cols = linea.split()
            star_name.append(cols[-1])
        elif i == 1:
            cols = linea.split()
            parent_process.append(cols[-1])
            del cols
        elif i == 2:
            cols = linea.split()
            process_id.append(cols[-1])
            del cols

        else:
            linea = linea.strip()
            m = re.search(r'.* run_iteration: Finished with star .*', linea)
            if m:
                finished_file = True

    h.close()

    finished.append(finished_file)

parent_process = np.array(parent_process)
star_name = np.array(star_name)
process_id = np.array(process_id)
finished = np.array(finished)

i_p = np.where(parent_process == parent_p)[0]

star_name = star_name[i_p]
process_id = process_id[i_p]
finished = finished[i_p]
#num_star = np.zeros(len(star_name))


#for i,s in enumerate(star_name):
#    i_s = int(np.where(list_stars == s)[0])
#    num_star[i] = i_s

p = np.unique(process_id)

#i_sort = np.argsort(num_star)
#star_name = star_name[i_sort]
#process_id = process_id[i_sort]
#finished = finished[i_sort]
#num_star = num_star[i_sort]

total_stars = 0
complete_processes = 0
running_processes = []

for process in p:
    print '\nProcess ' + process
    i_p = np.where(process_id == process)[0]
    print 'Number of stars in this process are ' + str(len(i_p))
    #print 'indices are:'
    #print num_star[i_p]
    finished_p = finished[i_p]
    non_finished_process = np.where(finished_p == False)[0]
    print 'Number of non-finished stars is ' + str(len(non_finished_process))
    print 'They are:'
    print star_name[i_p][non_finished_process]
    print '***************************************'

    total_stars += len(i_p)
    if len(non_finished_process) == 0:
        complete_processes += 1
    else:
        running_processes.append(process)


print '\nTotal number of stars completed or being processed is ' + str(total_stars)
print 'Completed processes: ' + str(complete_processes)
print 'Processes still running are '
print running_processes
print '***************************************'
