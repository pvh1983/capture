# To generate a nwe MNW2 package
import os
import flopy
import numpy as np
#os.system('cls')

# Load river file (all river locations)
d1=np.loadtxt('river_cells.dat')
d2=np.loadtxt('all_active_cells.dat')


nrows, ncols = d1.shape
dtype={'names':['f{}'.format(i) for i in range(ncols)],
       'formats':ncols * [d1.dtype]}

C = np.intersect1d(d1.view(dtype), d2.view(dtype))
np.savetxt('new_river_cells',C,delimiter=' ',fmt='%4d')

exit()
fid = open('river_cells_new.dat', 'w')
for i in range(C[0,:]):
	fid.write('%4d %4d \n' %(C[i,0], C[i,1]))
fid.close()
exit()

