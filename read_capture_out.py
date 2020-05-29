import numpy as np
import os
nlay = 4
nrow = 342
ncol = 252

# NOte: In dataset, data are in layer, row ...
dt=np.loadtxt('capture_all.out')
print dt[4,10]

c = np.empty([nlay,nrow,ncol])

for k in range(dt.shape[0]):
    #print dt[k,1],dt[k,2],dt[k,3],dt[k,10]
    c[int(dt[k,1])-1,int(dt[k,2])-1,int(dt[k,3])-1]=dt[k,10]
#tmp=np.empty([nrow*ncol,1])
#tmp1=np.empty([nrow*ncol,nlay])

#for ilay in range(nlay):
    #print ilay
tmp1=np.reshape(c,(nrow*ncol,nlay))  # Reshape from 2D nrowxncol to 1D
tmp2=np.reshape(tmp1,nlay*nrow*ncol)
np.savetxt('dataset_part2.inp',tmp2,fmt='%8.3f')
os.system('copy dataset_part1.inp+dataset_part2.inp+dataset_part3.inp dataset1.dat')