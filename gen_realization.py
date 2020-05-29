import numpy as np

real_ID = np.loadtxt('real_ID.dat')

f=open('rnum.dat','r')
lines=f.readlines()
print lines[int(real_ID)]
f.close()	    

## Write location of the new pumping well to file
fid = open('para.dat', 'w')
fid.write('%s \n' %(line_out))  # Write another line
fid.close()
