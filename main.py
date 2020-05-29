## NOTES:
## [1] Everytime when you switch feature_id, you have to run the code with 
## well_withdrawl=0 to get Qkt_before for all features. This is one time model run. 
## 

import os
import numpy as np
import numpy.ma as ma
import flopy.utils.binaryfile as bf





# input data
well_withdrawl=1  ## 0: before and 1: after
run_MODFLOW = 1  # 1: yes 0: No
inflow=1   # 1:yes, calculate flow to features (e.g rivers, lake); 0 or other numbers: flow from features
#feature_id=3   # 1: River, 2:GHB; 3: All river+GHB
print_notification_to_out_file=0  
osys = 2  ## 1: Windows; 2: Linux
Qwell=-100   # the rate of the new well withdrawal
total_well = 16  # that means only one new well is added
wberr = 1   # (L3/T) accepted error of (Total In - Total Out)
nRIV_cells = 3894
nGHB_cells = 507
if well_withdrawl==0:
	pumping_rate = 0  # Before well withdrawal
	print 'This run is for BEFORE the well withdrawal with pumping rate = 0.'
else:
	pumping_rate = Qwell  # Before well withdrawal
	#print 'This run is for AFTER well withdrawal with pumping rate = ', pumping_rate
print_river_flow_rates = 0  # 1: Print 0: No Print
print_Qkt_after = 1   # 1:yes, print 0: no
nlay = 4
nrow = 342
ncol = 252


## ================================================================
# BEFORE WELL WITHDRAWAL 
# Get these numbers by running the code with well_withdrawl=0 and
# check file Qkt_before.dat
## ================================================================
Q_entire_domain_before = 151021.0469
QMWN2_before = 22397.8691
Q_RVL_before = 100914.1016
Q_HDB_before = 19570.0654


### FUNCTION: Generate a new MWN2 file
def gen_mwn2(osys,total_well,pumping_rate,nw_lay,nw_row,nw_col):
    #print 'Pumping at row col top bot:',irow, icol, nw_top, nw_bot
	## Open file to write the first line
    fid = open('fw1.tmp', 'w')
    fid.write('%d %s \n' %(total_well, ' 740 0'))
    fid.close()

    ## Write information of the new well

    fid = open('fw3.tmp', 'w')
    fid.write('%s \n' %('\'NEW_WELL\' 1                                     # 2a.  WELLID,NNODES'))
    fid.write('%s \n' %('THIEM 0 0 0 0                                     # 2b.  LOSSTYPE,PUMPLOC,QLIMIT,PPFLAG,PUMPCAP'))
    fid.write('%s \n' %('1.0                                               # 2c.  Rw'))
    fid.write('%d %d %d %s \n' %(nw_lay,nw_row,nw_col,'                   # 2d2. LAY,ROW,COL'))
    fid.write('%d \n' %(total_well))  # Write the total number of wells
    fid.close()

    ## Write pumping rate of the new well
    fid = open('fw5.tmp', 'w')
    fid.write('%s %d \n' %('\'NEW_WELL\'',pumping_rate))  # Write another line
    fid.close()
    # Combine all files
    if osys == 1:  # for Windows
        os.system('copy fw1.tmp + fw2.inp + fw3.tmp + fw4.inp + fw5.tmp Future_SS.mnw2')
        os.system('del *.tmp')
    else:  # for Linux
        os.system('cat fw1.tmp fw2.inp fw3.tmp fw4.inp fw5.tmp > Future_SS.mnw2')
        os.system('rm -f *.tmp')

## FUNCTION: Read MODFLOW out file and extract information
def get_mfoutput():
    # Get TOTAL OUT from .out file
	with open("Future_SS.out") as f:
		c1 = 0   ## count variable
		c2 = 0
		c3 = 0
		aaa = np.empty([2,1])
		bbb = np.empty([2,1])
		ccc = np.empty([2,1])
		for line in f:
			if "           TOTAL OUT" in line:
				 line_with_TOTAL_OUT =  line.split()
#				 print 'line_with_TOTAL_OUT=',line_with_TOTAL_OUT
			elif "            IN - OUT =" in line:
				line_with_error =  line.split()
			elif "                MNW2 =" in line:
#				print 'c1=',c1			
				line_with_MWN2 =  line.split()
				aaa[c1]=float(line_with_MWN2[2])
#				print 'aaa=',aaa[c1]
				c1 = c1 + 1
			elif "     HEAD DEP BOUNDS =" in line:
#				print 'c2=',c2			
				line_with_HDB =  line.split()
				bbb[c2]=float(line_with_HDB[4])
#				print 'bbb=',bbb[c2]
				c2 = c2 + 1
			elif "       RIVER LEAKAGE =" in line:
#				print 'c3=',c3
				line_with_RVL =  line.split()
#				print 'line_with_RVL=',line_with_RVL
				ccc[c3]=float(line_with_RVL[3])
#				print 'ccc=',ccc[c3]
				c3 = c3 + 1

	QMWN2_=aaa[1]
	Q_entire_domain_ = float(line_with_TOTAL_OUT[3])
	err_ = float(line_with_error[4])
	Q_RVL_ = ccc[1]-ccc[0]
	Q_HDB_ = bbb[1]-bbb[0]
	#print 'err=',err	
	return Q_entire_domain_, err_, QMWN2_,Q_RVL_,Q_HDB_

	
	
### Load the new well information from para.dat
data = np.loadtxt('para.dat',delimiter=None)
# Notes: the numbers in para.dat are ilay,irow,icol,ib,nw_top,nw_bottom,flag
# print 'Pumping location is: ', data
ilay = data[1]
irow = data[2]
icol = data[3]
nw_top = data[5]
nw_bot = data[6]
### Generate a new well package
gen_mwn2(osys,total_well,pumping_rate,ilay,irow,icol)


## Run MODFLOW using the new well package
if run_MODFLOW==1:
	if osys == 1:  # for Windows
		os.system('MODFLOW-NWT_64.exe Future_SS.nam')
	else:  # for Linux
#	    print 'Running MODFLOW. Please wait ...'
#            os.system('ln -s /home/hpham/exe_files/mfnwt')
            os.system('./mfnwt Future_SS.nam > out.tmp')
	    Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB=get_mfoutput()
	    if abs(err) > wberr:
			print 'Warning: High residual between Total_In - Total_Out. Error (L3/T) = ', int(data[0]),int(data[1]),int(data[2]),int(data[3]), err
			os.system('rm -f Future_SS.nwt')
			os.system('rm -f Future_SS.out')
			os.system('cp -f /projects/DFN/sc/model_final/Future_SS_sml.nwt Future_SS.nwt')
			os.system('./mfnwt Future_SS.nam > out.tmp')
			Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB = get_mfoutput()	    	    
			print '                    Re-run MODFLOW and get (Total_In - Total_Out) = ', int(data[0]),int(data[1]),int(data[2]),int(data[3]), err
			#	    if print_notification_to_out_file==1:
#		print 'NEW MODFLOW run!'
else:
	print 'Warning: No NEW MODFLOW run! Please select run_MODFLOW=1 to run MODFLOW'

#print 'Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB = ', Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB


### Print information before the well withdrawal
if well_withdrawl==0:
	np.savetxt('Qkt_before.dat', [Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB], delimiter=',',fmt='%8.4f')
	print 'Done running the base scenario, no capture calculation. Existing ... '
	exit()

#if print_Qkt_after == 1:
#	np.savetxt('Qkt_after.dat', [Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB], delimiter=',',fmt='%8.4f')


Delta_Q_entire_domain = Q_entire_domain_before-Q_entire_domain
Delta_QMNW2=QMWN2_before-QMWN2
Delta_QRIV=-(Q_RVL_before-Q_RVL)
Delta_QGHB=-(Q_HDB_before-Q_HDB)
Delta_Qkt = Delta_QRIV+Delta_QGHB

## Write the final output to file
fo = np.empty([14])  # final outputs
fo[0:7]=data[0:7]  # Only 6 values, not 7
#print fo
#print data
fo[7]=Delta_Q_entire_domain
fo[8]=Delta_QMNW2      
fo[9]=Delta_QRIV       
fo[10]=Delta_QGHB      
fo[11]=Delta_Qkt       
if abs(Delta_QMNW2) > 1:
	fo[12]=Delta_Qkt/Delta_QMNW2   ## Capture from both rivers and lake
else:
	fo[12]=-999
	
fo[13]=abs(err)           #
#print 'CF_RIV | CF_GHD | CF_BOTH = ', Delta_QRIV/Delta_QMNW2, Delta_QGHB/Delta_QMNW2, Delta_Qkt/Delta_QMNW2
fid = open('capture.out', 'w')
fid.write('%6d %2d %3d %3d %2d %9.2f %9.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n' %(fo[0],fo[1],fo[2],fo[3],fo[4],fo[5],fo[6],fo[7],fo[8],fo[9],fo[10],fo[11],fo[12],fo[13]))  # Write another line
fid.close()
#np.savetxt('capture.out', final_ouput.transpose(), delimiter=',',fmt='%1.4f')

