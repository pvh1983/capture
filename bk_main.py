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
feature_id=3   # 1: River, 2:GHB; 3: All river+GHB
print_notification_to_out_file=0  
osys = 2  ## 1: Windows; 2: Linux
Qwell=-100   # the rate of the new well withdrawal
total_well = 16  # that means only one new well is added

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



## Function to generate a new MWN2 file
def gen_mwn2(osys,total_well,pumping_rate,nw_top,nw_bottom,nw_row,nw_col):
    print 'Pumping at row col top bot:',irow, icol, nw_top, nw_bot
	## Open file to write the first line
    fid = open('fw1.tmp', 'w')
    fid.write('%d %s \n' %(total_well, ' 740 0'))
    fid.close()

    ## Write information of the new well

    fid = open('fw3.tmp', 'w')
    fid.write('%s \n' %('\'NEW_WELL\' -1                                     # 2a.  WELLID,NNODES'))
    fid.write('%s \n' %('THIEM 0 0 0 0                                     # 2b.  LOSSTYPE,PUMPLOC,QLIMIT,PPFLAG,PUMPCAP'))
    fid.write('%s \n' %('1.0                                               # 2c.  Rw'))
    fid.write('%f %f %d %d %s \n' %(nw_top,nw_bottom,nw_row,nw_col,'                   # 2d2. Ztop,Zbotm,ROW,COL'))
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
		
## =====================================================
## Generate a new pumping file using data from para.dat
## =====================================================



# Load the new well information from para.dat
data = np.loadtxt('para.dat',delimiter=None)
# Notes: the numbers in para.dat are ilay,irow,icol,ib,nw_top,nw_bottom,flag
print 'Pumping location is: ', data
irow = data[2]
icol = data[3]
nw_top = data[5]
nw_bot = data[6]

gen_mwn2(osys,total_well,pumping_rate,nw_top,nw_bot,irow,icol)



## =====================================================
## Run MODFLOW using the new well package
## =====================================================
if run_MODFLOW==1:
	if osys == 1:  # for Windows
		os.system('MODFLOW-NWT_64.exe Future_SS.nam')
	else:  # for Linux
	    print 'Running MODFLOW. Please wait ...'
            os.system('./mfnwt Future_SS.nam > out.tmp')
	    os.system('rm -f fort.*')	    
	if print_notification_to_out_file==1:
		print 'NEW MODFLOW run!'
else:
	print 'Warning: No NEW MODFLOW run! Please select run_MODFLOW=1 to run MODFLOW'
 
## =====================================================
## Read MODFLOW ccf file using flopy
## =====================================================
	
	
	
	
## Construct a matrix of river conductance
dt2=np.loadtxt('river_info.dat')
cond = np.empty([nlay,nrow,ncol])
cond.fill(-999)
for k in range(dt2.shape[0]):
    cond[int(dt2[k,2])-1, int(dt2[k,0])-1, int(dt2[k,1])-1] = dt2[k,4]
	
	
	
# Load River Cells ID
if feature_id==1:
	rcid = np.loadtxt('river_cells.dat') # river cell id
elif feature_id==2:
	rcid = np.loadtxt('GHB_cells.dat') # river cell id
elif feature_id==3:
	rcid = np.loadtxt('river_n_GHB.dat') # river cell id
	
#print rcid.shape[0]
nFeatureCells=rcid.shape[0]

## Read MODFLOW ccf file using flopy
cbb = bf.CellBudgetFile('Future_SS.ccf')
#cbb.list_records()
#list_unique_records_ = cbb.list_unique_records()  # Print all RECORDS
#print 'nrecords=',cbb.get_nrecords()
#print '[.] Extract flow rate data from MODFLOW CCF file ====='
#CHD = cbb.get_data(totim=1, text='CONSTANT HEAD',  full3D=True)[0]
FRF = cbb.get_data(text='FLOW RIGHT FACE',full3D=True)[0]
FFF = cbb.get_data(text='FLOW FRONT FACE',full3D=True)[0]
FLF = cbb.get_data(text='FLOW LOWER FACE',full3D=True)[0]
#FBF = cbb.get_data(totim=1, text='FLOW BACK FACE',full3D=True)[0]
RLK = cbb.get_data(text='RIVER LEAKAGE',  full3D=True)[0]
HDB = cbb.get_data(text='HEAD DEP BOUNDS',full3D=True)[0]
RCH = cbb.get_data(text='RECHARGE',       full3D=True)[0]
MNW = cbb.get_data(text='MNW2',           full3D=True)[0]

# Convert masked element to zero
#FRF = np.where(FRF.mask, 0, FRF)
#FFF = np.where(FFF.mask, 0, FFF)
#FLF = np.where(FLF.mask, 0, FLF)
RLK = np.where(RLK.mask, 0, RLK)
HDB = np.where(HDB.mask, 0, HDB)
RCH = np.where(RCH.mask, 0, RCH)
MNW = np.where(MNW.mask, 0, MNW)

# ## Consider flow to feature only
# FRF[FRF>0]=0
# FFF[FFF>0]=0
# FLF[FLF>0]=0
# RLK[RLK>0]=0
# HDB[HDB>0]=0
# RCH[RCH>0]=0
# MNW[MNW>0]=0

#Qkt=FRF+FFF+FLF+RLK+HDB+RCH+MNW


## =====================================================
#  Start developing capture maps based on Leake et al.,[2010]
#  =====================================================

## [Step 1] Calculate Qkt and Q_entire_domain by running the model 
## without the added withdrawal

## [Step 2] Run the model with ADDED withdrawal

## [Step 3] Calculate Qkt, Q_entire_domain, Delta_Qkt and Delta_Q_entire_domain

# ## Select which features to simulate
if feature_id==1:
     Qkt = RLK   # Capture from rivers in the study polygon.
elif feature_id==2:
     Qkt = HDB
elif feature_id==3:
     Qkt = RLK+HDB



# Get TOTAL OUT from .out file
with open("Future_SS.out") as f:
    for line in f:
        if "TOTAL OUT" in line:
             line_with_TOTAL_OUT =  line.split()
	elif "PERCENT DISCREPANCY" in line:
	     line_with_error =  line.split()
	elif "                MNW2 =" in line:
	     line_with_MWN2 =  line.split()
	     aaa=float(line_with_MWN2[2])
QMWN2=aaa
Q_entire_domain = float(line_with_TOTAL_OUT[3])
err = float(line_with_error[3])


# xxx
QktEachCell = np.empty([nFeatureCells+1+1])  # Add two rows to write Q_entire_domain and QMWN2
QktEachCell[0] = Q_entire_domain
QktEachCell[1] = QMWN2
#if print_river_flow_rates == 1: 
if well_withdrawl==1:
	fid = open('budget_after.tmp', 'w')
else:
	fid = open('budget_before.tmp', 'w')
fid.write('%s \n' %('   i   ilay irow icol    FRFi      FFFi      FLFi       RLKi    HDBi    RCHi    MNWi    Qkti  Cond.'))
for i in range(nFeatureCells): # 
	#print(i,rcid[i-1,0]-1,rcid[i-1,1]-1, Qkt[irow,icol])
	irow = int(rcid[i,0]-1)  # Because of nparray's propteties
	icol = int(rcid[i,1]-1)
	ilay = int(rcid[i,2]-1)

	FRFi = FRF[ilay,irow,icol]
	FFFi = FFF[ilay,irow,icol]
	FLFi = FLF[ilay,irow,icol]
	RLKi = RLK[ilay,irow,icol]
	HDBi = HDB[ilay,irow,icol]
	RCHi = RCH[ilay,irow,icol]
	MNWi = MNW[ilay,irow,icol]	
	QktEachCell[i+1+1] = Qkt[ilay,irow,icol]	
	fid.write('%5d %4d %4d %4d %9.2f %9.2f %9.2f %9.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n' %(i,ilay+1, irow+1,icol+1,FRFi,FFFi,FLFi,RLKi,HDBi,RCHi,MNWi,QktEachCell[i+1+1],cond[ilay,irow,icol]))
#	QktEachCell[i+1] = Qkti # get Qkt given row and column

	#print Qkti
fid.close()
#print 'File budget.tmp is printed.'
#SumQktEachCell = QktEachCell.sum()
#print 'QktEachCell', QktEachCell
#print '[.] SumQktEachCell=',SumQktEachCell
#print '[.] Save QktEachCell in Qkt_before.dat'

if well_withdrawl==0:
	np.savetxt('Qkt_before.dat', QktEachCell, delimiter=',',fmt='%8.4f')
	print 'Done running the base scenario, no capture calculation. Existing ... '
	exit()
else:
	print 'One new well is added'

#if print_Qkt_after == 1:
#	np.savetxt('Qkt_after.dat', QktEachCell, delimiter=',',fmt='%8.4f')

Qafter = QktEachCell
#print 'Qafter=',Qafter

## Load Qkt_before.dat and calculate Delta_Q_entire_domain, Delta_Qkt
if well_withdrawl==1:
	Qbefore = np.loadtxt('../model_final/Qkt_before.dat')
#print 'Qbefore=',Qbefore[0:2]

Delta_Q_entire_domain = Qbefore[0]-Qafter[0]	
Delta_QMNW2 = Qbefore[1]-Qafter[1]	
print 'Delta_Q_entire_domain', Delta_Q_entire_domain

Delta_Qkt = Qbefore[2:]-Qafter[2:]		
#print 'Delta_Qkt=', Delta_Qkt.shape
#print 'Qbefore[1:]', Qbefore[1:]
#print 'Qafter=', Qafter[1:]

#Qwell_check = Delta_Q_entire_domain+Delta_Qkt.sum()  ## not CORRECT
#print 'Qwell_check=',Qwell_check

Capture = Delta_Qkt.sum()/Qwell
print 'Delta_Qkt.sum=', Delta_Qkt.sum()
print 'Capture=',Capture

#print 'Qwell=',Qwell

## Write the final output to file
fo = np.empty([12])  # final outputs
fo[0:7]=data[0:7]  # Only 6 values, not 7
#print fo
#print data
fo[7]=Delta_Q_entire_domain
fo[8]=Delta_Qkt.sum()
fo[9]=Delta_QMNW2
fo[10]=Capture
fo[11]=err
#print 'final_ouput=', fo
fid = open('capture.out', 'w')
fid.write('%6d %2d %3d %3d %2d %7.2f %7.2f %9.3f %9.3f %9.3f %9.3f %9.3f \n' %(fo[0],fo[1],fo[2],fo[3],fo[4],fo[5],fo[6],fo[7],fo[8],fo[9],fo[10],fo[11]))  # Write another line
fid.close()
#np.savetxt('capture.out', final_ouput.transpose(), delimiter=',',fmt='%1.4f')
