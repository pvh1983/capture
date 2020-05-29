import os
import numpy as np
import numpy.ma as ma
import flopy.utils.binaryfile as bf

## Function to generate a new MWN2 file
def gen_mwn2(osys,total_well,pumping_rate,nw_top,nw_bottom,nw_row,nw_col):
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

# input data
feature_id=1   # 1: River, 2:GHB; 3: All

osys = 2  ## 1: Windows; 2: Linux
run_MODFLOW = 1  # 1: yes 0: No
Qwell=-100   # the rate of the new well withdrawal
total_well = 16  # that means only one new well is added
well_withdrawl = 1  ## 0: before and 1: after
if well_withdrawl==0:
	pumping_rate = 0  # Before well withdrawal
	print 'This run is for BEFORE the well withdrawal with pumping rate = 0.'
else:
	pumping_rate = Qwell  # Before well withdrawal
	#print 'This run is for AFTER well withdrawal with pumping rate = ', pumping_rate
print_river_flow_rates = 0  # 1: Print 0: No Print
print_Qkt_after = 1   # 1:yes, print 0: no

# Load the new well information from para.dat
data = np.loadtxt('para.dat',delimiter=None)
# Notes: the numbers in para.dat are ilay,irow,icol,ib,nw_top,nw_bottom,flag

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
	    os.system('./mfnwt Future_SS.nam')
	    os.system('rm -f fort.*')
	    #print 'Done running MODFLOW.'
	print 'NEW MODFLOW run!'
else:
	print 'No NEW MODFLOW run!'
 
## =====================================================
## Read MODFLOW ccf file using flopy
## =====================================================
	
# Load River Cells ID
if feature_id==1:
	rcid = np.loadtxt('river_cells.dat') # river cell id
elseif feature_id==2:
	rcid = np.loadtxt('GHD_cells.dat') # river cell id
else:
	rcid = np.loadtxt('river_n_GHD.dat') # river cell id
	
#print rcid.shape[0]
nFeatureCells=rcid.shape[0]

## Read MODFLOW ccf file using flopy
cbb = bf.CellBudgetFile('Future_SS.ccf')
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

## =====================================================
#  Start developing capture maps based on Leake et al.,[2010]
#  =====================================================

## [Step 1] Calculate Qkt and QSt by running the model 
## without the added withdrawal

## [Step 2] Run the model with ADDED withdrawal

## [Step 3] Calculate Qkt, QSt, Delta_Qkt and Delta_QSt

#Qkt_tmp = FRF+FFF+FLF+RCH+RLK+MNW+HDB
Qkt = RLK[0, :, :]   # Capture from rivers in the study polygon.
#Qkt = np.flipud(Qkt_tmp[0, :, :])

# Get TOTAL OUT from .out file
with open("Future_SS.out") as f:
    for line in f:
        if "TOTAL OUT" in line:
             line_with_TOTAL_OUT =  line.split()
QSt = float(line_with_TOTAL_OUT[3])
#print 'QSt=', QSt


# xxx
QktEachCell = np.empty([nFeatureCells+1])  # Add one cell to write QSt
QktEachCell[0] = QSt
#if print_river_flow_rates == 1: 
fid = open('budget.tmp', 'w')
fid.write('%s \n' %(' i irow icol FRFi FFFi FLFi RLKi HDBi RCHi MNWi Qkti'))
for i in range(nFeatureCells): # xxx river cells
	#print(i,rcid[i-1,0]-1,rcid[i-1,1]-1, Qkt[irow,icol])
	irow = int(rcid[i,0]-1)  # Because of nparray's propteties
	icol = int(rcid[i,1]-1)
	FRFi = FRF[0,irow,icol]
	FFFi = FFF[0,irow,icol]
	FLFi = FLF[0,irow,icol]
	RLKi = RLK[0,irow,icol]
	HDBi = HDB[0,irow,icol]
	RCHi = RCH[0,irow,icol]
	MNWi = MNW[0,irow,icol]	
	Qkti = Qkt[irow,icol]
	fid.write('%5d %4d %4d %9.2f %9.2f %9.2f %9.2f %7.2f %7.2f %7.2f %7.2f\n' %(i,irow+1,icol+1,FRFi,FFFi,FLFi,RLKi,HDBi,RCHi,MNWi,Qkti))
	QktEachCell[i+1] = Qkti # get Qkt given row and column

	#print Qkti
fid.close()
#print 'File budget.tmp is printed.'
#SumQktEachCell = QktEachCell.sum()
#print 'QktEachCell', QktEachCell
#print '[.] SumQktEachCell=',SumQktEachCell
#print '[.] Save QktEachCell in Qkt_before.dat'

if well_withdrawl==0:
	np.savetxt('Qkt_before.dat', QktEachCell, delimiter=',',fmt='%8.4f')
	print 'Run the base scenario, no capture calculation. Existing ... '
	exit()
else:
	print 'One new well is added'

if print_Qkt_after == 1:
	np.savetxt('Qkt_after.dat', QktEachCell, delimiter=',',fmt='%8.4f')

Qafter = QktEachCell
#print 'Qafter=',Qafter

## Load Qkt_before.dat and calculate Delta_QSt, Delta_Qkt
if well_withdrawl==1:
	Qbefore = np.loadtxt('../model_before_withdrawl/Qkt_before.dat')
#print 'Qbefore=',Qbefore[0:2]

Delta_QSt = Qbefore[0]-Qafter[0]	
print 'Delta_QSt', Delta_QSt

Delta_Qkt = Qbefore[1:]-Qafter[1:]
#print 'Delta_Qkt=', Delta_Qkt.shape
#print 'Qbefore[1:]', Qbefore[1:]
#print 'Qafter=', Qafter[1:]

Qwell_check = Delta_QSt+Delta_Qkt.sum()
print 'Qwell_check=',Qwell_check

Capture = Delta_Qkt.sum()/Qwell
print 'Delta_Qkt.sum=', Delta_Qkt.sum()
print 'Capture=',Capture

#print 'Qwell=',Qwell

## Write the final output to file
fo = np.empty([11])  # final outputs
fo[0:7]=data[0:7]  # Only 6 values, not 7
#print fo
#print data
fo[7]=Delta_QSt
fo[8]=Delta_Qkt.sum()
fo[9]=Qwell_check
fo[10]=Capture
#print 'final_ouput=', fo
fid = open('capture.out', 'w')
fid.write('%6d %2d %3d %3d %2d %6.3f %6.3f %6.3f %6.3f %6.3f %9.5f \n' %(fo[0],fo[1],fo[2],fo[3],fo[4],fo[5],fo[6],fo[7],fo[8],fo[9],fo[10]))  # Write another line
fid.close()
#np.savetxt('capture.out', final_ouput.transpose(), delimiter=',',fmt='%1.4f')
