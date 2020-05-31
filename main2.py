import os
import numpy as np
import numpy.ma as ma
import flopy.utils.binaryfile as bf


'''
# Notes: 
- Use conda env b3new
- Check path to *.ccf *.out? [Check b4 submit]
- capture.out
- Run this from terminal: export source_path=$(pwd)

### Last updated:
12/20/2018
- After AGU, for the second round of UHR capture.
- This round, transient model with one ts of 2 years.
- 
05/28/2020
- Ro's model with 101 time steps
- 

### To run the codes 
- Prepare list of input files listed below
---- rnum_test.dat: List of all potential pmploc (submit.linux reads this file)
---- Modify fw2.inp and fw4.inp to generate a new .WEL file. 
     (These files are extracted from Capture_May2020.wel)
- Choose options to run
---- run the code once with well_withdrawl = 1 (by "python main_v31_trans.py")
---- submit the code (by "qsub submit.linux")



### If for a new project, DO:
- Change the name of model (search Capture_May2020 and replace it by YOUR_MODEL)
- 


# [] Input files 
- ID_Upper_Humboldt_River_Cells.dat to define the river location (cells).
- river_cells.dat to define the all river cells in the model domain.
- Change the name of Capture_May2020.wel

# Run sequence
- submit.linux
- main_v31_trans.py


'''

# Function to generate a new MWN2 file


def gen_mwn2(osys, total_well, pumping_rate, nw_lay, nw_row, nw_col):
    # Open file to write the first line
    fid = open('fw1.tmp', 'w')
    fid.write('%d  740 AUX IFACE AUX QFACT AUX CELLGRP \n' % (total_well))
    fid.write('%d           0         0 \n' % (total_well))
    fid.close()

    # Write information of the new well

    fid = open('fw3.tmp', 'w')
    fid.write('%d %4d %4d %4.2f %s \n' % (nw_lay, nw_row, nw_col,
                                          pumping_rate, ' 0 1.0 3745'))  # Write another line
    fid.close()

    # Combine all files
    if osys == 1:  # for Windows
        os.system('copy fw1.tmp + fw2.inp + fw3.tmp + fw4.inp TR_NWT_forward.wel')
        os.system('del *.tmp')
    else:  # for Linux
        os.system('cat fw1.tmp fw2.inp fw3.tmp fw4.inp > TR_NWT_forward.wel')
        os.system('rm -f *.tmp')

# FUNCTION: Read MODFLOW out file and extract information


def get_mfoutput(nts):
    # Get TOTAL OUT from .out file
    with open("TR_NWT_forward.out") as f:
        c0 = 0
        c1 = 0  # count variable
        c2 = 0
        c3 = 0
        aaa = np.empty([2*nts, 1])
        bbb = np.empty([2*nts, 1])
        ccc = np.empty([2*nts, 1])
        err_ = np.empty([2*nts, 1])
        for line in f:
            if "           TOTAL OUT" in line:
                line_with_TOTAL_OUT = line.split()
#				 print 'line_with_TOTAL_OUT=',line_with_TOTAL_OUT
            elif "            IN - OUT =" in line:
                line_with_error = line.split()
                # print 'line with error = ', line_with_error
                err_[c0] = float(line_with_error[9])
                # print 'line with error = ', err_[c0]
                c0 = c0 + 1
            elif "        WELLS =" in line:
                # print 'Lines with WEL = ',c1
                line_with_MWN2 = line.split()
                # print 'Lines with WEL = ', line_with_MWN2
                aaa[c1] = float(line_with_MWN2[5])
#				print 'WEL in/out =',aaa[c1]
                c1 = c1 + 1
#			elif "     HEAD DEP BOUNDS =" in line:
#				print 'c2=',c2
#				line_with_HDB =  line.split()
#				bbb[c2]=float(line_with_HDB[4])
#				print 'bbb=',bbb[c2]
#				c2 = c2 + 1
#			elif "       RIVER LEAKAGE =" in line:
                # print 'c3=',c3
#				line_with_RVL =  line.split()
#				print 'line_with_RVL=',line_with_RVL
#				ccc[c3]=float(line_with_RVL[3])
                # print 'ccc=',ccc[c3]
#				c3 = c3 + 1

#	QMWN2_=aaa[1]
    Q_entire_domain_ = float(line_with_TOTAL_OUT[3])

#	Q_RVL_ = ccc[1]-ccc[0]
#	Q_HDB_ = bbb[1]-bbb[0]
    # print 'err=',err_
    return Q_entire_domain_, err_


# =====================================================
# Generate a new pumping file using data from pmploc.dat
# =====================================================

print '   '
print '   '

# input data
feature_id = 1   # 1: River, 2:GHB; 3: All

run_MODFLOW = 1  # 1: yes 0: No
well_withdrawl = 1  # 0: before and 1: after

Qwell = -20000   # the rate of the new well withdrawal
total_well = 3745  # that means only one new well is added
osys = 2  # 1: Windows; 2: Linux

if well_withdrawl == 0:
    pumping_rate = 0  # Before well withdrawal
#	print 'This run is for BEFORE the well withdrawal with pumping rate = 0.'
else:
    pumping_rate = Qwell  # Before well withdrawal
    # print 'This run is for AFTER well withdrawal with pumping rate = ', pumping_rate
print_river_flow_rates = 0  # 1: Print 0: No Print
print_Qkt_after = 0   # 1:yes, print 0: no

# Load the new well information from pmploc.dat
data = np.loadtxt('pmploc.dat', delimiter=None)
# Notes: the numbers in pmploc.dat are ID,ilay,irow,icol for WEL package

# Load cell IDs of the Upper Humboldt River
id_tmp = np.loadtxt('ID_Upper_Humboldt_River_Cells.dat', delimiter=None)
# print ID_HRiv_cells
# Convert to int. Minus 1 b/c python array 0
ID_HRiv_cells = [int(y-1) for y in id_tmp]

rid = data[0]
ilay0 = data[1]
irow0 = data[2]
icol0 = data[3]

gen_mwn2(osys, total_well, pumping_rate, ilay0, irow0, icol0)

# =====================================================
# Run MODFLOW using the new well package
# =====================================================
if run_MODFLOW == 1:
    if osys == 1:  # for Windows
        os.system('MODFLOW-NWT_64.exe TR_NWT_forward.nam')
    else:  # for Linux
        os.system('./mfnwt TR_NWT_forward.nam')
        os.system('rm -f fort.*')
        # print 'Done running MODFLOW.'
#	print 'NEW MODFLOW run!'
else:
    print 'No NEW MODFLOW run!'


# =====================================================
# Read MODFLOW ccf file using flopy
# =====================================================

# Load River Cells ID
if feature_id == 1:
    rcid = np.loadtxt('river_cells.dat')  # river cell id
elif feature_id == 2:
    rcid = np.loadtxt('GHD_cells.dat')  # river cell id
else:
    rcid = np.loadtxt('river_n_GHB.dat')  # river cell id

# print rcid.shape[0]
nFeatureCells = rcid.shape[0]

# Read MODFLOW ccf file using flopy
cbb = bf.CellBudgetFile('TR_NWT_forward.ccf')

# list_unique_records = cbb.list_unique_records()  # Print all RECORDS
# print 'nrecords=',cbb.get_nrecords()

# Get a list of unique stress periods and time steps in the file
list_stress_periods = cbb.get_kstpkper()
# print 'list_stress_periods = ', cbb.get_kstpkper()

nts = len(list_stress_periods)
# print 'Number of stress periods = ', nts


Q_entire_domain, err = get_mfoutput(nts)
# print 'Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB = ', Q_entire_domain, err, QMWN2, Q_RVL, Q_HDB
#err_max = err
# print 'Error of water budget for all time steps (L3/T) = ', err[0:nts]


# list_unique_times = cbb.get_times()  # Get a list of unique times in the file
# print 'list_unique_times = ', cbb.get_times()

print ' '
print 'Extracting  flow rate data from MODFLOW CCF file ... '
n_out_rows = 9  # additional rows to print STO, CHD, DRN, RLV_all, EVT, HDB, RCH,HRiv, and ORiv
Qafter = np.empty([n_out_rows])  #
Qafter_all_ts = np.empty([n_out_rows, nts])  #
QktEachFeature_all_ts = np.empty([n_out_rows, nts])  #
STO_all_ts = np.empty([nts])
WEL_all_ts = np.empty([nts])
QktEachFeature = np.empty([n_out_rows])  #
QktEachFeature_tmp = np.empty([nFeatureCells])  #


if well_withdrawl == 1:
    Qbefore = np.loadtxt('../model_final/Qkt_before.dat', delimiter=',')
    #QStbefore = np.loadtxt('../model_final/QSt_before.dat',delimiter=',')
# print 'Qbefore = ', Qbefore.shape

print ' '
# print 'Testing ...'

#fid = open('budget.tmp', 'w')
fid1 = open('capture.out', 'w')
for ts in range(nts):  # number of stress periods
    #	print 'ts = ',ts
    if ts < nts-1:  # Because the last time step is set to steady state
        STO = cbb.get_data(kstpkper=(0, ts), text='STORAGE',  full3D=True)[0]
        STO_all_ts[ts] = STO.sum()
        QktEachFeature[0] = STO.sum()
    else:
        STO = 0
        STO_all_ts[ts] = 0
        QktEachFeature[0] = 0

    CHD = cbb.get_data(kstpkper=(0, ts), text='CONSTANT HEAD',  full3D=True)[0]
    FRF = cbb.get_data(
        kstpkper=(0, ts), text='FLOW RIGHT FACE', full3D=True)[0]
    FFF = cbb.get_data(
        kstpkper=(0, ts), text='FLOW FRONT FACE', full3D=True)[0]
    FLF = cbb.get_data(
        kstpkper=(0, ts), text='FLOW LOWER FACE', full3D=True)[0]
    WEL = cbb.get_data(kstpkper=(0, ts), text='WELLS',          full3D=True)[0]
    DRN = cbb.get_data(kstpkper=(0, ts), text='DRAINS',         full3D=True)[0]
    RLK = cbb.get_data(kstpkper=(0, ts), text='RIVER LEAKAGE',  full3D=True)[0]
    EVT = cbb.get_data(kstpkper=(0, ts), text='ET',             full3D=True)[0]
    HDB = cbb.get_data(
        kstpkper=(0, ts), text='HEAD DEP BOUNDS', full3D=True)[0]
    RCH = cbb.get_data(kstpkper=(0, ts), text='RECHARGE',       full3D=True)[0]

    # Convert masked element to zero
    #FRF = np.where(FRF.mask, 0, FRF)
    #FFF = np.where(FFF.mask, 0, FFF)
    #FLF = np.where(FLF.mask, 0, FLF)
    #STO = np.where(STO.mask, 0, STO)
    RLK = np.where(RLK.mask, 0, RLK)
    HDB = np.where(HDB.mask, 0, HDB)
    CHD = np.where(CHD.mask, 0, CHD)
    RCH = np.where(RCH.mask, 0, RCH)
    WEL = np.where(WEL.mask, 0, WEL)
    DRN = np.where(DRN.mask, 0, DRN)


#	print('ts = %3d, STORAGE CHANGE = %9.1f' %(ts+1, STO.sum()))

    WEL_all_ts[ts] = WEL.sum()
#	print('ts = %3d, WEL.sum = %9.1f' %(ts+1, WEL.sum()))
    # exit

    # =====================================================
    #  Start developing capture maps based on Leake et al.,[2010]
    #  =====================================================

    # [Step 1] Calculate Qkt and QSt by running the model
    # without the added withdrawal

    # [Step 2] Run the model with ADDED withdrawal

    # [Step 3] Calculate Qkt, QSt, Delta_Qkt and Delta_STO

    #Qkt_tmp = FRF+FFF+FLF+RCH+RLK+WEL+HDB
    # Qkt = RLK[0, :, :]   # Capture from rivers in the study polygon.
    # Qkt = RLK[0, :, :] + HDB[0, :, :]   # Capture from rivers + Lake layer 1 only
#	Qkt = RLK   # Capture from rivers + Lake layer 1 only
    #Qkt = np.flipud(Qkt_tmp[0, :, :])

    # Get TOTAL OUT from .out file
#	with open("TR_NWT_forward.out") as f:
#	    for line in f:
#	        if "TOTAL OUT" in line:
#	             line_with_TOTAL_OUT =  line.split()
#	QSt = float(line_with_TOTAL_OUT[3])
#	#print 'QSt=', QSt
    QktEachFeature[1] = CHD.sum()
    QktEachFeature[2] = WEL.sum()
    QktEachFeature[3] = DRN.sum()
    QktEachFeature[4] = RLK.sum()
    QktEachFeature[5] = EVT.sum()
    QktEachFeature[6] = RCH.sum()

#	print 'QktEachFeature[0:6] = ', QktEachFeature[0:6]
    # if print_river_flow_rates == 1:
#	fid.write('Time step = %s \n' %(ts+1))
#	fid.write('%s \n' %('i  irow icol ilay Right     Left      Front     Back     Lower      RLKi        HDBi    RCHi        WELi     Qkti '))
    for i in range(nFeatureCells):  # xxx river cells
        #print(i,rcid[i-1,0]-1,rcid[i-1,1]-1, Qkt[irow,icol])
        irow = int(rcid[i, 0]-1)  # Because of nparray's propteties
        icol = int(rcid[i, 1]-1)
        ilay = int(rcid[i, 2]-1)
#		Right = -FRF[ilay,irow,icol]  # RIGHT FACE
#		Left = FRF[ilay,irow,icol-1]  # LEFT FACE
#		Front = -FFF[ilay,irow,icol]  # FRONT FACE
#		Back = FFF[ilay,irow-1,icol]  # BACK FACE
#		Lower = -FLF[ilay,irow,icol]  # LOWER FACE
#		RLKi = RLK[ilay,irow,icol]    #
#		HDBi = HDB[ilay,irow,icol]
#		RCHi = RCH[ilay,irow,icol]
#		WELi = WEL[ilay,irow,icol]
#		Qkti = RLK[ilay,irow,icol]
#		fid.write('%5d %4d %4d %4d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n' %(i,irow+1,icol+1,ilay+1,Right,Left,Front,Back,Lower,RLKi,HDBi,RCHi,WELi,Qkti))
        QktEachFeature_tmp[i] = RLK[ilay, irow, icol]
        # print i, Qkti, i+n_out_rows

        # print Qkti

    Q_HRiv = QktEachFeature_tmp[ID_HRiv_cells].sum()
    # Sum of RLK from all other cells (not belong to the Humboldt river)
    Q_ORiv = QktEachFeature_tmp.sum() - Q_HRiv

    QktEachFeature[7] = Q_HRiv  # Humboldt river
    QktEachFeature[8] = Q_ORiv  # Other river
    QktEachFeature_all_ts[:, ts] = QktEachFeature

    # print 'File budget.tmp is printed.'
    #SumQktEachFeature = QktEachFeature.sum()
#	print 'QktEachFeature', QktEachFeature.shape
    # print '[.] SumQktEachFeature=',SumQktEachFeature
    # print '[.] Save QktEachFeature in Qkt_before.dat'


#	if print_Qkt_after == 1:
#		np.savetxt('Qkt_after.dat', QktEachFeature, delimiter=',',fmt='%9.5f')

    Qafter = QktEachFeature
    Qafter_all_ts[:, ts] = QktEachFeature
    # print 'Qafter=',Qafter

    # Load Qkt_before.dat and calculate Delta_STO, Delta_Qkt
    if well_withdrawl == 1:
        #Qbefore = np.loadtxt('../model_final/Qkt_before.dat',delimiter=',')
        #QStbefore = np.loadtxt('../model_final/QSt_before.dat',delimiter=',')
        # print 'Qbefore=',Qbefore[0:2]
        Delta_STO = Qbefore[0, ts]-Qafter[0]
        # print 'Change in storage Delta_STO = ', Delta_STO

        Delta_CHD = Qafter[1] - Qbefore[1, ts]
        Delta_WEL = Qafter[2] - Qbefore[2, ts]
        Delta_DRN = Qafter[3] - Qbefore[3, ts]
        Delta_RLK = Qafter[4] - Qbefore[4, ts]
        Delta_EVT = Qafter[5] - Qbefore[5, ts]
        Delta_RCH = Qafter[6] - Qbefore[6, ts]
        Delta_HRi = Qafter[7] - Qbefore[7, ts]  # Humboldt river
        Delta_ORi = Qafter[8] - Qbefore[8, ts]  # Other rivers
#		print 'Change in pumping rate Delta_WEL=', Delta_WEL
#		print 'Qafter = ', Qafter.shape
#		Delta_Qkt = Qbefore[n_out_rows:,ts]-Qafter[n_out_rows:]
        # print 'Sum of changes in all other features Delta_Qkt=', Delta_Qkt.sum()

        # print 'Qbefore[1:]', Qbefore[1:]
        # print 'Qafter=', Qafter[1:]

        #Qwell_check = Delta_STO+Delta_Qkt.sum()
#		print('Qwell_check=%5.2f. (Must equal %5.2f)'%(Qwell_check,pumping_rate))

#		Capture = Delta_STO/Qwell + Delta_Qkt.sum()/Qwell
        # print 'Delta_Qkt.sum=', Delta_Qkt.sum()
        # print 'Capture=',Capture
#		if ts==0:
#			print ' '
#			print 'Estimated capture fractions:'
#			print '  TS   Qwell_ask QWel_act STO_befor STO_after   D_QSt      D_Qkt    Capture   error  '
#		if ts < nts-1: # Because the last time step is set to steady state
#			print('%4d %9.1f %9.1f %9.1f %9.1f %9.2f %9.2f %9.2f %7.2f ' %(ts+1, Qwell, Delta_WEL, Qbefore[0,ts], STO.sum(), Delta_STO, Delta_Qkt.sum(), Capture, err[ts]))
#		else:
#			print('%4d %9.1f %9.1f %9.1f %9.1f %9.2f %9.2f %9.2f %7.2f ' %(ts+1, Qwell, Delta_WEL, Qbefore[0,ts], 0, Delta_STO, Delta_Qkt.sum(), Capture, err[ts]))
        # print 'Qwell=',Qwell

        # Write the final output to file
#		fo = np.empty([8])  # final outputs
#		fo[0:4]=data[0:4]  # 0:4 to get 4 columns
        # print fo[0:4]
        # print data
#		fo[4]=Delta_STO
#		fo[5]=Delta_Qkt.sum()
#		fo[6]=Capture
#		fo[7]=err[ts]
        # print 'final_ouput=', fo
#		if ts==0:
#			fid1.write('Qwell = %9.1f (L3/T)\n' %( Qwell ))
#			fid1.write('    ID   TS  Lay  Row  Col Capture err[ts] D_STO D_all_features D_WEL D_DRN D_RLK D_EVT D_RCH \n')
        # fid1.write('%6d %4d %4d %4d %4d %9.1f %9.1f %12.3f %12.3f %9.1f %9.3f %9.2f \n' %(fo[0],ts+1,fo[1],fo[2],fo[3],Qbefore[0,ts], STO.sum(), fo[4],fo[5], Delta_WEL, fo[6],fo[7]))  # Write another line
        fid1.write('%6d, %4d, %4d, %4d, %4d, %9.3f, %9.3f, %9.3f, %9.3f, %9.3f, %9.3f, %9.3f, %9.3f, %9.3f \n'
                   % (rid, ts+1, ilay0, irow0, icol0, err[ts], Delta_STO, Delta_WEL, Delta_DRN, Delta_RLK, Delta_HRi, Delta_ORi, Delta_EVT, Delta_RCH))  # Write another line
# fid.close()
fid1.close()
#np.savetxt('capture.out', final_ouput.transpose(), delimiter=',',fmt='%1.4f')


if well_withdrawl == 0:
    np.savetxt('Qkt_before.dat', QktEachFeature_all_ts,
               delimiter=',', fmt='%9.5f')
#	np.savetxt('QSt_before.dat', STO_all_ts, delimiter=',',fmt='%8.1f')
    print 'Run the base scenario, no capture calculation. Existing ... '
    exit()
else:
    print 'One new well is added'

#np.savetxt('Qkt_after_all.dat', Qafter_all_ts, delimiter=',',fmt='%8.4f')
print 'Check file capture.out for the outputs.'
