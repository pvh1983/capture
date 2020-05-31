# CHECK if this is an older verion.

from __future__ import division
import flopy
import numpy as np
import os
import matplotlib.pyplot as plt

from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
from shutil import copyfile


os.system('set GDAL_DATA=C:\Miniconda2\pkgs\libgdal-2.1.0-vc9_0\Library\share\gdal')

ifile = 'capture_all.out'   # the output file after runing simulations on Oasis

nlay = 3  # must equal number of MODFLOW layer for dataset
nrow = 770
ncol = 405
nts = 7

xorg = 1902510.0
yorg = 14487900.0
cell_size = 900   # feet


thre_Qact = 19775  # Just consider locations where can pump more than this theshold
nts_to_run = 7  # Time step to get outputs 1-7


# Load cell IDs of the Upper Humboldt River
ID_HRiv_cells = np.loadtxt('ID_Upper_Humboldt_River_Cells.dat', delimiter=None)
# print ID_HRiv_cells

# Convert UTM to decimal degree

# def transform_utm_to_wgs84(easting, northing, zone):
#     utm_coordinate_system = osr.SpatialReference()
#     utm_coordinate_system.SetWellKnownGeogCS("WGS84")  # Set geographic coordinate system to handle lat/lon
#     is_northern = northing > 0
#     utm_coordinate_system.SetUTM(zone, is_northern)
#
#     wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS()  # Clone ONLY the geographic coordinate system
#
#     # create transform component
#     utm_to_wgs84_geo_transform = osr.CoordinateTransformation(utm_coordinate_system, wgs84_coordinate_system)  # (, )
#     return utm_to_wgs84_geo_transform.TransformPoint(easting, northing, 0)  # returns lon, lat, altitude
#

# def transform_wgs84_to_utm(lon, lat):
#     def get_utm_zone(longitude):
#         return (int(1 + (longitude + 180.0) / 6.0))
#
#     def is_northern(latitude):
#         """
#         Determines if given latitude is a northern for UTM
#         """
#         if (latitude < 0.0):
#             return 0
#         else:
#             return 1

# utm_coordinate_system = osr.SpatialReference()
# utm_coordinate_system.SetWellKnownGeogCS("WGS84")  # Set geographic coordinate system to handle lat/lon
# utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))
#
# wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS()  # Clone ONLY the geographic coordinate system
#
# # create transform component
# wgs84_to_utm_geo_transform = osr.CoordinateTransformation(wgs84_coordinate_system, utm_coordinate_system)  # (, )
# return wgs84_to_utm_geo_transform.TransformPoint(lon, lat, 0)  # returns easting, northing, altitude


def print_2D_data():
    # Print x y and CF to a file
    for k2 in range(nlay):
        fid3 = open('CF.tmp', 'w')
        fid3.write('X, Y, Z \n')
        ofile_cf = 'CF2D_fID_' + \
            str(feature_id) + '_ts_' + str(ts + 1) + \
            '_lay_' + str(k2 + 1) + '.txt'
        for i2 in range(nrow):
            for j2 in range(ncol):
                # print k,i,j
                fid3.write('%9.1f, %9.1f, %4.3f \n' %
                           (x_coor[j2], y_coor[i2], c[k2, i2, j2]))
        fid3.close()
        copyfile('CF.tmp', ofile_cf)

    # Interpolation


def interpolation_cell(c, nlay, nrow, ncol):
    c2 = c
    d9 = np.empty([9, 1])
    for k2 in range(nlay):
        # print 'k2 = ', k2
        for i2 in range(1, nrow-1, 1):
            for j2 in range(1, ncol-1, 1):
                d9 = [c[k2, i2-1, j2-1], c[k2, i2-1, j2], c[k2, i2-1, j2+1],
                      c[k2, i2, j2-1], c[k2, i2, j2], c[k2, i2, j2+1],
                      c[k2, i2+1, j2-1], c[k2, i2+1, j2], c[k2, i2+1, j2+1]]
                # print 'd9 = ',d9
                d9n = [i for i in d9 if i >= 0.0]  # Find value -999
                if d9n:  # If d9 is NOT empty
                    # print feature_id,ts+1,k2,i2,j2,len(d9n),sum(d9n)/len(d9n)
                    c2[k2, i2, j2] = sum(d9n)/len(d9n)


# Generate lat long
x_coor = np.empty([ncol])
y_coor = np.empty([nrow])
# array=c[0,:,:]
# array[array==-999]=-1

# print 'MIN_array:', array.min()
# print 'MAX_array:', array.max()

#lat  = np.empty([nrow,ncol])
#long = np.empty([nrow,ncol])

for i in range(nrow):
    y_coor[i] = (yorg+nrow*cell_size)-0.5*cell_size-(i-1+1)*cell_size

for j in range(ncol):
    x_coor[j] = (j-1+1)*cell_size+0.5*cell_size+xorg

#    print aaa
    #aaa=(j - 1 + 1) * cell_size + 0.5 * cell_size + xorg
    #x_coor[0,j],tmp_x,z_tmpx = transform_utm_to_wgs84(aaa,4322050,11)
# print x_coor[0,0]


#    print  bbb
    #bbb = (yorg + nrow * cell_size) - 0.5 * cell_size - (i - 1 + 1) * cell_size
    #tmp_x, y_coor[0,i], z_tmpy = transform_utm_to_wgs84(745850,bbb,11)
# print x_coor[0,0], y_coor[0,0]
# print x_coor, y_coor

lat = np.repeat(x_coor, nrow, axis=0)
long_ = (np.repeat(y_coor, ncol, axis=0))
lon = long_.transpose()
# print lat[0,0]
# print lon[0,0]



# Calculate model run times
all_pmploc = 24532*2
time_ = np.empty([all_pmploc, 1])
with open("hpham.out") as f:
    ct = 0
    for line in f:
        if "Elapsed run time:" in line:
            line_with_time = line.split()
#            print 'line_with_time=',line_with_time
            time_min = float(line_with_time[3])
            time_sec = float(line_with_time[5])
            total_time = time_min + time_sec/60
            # print total_time
#            if total_time >= 3:
            time_[ct] = total_time
            ct = ct + 1
# time_ = time_[~np.all(time_ == 0, axis=0)] # Remove all zero rows
time_ = np.delete(time_, np.where(~time_.any(axis=1))[0], axis=0)
np.savetxt('time.dat', time_, delimiter='None', fmt='%3.1f')
print('n model runs, mean and 1-std model run time = ',
      ct-1, time_.mean(), time_.std())

# the histogram of the data
plt.hist(time_, 50, density=True, facecolor='g', alpha=0.75)
plt.xlabel('MODFLOW model run time (mins)')
plt.ylabel('Probability')
plt.title('Histogram of model run times. Nsamples = %d' % (time_.shape[0]))
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([8, 10, 0, 2])
plt.grid(True)
# plt.show()
plt.savefig('mruntime.png')

# the histogram of the data


def plot_hist(data1, time_step):
    # print abs(data)
    plt.hist(abs(data), 10, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('Withdrawal Rate (ft^3/day)')
    plt.ylabel('Probability')
    plt.title('Time step = %2d. Nsamples = %d \n.' %
              (time_step+1, data1.shape[0]))
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([8, 10, 0, 2])
    plt.grid(True)
    ofile2 = 'Qact_TS_' + str(time_step) + '_.png'
    plt.savefig(ofile2)
    # plt.show()


def plot_hist_err(data2):
    # print abs(data)
    plt.hist(abs(data2), 10, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('Water Balance Errors (ft^3/day)')
    plt.ylabel('Probability')
    plt.title('Nsamples = %d \n.' % (data2.shape[0]))
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([8, 10, 0, 2])
    # plt.grid(True)
    plt.savefig('Err.png')
    # plt.show()



fid1 = open('outputs.txt', 'w')
fid1.write('Capture Fraction Calculation \n')
fid1.write('Input file = %s \n' % (ifile))
fid1.write('Q_threshold =  %d (ft3/day)\n' % (thre_Qact))
fid1.write('Mean of model run times =  %3.1f (mins)\n' % (time_.mean()))
fid1.write('1std of model run times =  %3.1f (mins)\n' % (time_.std()))



ofile3 = 'Qthre_' + str(thre_Qact) + '.txt'
fid3 = open(ofile3, 'w')

# for feature_id in range(1,5,1):   ## 4: error
for feature_id in range(1, 7, 1):
    print('feature_id = {feature_id} \n')
    # 1: Humboldt, 2: Other rivers, 3: all rivers, 4: all features
    # 5: Water budget error; 6: Actual pumping rate

    if feature_id == 1:  # ONLY the Upper Humboldt River
        feature_name = 'CF Humboldt'
    elif feature_id == 2:  # From other rivers (not the Upper Humboldt River)
        feature_name = 'CF Tributaries'
    elif feature_id == 3:
        feature_name = 'CF Main Stem and Tributaries'
    elif feature_id == 4:
        feature_name = 'CF_AFea'
    if feature_id == 5:
        feature_name = 'WB_err'
    elif feature_id == 6:
        feature_name = 'Qwell_act'

    fid1.write('Feature ID = %2d \n' % (feature_id))
    # Construct a matrix of river conductance
    dt2 = np.loadtxt('river_info.dat')
    cond = np.empty([nlay, nrow, ncol])
    cond.fill(-999)
    for k in range(dt2.shape[0]):
        cond[int(dt2[k, 2])-1, int(dt2[k, 0])-1, int(dt2[k, 1])-1] = dt2[k, 4]

    # NOte: In dataset, data are in layer, row ...
    # dt=np.loadtxt('capture_river.out')

    # dt=np.loadtxt('capture_GHB.out')
    # dt=np.loadtxt('capture_river_NEW.out')
    #Qwell = dt[:,9]
    # Qwell[Qwell==0]=-1e12
    # print Qwell
#    dt=np.loadtxt('capture_river_n_GHB.out') # Tahoe
    data = np.loadtxt(ifile, delimiter=',')  # Tahoe
    [dt_row, dt_col] = data.shape
    print('Number of rows in the ifile = ', dt_row)
    list_of_ts = data[:, 1]
    if feature_id == 1:
        plot_hist_err(data[:, 5])
    # Extract data for each stress period
    for ts in range(nts_to_run):
        # for ts in range(6, 7, 1): # max is 7 but 0-8
        #        print 'fID/TS = ',feature_id,ts+1
        fid1.write('Time step = %2d \n' % (ts+1))
        id_tmp = np.where(list_of_ts == ts+1)
        # print id_tmp
        dt = data[id_tmp]
#        if feature_id == 1:
#            plot_hist(dt[:,7], ts)
        [npmploc, tmp] = dt.shape
#        print 'Number of pumping locations = ', npmploc
        print('fID/TS/Qthred/nloc = ', feature_id, ts + 1, thre_Qact, npmploc)
        c = np.empty([nlay, nrow, ncol])  # Capture Fraction
        c[:, :, :] = -999
        c.fill(-999)
        ctQ = 1
        ctQ2 = 1
        ctQ3 = 1
        for k in range(npmploc):
            # print 'k = ', k+1
            if k == 0:  # Print title row
                fid1.write(
                    'ts    id      D_STOR        D_WEL      D_DRN   D_RLK     D_EVT     D_CHD     CF    Err(%) \n')

            if abs(dt[k, 7]) >= thre_Qact:   # Just consider cells that Q_wellwithdrawl > thre_Qact
                if feature_id == 1:   # ONLY the Upper Humboldt River
                    CF = -(dt[k, 10]) / dt[k, 7]
                # From other rivers (not the Upper Humboldt River)
                elif feature_id == 2:
                    CF = -(dt[k, 11]) / dt[k, 7]
                elif feature_id == 3:
                    CF = -(dt[k, 9]) / dt[k, 7]
                elif feature_id == 4:
                    CF = -(dt[k, 8]+dt[k, 9]+dt[k, 12]+dt[k, 13] -
                           dt[k, 6]) / dt[k, 7]  # all features
                    if CF > 1.01:
                        ctQ3 = ctQ3 + 1
                elif feature_id == 5:
                    # percent of Q actual pumped out/20000.
                    CF = abs((dt[k, 5] / dt[k, 7]) * 100)
                ctQ2 = ctQ2 + 1
            else:
                if feature_id <= 5:   # ONLY the Upper Humboldt River
                    CF = -888
                ctQ = ctQ + 1

            if feature_id == 6:
                CF = abs(dt[k, 7])  # Q actual pumped out

            # Assign CF values to MODFLOW cells
            c[int(dt[k, 2]) - 1, int(dt[k, 3]) - 1, int(dt[k, 4]) - 1] = CF

            # Calculate error in water balance
            Err_Bud = abs((dt[k, 5] / dt[k, 7]) * 100)

            # Print outputs to file
            #                fid1.write('%2d %6d %12.1f %12.1f %8.1f %8.1f %8.1f %8.1f %8.2f %8.1f \n' \
            #                           % (ts + 1, k, dt[k, 7], dt[k, 9], dt[k, 10], dt[k, 11], dt[k, 12], dt[k, 13], CF, Err_Bud))
        print 'Nloc not get enough Q and percent ', ctQ, npmploc, (ctQ/npmploc)*100
        print 'Nloc GOT enough Q and percent ', ctQ2, npmploc, (ctQ2 / npmploc) * 100
        print 'Nloc with CF > 1.01 ', ctQ3, npmploc, (ctQ3 / npmploc) * 100
        fid3.write('%3d %3d %6d %8d %5.2d \n' %
                   (feature_id, ts+1, ctQ, npmploc, (ctQ/npmploc)*100))
        fid1.write('Warnings: \n')
        fid1.write('[1] thre_Qact = %d (L3/T) \n' % (thre_Qact))
        fid1.write('[2] Total locations that have Q <= thre_Qact is %d (which is %3.1f %%). \n'
                   % (ctQ, (ctQ / npmploc) * 100))
        fid1.write('[3] These locations were not printed above. \n')


        # # Avoid river and GHB cells
        # fc=np.loadtxt('river_n_GHB.dat')  # location of the feature cells
        # for k in range(fc.shape[0]):
        #     #print dt[k,1],dt[k,2],dt[k,3],dt[k,10]
        #     c[int(fc[k,2])-1,int(fc[k,0])-1,int(fc[k,1])-1]=-999

        # Avoid cells inside a small polygon in the model domain
        # file I08_B118_CA_GroundwaterBasins.shp
        # fc=np.loadtxt('id_of_small_polygon_inside_model_domain.dat')  # location of the feature cells
        # for k in range(fc.shape[0]):
        #     c[int(fc[k,2])-1,int(fc[k,0])-1,int(fc[k,1])-1]=-999

        print 'MIN_C:', c.min()
        print 'MAX_C:', c.max()
        # tmp=np.empty([nrow*ncol,1])
        # tmp1=np.empty([nrow*ncol,nlay])
        #cc = c
        #cc[cc==-999] = 0
        # Print statistic
        # print 'MIN_CC:', cc.min()
        # print 'MAX_CC:', cc.max()
        # for ilay in range(nlay):
        # print ilay


        # use c if no interpolation; use c2 if using interpolation
        # Reshape from 2D nrowxncol to 1D
        tmp1 = np.reshape(c, (nrow*ncol, nlay))
        tmp2 = np.reshape(tmp1, nlay*nrow*ncol)
        np.savetxt('dataset_part2.inp', tmp2, fmt='%8.3f')


        # Generate dataset file
        if ts == 0:
            # Create file dataset_part0.inp
            os.system('del dataset_.tmp')
            fid2 = open('dataset_part0.inp', 'w')
            fid2.write('DATASET \n')
            fid2.write('OBJTYPE "grid3d" \n')
            fid2.write('BEGSCL \n')
            fid2.write('ND 935550 \n')
            fid2.write('NC 935550 \n')
            d_name = feature_name
            fid2.write('NAME "%s" \n' % (d_name))
            fid2.write('TS 1 %d \n' % (ts+1))
            fid2.close()
            os.system(
                'type dataset_part0.inp dataset_part1.inp dataset_part2.inp >> dataset_.tmp')
        elif ts < nts:
            fid2 = open('dataset_part0.inp', 'w')
            fid2.write('TS 1 %d \n' % (ts + 1))
            fid2.close()
            os.system(
                'type dataset_part0.inp dataset_part1.inp dataset_part2.inp >> dataset_.tmp')
        else:
            fid2 = open('dataset_part0.inp', 'w')
            fid2.write('TS 1 %d \n' % (ts + 1))
            fid2.close()
            os.system(
                'type dataset_part0.inp dataset_part1.inp dataset_part2.inp dataset_part3.inp >> dataset_.tmp')

    if feature_id == 1:
        ofile = 'dataset_Humboldt_river_Qthre' + str(thre_Qact) + '.dat'
        copyfile('dataset_.tmp', ofile)
    elif feature_id == 2:
        ofile = 'dataset_tributaries_Qthre' + str(thre_Qact) + '.dat'
        copyfile('dataset_.tmp', ofile)
    elif feature_id == 3:
        ofile = 'dataset_all_rivers_Qthre' + str(thre_Qact) + '.dat'
        copyfile('dataset_.tmp', ofile)
    elif feature_id == 4:
        ofile = 'dataset_all_features_Qthre' + str(thre_Qact) + '.dat'
        copyfile('dataset_.tmp', ofile)
    elif feature_id == 5:
        ofile = 'dataset_Err' + '.dat'
        copyfile('dataset_.tmp', ofile)
    elif feature_id == 6:
        ofile = 'dataset_Qact' + '.dat'
        copyfile('dataset_.tmp', ofile)
fid1.close()
fid3.close()
print 'Done! - Check output file: datase_*.dat'


exit()

# for i in range(nrow):
#    for j in range(ncol):

# print transform_wgs84_to_utm(745850,4322050)
# print transform_utm_to_wgs84(745850,4322050,11)


xmin, ymin, xmax, ymax = [lon.min(), lat.min(), lon.max(), lat.max()]
print 'xmin,ymin,xmax,ymax', xmin, ymin, xmax, ymax

xres = (xmax-xmin)/float(ncol)
yres = (ymax-ymin)/float(nrow)

geotransform = (xmin, xres, 0, ymax, 0, -yres)
# That's (top left x, w-e pixel resolution, rotation (0 if North is up),
#         top left y, rotation (0 if North is up), n-s pixel resolution)
# I don't know why rotation is in twice???

output_raster = gdal.GetDriverByName('GTiff').Create(
    'myraster.tif', ncol, nrow, 1, gdal.GDT_Float32)  # Open the file

output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
# Establish its coordinate encoding
utm_sr = osr.SpatialReference()
utm_sr.ImportFromEPSG(3423)                    #

# utm_sr.ImportFromEPSG(26911)                    # NAD83 / UTM zone 11N

exit()

# Exports the coordinate system
output_raster.SetProjection(utm_sr.ExportToWkt())
# to the file
output_raster.GetRasterBand(1).WriteArray(
    array)   # Writes my array to the raster

# plt.show()


print(utm_sr)  # WKT
# Get the projection name.
print(utm_sr.GetAttrValue('PROJCS'))

# Get the authority.
print(utm_sr.GetAttrValue('AUTHORITY'))
print(utm_sr.GetAttrValue('AUTHORITY', 1))

# Get the datum code.
print(utm_sr.GetAuthorityCode('DATUM'))

# Get the false easting.
print(utm_sr.GetProjParm(osr.SRS_PP_FALSE_EASTING))
