# -*- coding: utf-8 -*-
"""
Testing SpaFHy-input reading from asciigrid-files

Created on Tue Aug 25 10:52:17 2020

@author: 03081268
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import re
from scipy.interpolate import griddata


""" ************ Reading and writing Ascii -grids ********* """

def read_AsciiGrid(fname, setnans=True):

    """ reads AsciiGrid format in fixed format as below:

        ncols         750
        nrows         375
        xllcorner     350000
        yllcorner     6696000
        cellsize      16
        NODATA_value  -9999
        -9999 -9999 -9999 -9999 -9999
        -9999 4.694741 5.537514 4.551162
        -9999 4.759177 5.588773 4.767114
    IN:
        fname - filename (incl. path)
    OUT:
        data - 2D numpy array
        info - 6 first lines as list of strings
        (xloc,yloc) - lower left corner coordinates (tuple)
        cellsize - cellsize (in meters?)
        nodata - value of nodata in 'data'
    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np
    fid = open(fname, 'r')
    info = fid.readlines()[0:6]
    fid.close()

    # print info
    # conversion to float is needed for non-integers read from file...
    if len(info[2]) > 29:
        xloc = float(re.split(' |\.', info[2])[-2])
        yloc = float(re.split(' |\.', info[2])[-2])
        cellsize = float(info[4].split(' ')[-1])
        nodata = float(info[5].split(' ')[-1])
    else:
        xloc = float(info[2].split(' ')[-1])
        yloc = float(info[3].split(' ')[-1])
        cellsize = float(info[4].split(' ')[-1])
        nodata = float(info[5].split(' ')[-1])

    # read rest to 2D numpy array
    data = np.loadtxt(fname, skiprows=6)

    if setnans is True:
        data[data == nodata] = np.NaN
        nodata = np.NaN
    return data, info, (xloc, yloc), cellsize, nodata


def write_AsciiGrid(fname, data, info, fmt='%.18e'):
    """ writes AsciiGrid format txt file
    IN:
        fname - filename
        data - data (numpy array)
        info - info-rows (list, 6rows)
        fmt - output formulation coding

    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np

    # replace nans with nodatavalue according to info
    nodata = int(info[-1].split(' ')[-1])
    data[np.isnan(data)] = nodata
    # write info
    fid = open(fname, 'w')
    fid.writelines(info)
    fid.close()

    # write data
    fid = open(fname, 'a')
    np.savetxt(fid, data, fmt=fmt, delimiter=' ')
    fid.close()

""" ********* Flatten 2d array with nans to dense 1d array ********** """

"""
***** SVE -valuma-alueet -- get gis data to create catchment ******
"""
def create_catchment(pgen, fpath, plotgrids=False, plotdistr=False):
    """
    reads gis-data grids from selected catchments and returns numpy 2d-arrays
    IN:
        ID - SVE catchment ID (int or str)
        fpath - folder (str)
        psoil - soil properties
        plotgrids - True plots
    OUT:
        GisData - dictionary with 2d numpy arrays and some vectors/scalars.

        keys [units]:'dem'[m],'slope'[deg],'soil'[coding 1-4], 'cf'[-],'flowacc'[m2], 'twi'[log m??],
        'vol'[m3/ha],'ba'[m2/ha], 'age'[yrs], 'hc'[m], 'bmroot'[1000kg/ha],'LAI_pine'[m2/m2 one-sided],'LAI_spruce','LAI_decid',
        'info','lat0'[latitude, euref_fin],'lon0'[longitude, euref_fin],loc[outlet coords,euref_fin],'cellsize'[cellwidth,m],
        'peatm','stream','cmask','rockm'[masks, 1=True]

    TODO (6.2.2017 Samuli):
        mVMI-datan koodit >32766 ovat vesialueita ja ei-metsäalueita (tiet, sähkölinjat, puuttomat suot) käytä muita maskeja (maastotietokanta, kysy
        Auralta tie + sähkölinjamaskit) ja IMPOSE LAI ja muut muuttujat ko. alueille. Nyt menevät no-data -luokkaan eikä oteta mukaan laskentaan.
    """
    ID=pgen['catchment_id']
    # fpath = os.path.join(fpath, str(ID) + '\\sve_' + str(ID) + '_')
    #fpath = os.path.join(fpath, str(ID))
    #bname = 'sve_' + str(ID) + '_'
    bname = 'sve_'
    #bname = pgen['prefix'] + '_'+ str(ID) + '_'
    #bname = ' '
    #print os.path.join(fpath, bname + 'dem_16m_aggr.asc')
    # specific leaf area (m2/kg) for converting leaf mass to leaf area
    # SLA = {'pine': 5.54, 'spruce': 5.65, 'decid': 18.46}  # m2/kg, Kellomäki et al. 2001 Atm. Env.
    SLA = {'pine': 6.8, 'spruce': 4.7, 'decid': 14.0}  # Härkönen et al. 2015 BER 20, 181-195

    # values to be set for 'open peatlands' and 'not forest land'
    nofor = {'vol': 0.1, 'ba': 0.01, 'height': 0.1, 'cf': 0.01, 'age': 0.0,
             'LAIpine': 0.01, 'LAIspruce': 0.01, 'LAIdecid': 0.01, 'bmroot': 0.01}
    opeatl = {'vol': 0.01, 'ba': 0.01, 'height': 0.1, 'cf': 0.1, 'age': 0.0,
              'LAIpine': 0.01, 'LAIspruce': 0.01, 'LAIdecid': 0.1, 'bmroot': 0.01}

    # catchment mask cmask ==1, np.NaN outside
    cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'cmask_2023.asc'))
    cmask[np.isfinite(cmask)] = 1.0

    # dem, set values outside boundaries to NaN
    dem, info, pos, cellsize, nodata = read_AsciiGrid(os.path.join(fpath, bname + 'dem_infl_2023.asc'))
    #dem = dem*cmask
    # latitude, longitude arrays
    nrows, ncols = np.shape(dem)
    lon0 = np.arange(pos[0], pos[0] + cellsize*ncols, cellsize)
    lat0 = np.arange(pos[1], pos[1] + cellsize*nrows, cellsize)
    lat0 = np.flipud(lat0)  # why this is needed to get coordinates correct when plotting?

    # flowacc, D-infinity, nr of draining cells
    flowacc, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'flow_accum_d8_2023.asc'))
    #flowacc = flowacc*cmask

    flowpoint, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'flowpoint.asc'))
    #flowpoint = flowpoint*cmask

    # slope, degrees
    slope, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'slope_demfp0_001.txt'))
    #slope = slope*cmask
    # twi
    twi, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'twi.asc'))
    #twi = twi*cmask
    """
    Create soiltype grid and masks for waterbodies, streams, peatlands and rocks
    """
    # Maastotietokanta water bodies: 1=waterbody
    stream, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'virtavesikapea_mtk.asc'))
    stream[np.isfinite(stream)] = -1.0
    #stream = stream*cmask
    # maastotietokanta peatlandmask
    peatm, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'suo_mtk.asc'))
    peatm[np.isfinite(peatm)] = 1.0
    #peatm = peatm*cmask
    # maastotietokanta kalliomaski
    rockm, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'kallioalue_mtk.asc'))
    rockm[np.isfinite(rockm)] = 1.0
    #rockm = rockm*cmask

    peatsoils, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'peatsoils.asc'))
    #peatsoils[np.isfinite(peatsoils)] = 1.0
    #print(np.where(peatsoils != 1.0))

    """
    gtk soilmap: read and re-classify into 4 texture classes
    #GTK-pintamaalaji grouped to 4 classes (Samuli Launiainen, Jan 7, 2017)
    #Codes based on maalaji 1:20 000 AND ADD HERE ALSO 1:200 000
    """
    CoarseTextured = [195213, 195314, 19531421, 195313, 195310]
    MediumTextured = [195315, 19531521, 195215, 195214, 195601, 195411, 195112,
                      195311, 195113, 195111, 195210, 195110, 195312]
    FineTextured = [19531521, 195412, 19541221, 195511, 195413, 195410,
                    19541321, 195618]
    Peats = [195512, 195513, 195514, 19551822, 19551891, 19551892]
    Water = [195603]
    Peatsoils = [0, 35411]

    #CoarseTextured = [195213, 195314, 19531421, 195313, 195310]
    #MediumTextured = [195315, 19531521, 195215, 195214, 195601, 195411, 195112,
    #                  195311, 195113, 195111, 195210, 195110, 195312]
    #FineTextured = [19531521, 195412, 19541221, 195511, 195413, 195410,
    #                19541321, 195618]
    #Peats = [195512, 195513, 195514, 19551822, 19551891, 19551892]
    #Water = [195603]

    gtk_s, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'gtk_pintamaa_2023_16.asc'))

    #print(np.unique(gtk_s)
    r, c = np.shape(gtk_s)
    soil = np.ravel(gtk_s)

    #del gtk_s

    soil[np.in1d(soil, CoarseTextured)] = 1.0  # ; soil[f]=1; del f
    soil[np.in1d(soil, MediumTextured)] = 2.0
    soil[np.in1d(soil, FineTextured)] = 3.0
    soil[np.in1d(soil, Peats)] = 4.0
    soil[np.in1d(soil, Water)] = -1.0
    print('soil uniques', np.unique(soil.flatten()))
    #soil[soil == -1.0] = 2.0
    soil[(np.where(soil == 4.0)) and (np.where(peatsoils.flatten() != Peatsoils))] = 2.0
    #print(np.unique(soil)
    print('soil uniques', np.unique(soil.flatten()))

    # reshape back to original grid
    #pix = np.where(peatsoils.flatten() == 1.0)
    #print(pix)
    soil[np.in1d(peatsoils, Peatsoils)] = 4.0
    #soil[soil == pix] = 4.0
    #soil[np.in1d(peatsoils, 1.0)] = 4.0
    soil = soil.reshape(r, c)
    del r, c
    #plt.imshow(soil);plt.colorbar()
    #soil[np.in1d(soil, 4.0) and ~np.in1d(peatsoils, 1.0)] = 2.0
    #soil[(np.where(peatsoils != 1.0)) and (np.where(soil == 4.0))] = 2.0
    #plt.imshow(soil);plt.colorbar()

    # update waterbody mask
    stream[np.where(soil == -1.0)] = -1.0
    # update lakes from soil to mineral, will be masked later
    soil[np.where(soil == -1.0)] = 2.0
    
    # update cmask (lakes)
    #cmask[np.where(soil == -1.0)] = np.nan
    #stream = stream*cmask

    # virtual ditches
    #stream[flowacc > 60000] = -1.0


    # update catchment mask so that water bodies are left out (SL 20.2.18)
    #cmask[soil == -1.0] = np.NaN

    #Warn, remove this
    #cmask[soil <= 0] = np.NaN
    #soil = soil * cmask
    plt.imshow(soil);plt.colorbar()
    #plt.imshow(peatsoils*cmask)

    """ stand data (MNFI)"""
    # stand volume [m3ha-1]
    vol, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'tilavuus.asc'), setnans=False)
    #vol = vol*cmask

    # indexes for cells not recognized in mNFI
    ix_n = np.where((vol >= 32727) | (vol == -9999) )  # no satellite cover or not forest land: assign arbitrary values
    ix_p = np.where((vol >= 32727) & (peatm == 1))  # open peatlands: assign arbitrary values
    ix_w = np.where((vol >= 32727) & (stream == -1))  # waterbodies: leave out
    cmask[ix_w] = np.NaN  # NOTE: leaves waterbodies out of catchment mask

    stream[~np.isfinite(stream)] = 0.0
    #print(stream)

    vol[ix_n] = nofor['vol']
    vol[ix_p] = opeatl['vol']
    vol[ix_w] = np.NaN

    #pine volume [m3 ha-1]
    p_vol, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname +'manty.asc'))
    #p_vol = p_vol*cmask
    #spruce volume [m3 ha-1]
    s_vol, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname +'kuusi.asc'))
    #s_vol = s_vol*cmask
    #birch volume [m3 ha-1]
    b_vol, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname +'koivu.asc'))
    #b_vol = b_vol*cmask


    # basal area [m2 ha-1]
    ba, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'ppa.asc') )
    ba[ix_n] = nofor['ba']
    ba[ix_p] = opeatl['ba']
    ba[ix_w] = np.NaN
    #ba = ba*cmask

    # tree height [m]
    height, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'keskipituus.asc'))
    height = 0.1*height  # m
    height[ix_n] = nofor['height']
    height[ix_p] = opeatl['height']
    height[ix_w] = np.NaN
    #height = height*cmask

    # canopy closure [-]
    cf, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'latvuspeitto.asc'))
    cf = 1e-2*cf
    cf[ix_n] = nofor['cf']
    cf[ix_p] = opeatl['cf']
    cf[ix_w] = np.NaN
    #cf = cf*cmask
    # cfd, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'lehtip_latvuspeitto.asc'))
    # cfd = 1e-2*cfd  # percent to fraction

    # stand age [yrs]
    age, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname+'ika.asc'))
    age[ix_n] = nofor['age']
    age[ix_p] = opeatl['age']
    age[ix_w] = np.NaN
    #age = age*cmask

    # leaf biomasses and one-sided LAI
    bmleaf_pine, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_manty_neulaset.asc'))
    #bmleaf_pine = bmleaf_pine*cmask
    bmleaf_spruce, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_kuusi_neulaset.asc'))
    #bmleaf_spruce = bmleaf_spruce*cmask
    bmleaf_decid, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_lehtip_neulaset.asc'))
    #bmleaf_decid = bmleaf_decid*cmask
    bmleaf_pine[ix_n]=np.NaN; bmleaf_spruce[ix_n]=np.NaN; bmleaf_decid[ix_n]=np.NaN;


    LAI_pine = 1e-3*bmleaf_pine*SLA['pine']  # 1e-3 converts 10kg/ha to kg/m2
    LAI_pine[ix_n] = nofor['LAIpine']
    LAI_pine[ix_p] = opeatl['LAIpine']
    LAI_pine[ix_w] = np.NaN

    LAI_spruce = 1e-3*bmleaf_spruce*SLA['spruce']
    LAI_spruce[ix_n] = nofor['LAIspruce']
    LAI_spruce[ix_p] = opeatl['LAIspruce']
    LAI_spruce[ix_w] = np.NaN

    LAI_decid = 1e-3*bmleaf_decid*SLA['decid']
    LAI_decid[ix_n] = nofor['LAIdecid']
    LAI_decid[ix_p] = opeatl['LAIdecid']
    LAI_decid[ix_w] = np.NaN

    LAI_decid_shrub, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'lai_decid_shrub.asc'))
    #LAI_decid_shrub = LAI_decid_shrub*cmask
    LAI_decid_shrub[ix_w] = np.NaN

    LAI_everg_shrub, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'lai_everg_shrub.asc'))
    #LAI_everg_shrub = LAI_everg_shrub*cmask
    LAI_everg_shrub[ix_w] = np.NaN
    LAI_shrub = LAI_decid_shrub + LAI_everg_shrub
    LAI_shrub[LAI_shrub < 0] = 0.001

    LAI_grass, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'lai_graminoids.asc'))
    #LAI_grass = LAI_grass*cmask
    LAI_grass[ix_w] = np.NaN
    LAI_grass[LAI_grass < 0] = 0.001


    bmroot_pine, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_manty_juuret.asc'))
    bmroot_spruce, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_kuusi_juuret.asc'))
    bmroot_decid, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'bm_lehtip_juuret.asc'))
    bmroot = 1e-2*(bmroot_pine + bmroot_spruce + bmroot_decid)  # 1000 kg/ha
    bmroot[ix_n] = nofor['bmroot']
    bmroot[ix_p] = opeatl['bmroot']
    bmroot[ix_w] = np.NaN
    #bmroot = bmroot*cmask
    # site types
    maintype, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'paatyyppi.asc'))
    #maintype = maintype*cmask

    # interpolating maintype to not have nan on roads
    x = np.arange(0, maintype.shape[1])
    y = np.arange(0, maintype.shape[0])
    #mask invalid values
    array = np.ma.masked_invalid(maintype)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]

    maintype = griddata((x1, y1), newarr.ravel(),
                          (xx, yy),
                             method='nearest')
    #maintype = maintype*cmask


    sitetype, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname + 'kasvupaikka.asc'))
    #sitetype = sitetype*cmask
#    #HERE DUPLICATE, REMOVE LATER
#    #site main class: 1- mineral soil; 2-spruce fen; 3-pine bog; 4-open peatland
#    smc, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname +'paatyyppi.asc'))
#    smc = smc*cmask
#    #site fertility class: 1- ... 8 (1 Lehdot, 2 OMT, 3 MT, 4...)
#    sfc, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, bname +'kasvupaikka.asc'))
#    sfc = sfc*cmask
#    # downslope distance
#    print ('Warning: dsdistance decativated!')
#    #dsdist, _, _, _, _  = read_AsciiGrid(os.path.join(fpath, bname +'downslope_distance.asc'))
#    dsdist=None



    # catchment outlet location and catchment mean elevation
    (iy, ix) = np.where(flowacc == np.nanmax(flowacc))
    loc = {'lat': lat0[iy], 'lon': lon0[ix], 'elev': np.nanmean(dem)}
    # dict of all rasters
    GisData = {'cmask': cmask, 'dem': dem, 'flowacc': flowacc, 'flowpoint': flowpoint, 'slope': slope,
               'twi': twi, 'gtk_soilcode': gtk_s, 'soilclass': soil, 'peatm': peatm, 'stream': stream,
               'rockm': rockm, 'LAI_pine': LAI_pine, 'LAI_spruce': LAI_spruce,
               'LAI_conif': LAI_pine + LAI_spruce,'soil':soil,
               'LAI_decid': LAI_decid, 'LAI_shrub': LAI_shrub, 'LAI_grass': LAI_grass, 'bmroot': bmroot, 'ba': ba, 'hc': height,
               'vol': vol, 'p_vol':p_vol,'s_vol':s_vol,'b_vol':b_vol,'cf': cf, 'age': age, 'maintype': maintype, 'sitetype': sitetype,
               'cellsize': cellsize, 'info': info, 'lat0': lat0, 'lon0': lon0, 'loc': loc}

    if plotgrids is True:
        # %matplotlib qt
        # xx, yy = np.meshgrid(lon0, lat0)
        plt.close('all')

        plt.figure()

        plt.subplot(221)
        plt.imshow(dem); plt.colorbar(); plt.title('DEM (m)')
        plt.plot(ix, iy,'rs')
        plt.subplot(222)
        plt.imshow(twi); plt.colorbar(); plt.title('TWI')
        plt.subplot(223)
        plt.imshow(slope); plt.colorbar(); plt.title('slope(deg)')
        plt.subplot(224)
        plt.imshow(flowacc); plt.colorbar(); plt.title('flowacc (m2)')

        plt.figure(figsize=(6, 14))

        plt.subplot(221)
        plt.imshow(soil); plt.colorbar(); plt.title('soiltype')
        mask = cmask.copy()*0.0
        mask[np.isfinite(peatm)] = 1
        mask[np.isfinite(rockm)] = 2
        mask[np.isfinite(stream)] = 3

        plt.subplot(222)
        plt.imshow(mask); plt.colorbar(); plt.title('masks')
        plt.subplot(223)
        plt.imshow(LAI_pine+LAI_spruce + LAI_decid); plt.colorbar(); plt.title('LAI (m2/m2)')
        plt.subplot(224)
        plt.imshow(cf); plt.colorbar(); plt.title('cf (-)')


        plt.figure(figsize=(6,11))
        plt.subplot(321)
        plt.imshow(vol); plt.colorbar(); plt.title('vol (m3/ha)')
        plt.subplot(323)
        plt.imshow(height); plt.colorbar(); plt.title('hc (m)')
        #plt.subplot(223)
        #plt.imshow(ba); plt.colorbar(); plt.title('ba (m2/ha)')
        plt.subplot(325)
        plt.imshow(age); plt.colorbar(); plt.title('age (yr)')
        plt.subplot(322)
        plt.imshow(1e-3*bmleaf_pine); plt.colorbar(); plt.title('pine needles (kg/m2)')
        plt.subplot(324)
        plt.imshow(1e-3*bmleaf_spruce); plt.colorbar(); plt.title('spruce needles (kg/m2)')
        plt.subplot(326)
        plt.imshow(1e-3*bmleaf_decid); plt.colorbar(); plt.title('decid. leaves (kg/m2)')

    if plotdistr is True:
        twi0 = twi[np.isfinite(twi)]
        vol = vol[np.isfinite(vol)]
        lai = LAI_pine + LAI_spruce + LAI_decid
        lai = lai[np.isfinite(lai)]
        soil0 = soil[np.isfinite(soil)]

        plt.figure(100)
        plt.subplot(221)
        plt.hist(twi0, bins=100, color='b', alpha=0.5)
        plt.ylabel('f');plt.ylabel('twi')

        s = np.unique(soil0)
        colcode = 'rgcym'
        for k in range(0,len(s)):
            # print k
            a = twi[np.where(soil==s[k])]
            a = a[np.isfinite(a)]
            plt.hist(a, bins=50, alpha=0.5, color=colcode[k], label='soil ' +str(s[k]))
        plt.legend()
        plt.show()

        plt.subplot(222)
        plt.hist(vol, bins=100, color='k'); plt.ylabel('f'); plt.ylabel('vol')
        plt.subplot(223)
        plt.hist(lai, bins=100, color='g'); plt.ylabel('f'); plt.ylabel('lai')
        plt.subplot(224)
        plt.hist(soil0, bins=5, color='r'); plt.ylabel('f');plt.ylabel('soiltype')

    return GisData





