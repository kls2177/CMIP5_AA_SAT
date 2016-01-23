import numpy as np

class averages(object):

    def __init__(self,area=None):
        self.area = area

    def AAavg(self, var='var', lat='lat', lon='lon', mask='AA_mask'):
        """Returns the area-averged variable for a certain masked region"""
        y = lat*np.pi/180
        coslat = np.cos(y)
        coslat = np.tile(coslat,(lon.size,1)).T

        nan = float('nan')
        AA_tmp = np.zeros([var.shape[0], lat.size, lon.size])
        coslat_tmp = np.zeros([lat.size, lon.size])
        for i in range(lat.size):
            for j in range(lon.size):
                if mask[i,j] == 0:
                    AA_tmp[:,i,j] = var[:,i,j]
                    coslat_tmp[i,j] = coslat[i,j]
                else:
                    AA_tmp[:,i,j] = nan
                    coslat_tmp[i,j] = nan

        #area weighting (weight by cos(lat) because area of grid boxes get smaller closer to the pole)
        AA_tmp_cos = np.zeros([var.shape[0], lat.size, lon.size])
        for i in range(lat.size):
            for j in range(lon.size):
                AA_tmp_cos[:,i,j] = AA_tmp[:,i,j] * coslat_tmp[i,j]

        #sum over all gridpoints
        cos_sum1 = np.nansum(coslat_tmp, axis=0)
        cos_sum = np.nansum(cos_sum1);

        AA_sum1 = np.nansum(AA_tmp_cos, axis=1)
        AA_sum = np.nansum(AA_sum1,axis=1)

        #normalize by area
        AAavg = AA_sum / cos_sum

        return AAavg
