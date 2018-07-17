import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


dat=np.genfromtxt('OCGAS_IOCDFGC_EU_info0.6cert0.8.maf0.05.geno0.02.MHC.inv.removed.pruned2.outlier_removed.evec', dtype=None, skip_header=1)

fam=np.genfromtxt('1_2_test_1.0.txt.fam', dtype=None)

sex=fam['f4']
stat=fam['f5']

from scipy import stats

import matplotlib.pyplot as plt

def plot_dens(m1,m2,filt_m,m1_x,m2_x,filt_mx, name1='PC1', name2='PC2'):
    xmin = np.min([m1.min(),m1_x.min()])
    xmax = np.max([m1.max(),m1_x.max()])
    ymin = np.min([m2.min(),m2_x.min()])
    ymax = np.max([m2.max(),m2_x.max()])
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1[filt_m], m2[filt_m]])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    values = np.vstack([m1_x[filt_mx], m2_x[filt_mx]])
    kernel = stats.gaussian_kde(values)
    Zx = np.reshape(kernel(positions).T, X.shape)
    plt.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.subplot(131)
    plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
               extent=[xmin, xmax, ymin, ymax])
    plt.plot(m1[filt_m], m2[filt_m], 'k.', markersize=2,alpha=0.2)
    plt.plot(m1[filt_m==False], m2[filt_m==False], 'rx', markersize=5,alpha=0.25)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.title('Female Controls')
    plt.subplot(132)
    plt.imshow(np.rot90(Zx), cmap=plt.cm.gist_earth_r,
               extent=[xmin, xmax, ymin, ymax])
    plt.plot(m1_x[filt_mx], m2_x[filt_mx], 'k.', markersize=2,alpha=0.2)
    plt.plot(m1_x[filt_mx==False], m2_x[filt_mx==False], 'rx', markersize=5,alpha=0.25)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(name1)
    plt.title('Female Cases')
    plt.subplot(133)
    temp = Z*np.nan
    temp[(Z+Zx)>20.0] = 1
    plt.imshow(np.rot90((Z-Zx)/(Z+Zx)*temp), cmap=plt.cm.RdYlBu,
               extent=[xmin, xmax, ymin, ymax], vmin = -1, vmax = 1)
    # plt.contour(((Z-Zx)/(Z+Zx)*temp).T, levels=[-0.5,0.5], extent=[xmin, xmax, ymin, ymax], colors='k')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(name1)
    plt.title('Densities Ratio')
    plt.tight_layout()
    # plt.show()

def remove_one(m1,m2,filt_m,m1_x,m2_x,filt_mx,thresh,from_1=True):
    values = np.vstack([m1[filt_m], m2[filt_m]])
    kernel = stats.gaussian_kde(values)
    valuesx = np.vstack([m1_x[filt_mx], m2_x[filt_mx]])
    kernelx = stats.gaussian_kde(valuesx)
    Z = (kernel(values) - kernelx(values)) / (kernel(values) + kernelx(values))
    Zx = (kernel(valuesx) - kernelx(valuesx)) / (kernel(valuesx) + kernelx(valuesx))
    # if np.max(np.abs(Z))>np.max(np.abs(Zx)):
    if from_1:
        temp = np.where(Z == np.max(Z))[0][0]
        if (Z[temp] > 0) & (Z[temp] > thresh):
            filt_m[ np.where(filt_m)[0][temp] ] = False
            # print np.where(filt_m)[0][temp]
        return Z[temp]
    else:
        temp = np.where(Zx == np.min(Zx))[0][0]
        # print Zx[temp]
        if (Zx[temp] < 0) & (np.abs(Zx[temp]) > thresh):
            filt_mx[ np.where(filt_mx)[0][temp] ] = False
            # print np.where(filt_mx)[0][temp]
        return Zx[temp]




def remove_by_threshold(m1,m2,filt_m,m1_x,m2_x,filt_mx,thresh=0.3):
    i = 0
    while (remove_one(m1,m2,filt_m,m1_x,m2_x,filt_mx, thresh) > thresh):
        # print i
        i=i+1
    while remove_one(m1,m2,filt_m,m1_x,m2_x,filt_mx, thresh, False) < -thresh:
        # print i
        i=i+1

# DONT FORGET TO CHANGE NAMES IN THE PLOTS. FUNCTION ABOVE.
sexi=2

m1=dat['f1'][(sex==sexi) & (stat==1)]
m2=dat['f2'][(sex==sexi) & (stat==1)]
m3=dat['f3'][(sex==sexi) & (stat==1)]
m4=dat['f4'][(sex==sexi) & (stat==1)]
fam_copy=fam[(sex==sexi) & (stat==1)]
m1_x=dat['f1'][(sex==sexi) & (stat==2)]
m2_x=dat['f2'][(sex==sexi) & (stat==2)]
m3_x=dat['f3'][(sex==sexi) & (stat==2)]
m4_x=dat['f4'][(sex==sexi) & (stat==2)]
fam_copyx=fam[(sex==sexi) & (stat==2)]

filt_m = np.ones(len(m1))==1
filt_mx = np.ones(len(m1_x))==1

#tresh = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
tresh = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
#tresh = [0.4, 0.3, 0.2, 0.1]

for tr in tresh:
    remove_by_threshold(m1,m2,filt_m,m1_x,m2_x,filt_mx,thresh=tr)
    remove_by_threshold(m2,m3,filt_m,m2_x,m3_x,filt_mx,thresh=tr)
    remove_by_threshold(m3,m4,filt_m,m3_x,m4_x,filt_mx,thresh=tr)
    print filt_m.sum(), filt_mx.sum()
    plot_dens(m1,m2,filt_m,m1_x,m2_x,filt_mx)
    plt.savefig(str(sexi)+'_PC_12_'+str(tr)+'.png',dpi=200)
    plot_dens(m2,m3,filt_m,m2_x,m3_x,filt_mx,'PC2','PC3')
    plt.savefig(str(sexi)+'_PC_23_'+str(tr)+'.png',dpi=200)
    plot_dens(m3,m4,filt_m,m3_x,m4_x,filt_mx,'PC3','PC4')
    plt.savefig(str(sexi)+'_PC_34_'+str(tr)+'.png',dpi=200)
    plt.close('all')
    np.savetxt(str(sexi)+'_test_'+str(tr)+'.txt', np.hstack([fam_copy[filt_m][['f0','f1']], fam_copyx[filt_mx][['f0','f1']]]), fmt=('%s', '%s'))