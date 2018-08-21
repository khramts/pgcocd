import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# Load eigenvec file and fam file
dat=np.genfromtxt('IOEU.pca.noLD.maf0.05.pruned.eigenvec', dtype=None)
fam=np.genfromtxt('IOEU.pca.noLD.maf0.05.pruned.fam', dtype=None)
sex=fam['f4']
stat=fam['f5']


### Run these functions definitions without changing anything
def plot_dens(m1,m2,filt_m,m1_x,m2_x,filt_mx, name1='PC1', name2='PC2'):
    """

    :param m1: array of first PC for controls (can be any PC that is being compared)
    :param m2: array of second PC for controls(can be any PC that is being compared)
    :param filt_m:
    :param m1_x: array of first PC for cases (can be any PC that is being compared)
    :param m2_x: array of second PC for controls(can be any PC that is being compared)
    :param filt_mx:
    :param name1: name of first PC
    :param name2: name of second PC
    :return:
    """
    # creating grid 100 by 100
    xmin = np.min([m1.min(),m1_x.min()])
    xmax = np.max([m1.max(),m1_x.max()])
    ymin = np.min([m2.min(),m2_x.min()])
    ymax = np.max([m2.max(),m2_x.max()])
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    # done creating grid
    # create vector with all pixel positions in a grid
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1[filt_m], m2[filt_m]])
    # Scott's rule used to automatically decide kernel size
    kernel = stats.gaussian_kde(values)
    # Z density matrix (100x100)
    Z = np.reshape(kernel(positions).T, X.shape)
    #print(Z.sum())
    values = np.vstack([m1_x[filt_mx], m2_x[filt_mx]])
    kernel = stats.gaussian_kde(values)
    Zx = np.reshape(kernel(positions).T, X.shape)
    #print(Zx.sum())
    plt.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.subplot(131)
    plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
               extent=[xmin, xmax, ymin, ymax])
    plt.plot(m1[filt_m], m2[filt_m], 'k.', markersize=2,alpha=0.2)
    plt.plot(m1[filt_m==False], m2[filt_m==False], 'rx', markersize=8,alpha=0.75)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(name1)
    plt.ylabel(name2)
    #plt.title('Female Controls')
    plt.title('Controls')
    plt.subplot(132)
    plt.imshow(np.rot90(Zx), cmap=plt.cm.gist_earth_r,
               extent=[xmin, xmax, ymin, ymax])
    plt.plot(m1_x[filt_mx], m2_x[filt_mx], 'k.', markersize=2,alpha=0.2)
    plt.plot(m1_x[filt_mx==False], m2_x[filt_mx==False], 'rx', markersize=8,alpha=0.75)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(name1)
    #plt.title('Female Cases')
    plt.title('Cases')
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

### Run these functions definitions without changing anything
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

### Run these functions definitions without changing anything
def remove_by_threshold(m1,m2,filt_m,m1_x,m2_x,filt_mx,thresh=0.3):
    i = 0
    while (remove_one(m1,m2,filt_m,m1_x,m2_x,filt_mx, thresh) > thresh):
        # print i
        i=i+1
    while remove_one(m1,m2,filt_m,m1_x,m2_x,filt_mx, thresh, False) < -thresh:
        # print i
        i=i+1

### Here you can specify which PC to use. E.g. dat['f2'] represents PC1, dat['f3'] represents PC2 and so on.
m1=dat['f2'][(stat==1)]
m2=dat['f3'][(stat==1)]
m3=dat['f4'][(stat==1)]
m4=dat['f5'][(stat==1)]
fam_copy=fam[(stat==1)]
m1_x=dat['f2'][(stat==2)]
m2_x=dat['f3'][(stat==2)]
m3_x=dat['f4'][(stat==2)]
m4_x=dat['f5'][(stat==2)]
fam_copyx=fam[(stat==2)]

filt_m = np.ones(len(m1))==1
filt_mx = np.ones(len(m1_x))==1

# set thresholds
tresh = [1.0, 0.99, 0.95, 0.90, 0.85, 0.8, 0.75, 0.7]

## This code makes PC plots and a list of individuals retained at each threshold
for tr in tresh:
    remove_by_threshold(m3, m4, filt_m, m3_x, m4_x, filt_mx, thresh=tr)
    remove_by_threshold(m1,m2,filt_m,m1_x,m2_x,filt_mx,thresh=tr)
    remove_by_threshold(m2,m3,filt_m,m2_x,m3_x,filt_mx,thresh=tr)
    print filt_m.sum(), filt_mx.sum()
    plot_dens(m1,m2,filt_m,m1_x,m2_x,filt_mx)
    plt.savefig('PC_12_'+str(tr)+'.png',dpi=200)
    plot_dens(m2,m3,filt_m,m2_x,m3_x,filt_mx,'PC2','PC3')
    plt.savefig('PC_23_'+str(tr)+'.png',dpi=200)
    plot_dens(m3,m4,filt_m,m3_x,m4_x,filt_mx,'PC3','PC4')
    plt.savefig('PC_34_'+str(tr)+'.png',dpi=200)
    plt.close('all')
    np.savetxt('test_'+str(tr)+'.txt', np.hstack([fam_copy[filt_m][['f0','f1']], fam_copyx[filt_mx][['f0','f1']]]), fmt=('%s', '%s'))
