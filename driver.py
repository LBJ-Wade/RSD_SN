import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import quad
from scipy.interpolate import interp1d
import os
import sys

if len(sys.argv)!=3 :
    print "Usage: driver.py model (0,1) nsample"
    exit(1)
mod_int=int(sys.argv[1])
nsample_sz=float(sys.argv[2])

if mod_int==0 :
    model='2MPZ'
    zmax=0.2
    sigmaz=0.013
    fsky=0.65
    ell_max=100
elif mod_int==1 :
    model='WISExSCOS'
    zmax=0.5
    sigmaz=0.033
    fsky=0.65
    ell_max=300
else :
    print "Unknown model %d"%mod_int
    exit(1)

def nz_f(z) :
    #N(z) in deg^-2
    if mod_int==0 :
        z0=0.053
        a=561.4
        b=2.21
        c=1.43
    if mod_int==1 :
        z0=0.26070073
        a=5735.40881666
        b=1.08776304
        c=3.00661052
    return a*(z/z0)**b*np.exp(-(z/z0)**c)

def pdf_photo(z,z0,zf,sz) :
    denom=1./np.sqrt(2*sz*sz)
    return 0.5*(erf((zf-z)*denom)-erf((z0-z)*denom))

def get_zarr(s_photoz,n_width,z_max) :
    dz_half=s_photoz*n_width*0.5
    zc_arr=[dz_half/(1-dz_half)]
    zw_arr=[dz_half*(1+zc_arr[0])]
    z_last=zc_arr[0]
    while z_last<=z_max :
        zc_new=(z_last+dz_half*(2+z_last))/(1-dz_half)
        zw_new=dz_half*(1+zc_new)
        zc_arr.append(zc_new)
        zw_arr.append(zw_new)
        z_last=zc_new
    zc_arr=np.array(zc_arr)
    zw_arr=np.array(zw_arr)
    indices=np.where(zc_arr<z_max)
    zc_arr_out=zc_arr[indices]
    zw_arr_out=zw_arr[indices]
    sz_arr_out=s_photoz*(1+zc_arr_out)

    return zc_arr_out,zw_arr_out,sz_arr_out

numz=512
zarr=zmax*np.arange(numz)/(numz-1.)
nzarr=nz_f(zarr)

zca,zwa,sza=get_zarr(sigmaz,nsample_sz,zmax-nsample_sz*sigmaz*(1+zmax))
n_bins=len(zca)
nz_bins=[]
for zc,zw,sz in zip(zca,zwa,sza) :
    nz_bins.append(pdf_photo(zarr,zc-zw,zc+zw,sz)*nzarr)

nzt=np.zeros_like(nzarr)
n=0
for nz in nz_bins :
    nzf=interp1d(zarr,nz,bounds_error=False,fill_value=0);
    n_here=quad(nzf,0,0.6)[0]
    n+=n_here
    print n_here
    nzt+=nz
    plt.plot(zarr,nz);
print n
plt.plot(zarr,nzarr);
plt.plot(zarr,nzt);
plt.show()

def limberjack_driver(fname_nz1,fname_nz2,fname_bias,prefix_out,w_rsd=False,lmax=300,lmin_limber=0) :
    stout="omega_m= 0.3\n"
    stout+="omega_l= 0.7\n"
    stout+="omega_b= 0.05\n"
    stout+="w0= -1.0\n"
    stout+="wa= 0.0\n"
    stout+="h= 0.67\n"
    stout+="ns= 0.96\n"
    stout+="s8= 0.83\n"
    stout+="l_limber_min= %d\n"%lmin_limber
    stout+="d_chi= 3.\n"
    stout+="z_kappa= 10.\n"
    stout+="z_isw= 10.\n"
    stout+="r_smooth= 0.0001\n"
    stout+="l_max= %d\n"%lmax
    stout+="do_nc= 1\n"
    stout+="has_nc_dens= 1\n"
    if w_rsd :
        stout+="has_nc_rsd= 1\n"
    else :
        stout+="has_nc_rsd= 0\n"
    stout+="has_nc_lensing= 0\n"
    stout+="do_shear= 0\n"
    stout+="has_sh_intrinsic= 0\n"
    stout+="do_cmblens= 0\n"
    stout+="do_isw= 0\n"
    stout+="do_w_theta= 0\n"
    stout+="use_logbin= 0\n"
    stout+="theta_min= 0\n"
    stout+="theta_max= 0\n"
    stout+="n_bins_theta= 0\n"
    stout+="n_bins_decade= 0\n"
    stout+="window_1_fname= "+fname_nz1+"\n"
    stout+="window_2_fname= "+fname_nz2+"\n"
    stout+="bias_fname= "+fname_bias+"\n"
    stout+="sbias_fname= bullcrap\n"
    stout+="abias_fname= bullcrap\n"
    stout+="pk_fname= Pk_CAMB_test.dat\n"
    stout+="prefix_out= "+prefix_out+"_lj\n"
    f=open(prefix_out+"_lj_param.ini","w")
    f.write(stout)
    f.close()
    os.system('./LimberJack_dir/LimberJack '+prefix_out+'_lj_param.ini >> '+prefix_out+'_log')
    l,cl=np.loadtxt(prefix_out+"_lj_cl_dd.txt",unpack=True)
    os.system('rm '+prefix_out+'_lj*')
    return cl

def get_cls(z,nzs,lmax=300,lmin_limber=0,prefix='out',w_rsd=False) :
    if os.path.isfile(prefix+"_cls.npy") :
        return np.load(prefix+"_cls.npy")
    else :
        nbins=len(nzs)
        cls=np.zeros([lmax+1,nbins,nbins])
        np.savetxt(prefix+"bias.txt",np.transpose([z,1+z]))
        for i1 in np.arange(nbins) :
            np.savetxt(prefix+"nz1.txt",np.transpose([z,nzs[i1]]))
            for i2 in np.arange(nbins-i1)+i1 :
                print i1,i2
                np.savetxt(prefix+"nz2.txt",np.transpose([z,nzs[i2]]))
                cls[:,i1,i2]=limberjack_driver(prefix+"nz1.txt",prefix+"nz2.txt",prefix+"bias.txt",prefix,
                                               w_rsd=w_rsd,lmax=lmax,lmin_limber=lmin_limber)
                if i1!=i2 :
                    cls[:,i2,i1]=cls[:,i1,i2]

        np.save(prefix+"_cls",cls)
        os.system('rm '+prefix+'bias.txt '+prefix+'nz1.txt '+prefix+'nz2.txt')
        return cls

cls_norsd=get_cls(zarr,nz_bins,prefix='data/'+model+'_ns%.1lf_norsd'%nsample_sz,w_rsd=False,lmax=ell_max)
cls_wrsd=get_cls(zarr,nz_bins,prefix='data/'+model+'_ns%.1lf_wrsd'%nsample_sz,w_rsd=False,lmax=ell_max)
print np.shape(cls_norsd)
print np.shape(cls_wrsd)


'''
cls_limber_norsd=get_cls(zarr,nz_bins,lmax=300,prefix='no_rsd_limber',w_rsd=False)
cls_limber_wrsd=get_cls(zarr,nz_bins,lmax=300,prefix='w_rsd_limber',w_rsd=True)
cls_nolimber_norsd=get_cls(zarr,nz_bins,lmax=300,prefix='no_rsd_nolimber',w_rsd=False)
cls_nolimber_wrsd=get_cls(zarr,nz_bins,lmax=300,prefix='w_rsd_nolimber',w_rsd=True)

ls=np.arange(301)
for i1 in np.arange(n_bins) :
    for i2 in np.arange(n_bins-i1)+i1 :
        plt.figure()
        plt.title("%d "%i1+"%d "%i2)
        plt.plot(ls,cls_limber_norsd[:,i1,i2],'r--')
        plt.plot(ls,cls_limber_wrsd[:,i1,i2],'g--')
        plt.plot(ls,cls_nolimber_norsd[:,i1,i2],'r-')
        plt.plot(ls,cls_nolimber_wrsd[:,i1,i2],'g-')
        plt.loglog()
plt.show()
'''
