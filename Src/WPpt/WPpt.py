#!/usr/bin/env python

import numpy as np

import os
import sys
import time

import kappa_L
import kappa_accum
import kappa_omega
import fitting_models
 
#----------------------------------------title-----------------------------------------

print('------------------------------------------------------------------------')
print('----------***----******---***---*******---*******--********-------------')
print('-----------***---******---***--***---**--***---**-----**----------------')
print('-----------***--***-***--***---********--********-----**----------------')
print('------------***-***--***-***---***-------***----------**----------------')
print('------------*******--******----***-------***----------**----------------')
print('------------*******--******----**--------**-----------**---V.3.1--------')
print('------------------------------------------------------------------------')
print('---------The WPpt module for coherence transport in WPPT package--------')
print('------------------------------python code-------------------------------')
print('---------------------------By Zhongwei Zhang----------------------------')
print('---------------------Email: zhongwei@tongji.edu.cn----------------------')
print('---------------------------Date: Dec. 20 2021---------------------------')
print('------------------------------------------------------------------------')
               

#------------------------------ setting preparation --------------------------------
               
print('######## I. Commands are reading from CONTROL.WPPT...')
print('------------------------------------------------------------------------')
print(' ')

filename='CONTROL.WPPT'

file_have=os.path.exists(filename)

if file_have==False:
    
    print('ERROR: The file of CONTROL.WPPT does not exist...')
    sys.exit()
    
else:
                    
    f=open(filename,'r')
    fm=f.readlines()
    f.close()
    
    Nbranch=int(fm[3].split()[1])*3
    Nkpoints=int(fm[4].split()[1])
    
    calculation_mode= fm[29].split()[1]  # 1: only times calculations
                                 # 2: times+kappa calculations
                                 # 3: kappa calculations
    
    times_mode= fm[30].split()[1]    # 1: conventional lifetimes from SED and coherence times
                                 #   from coherence distribution function (wavelet transform)
                                 # 2: lifetimes and coherence times from phonon decay
                                 # 3: conventional lifetimes (SED) and coherence times from
                                 #   from spectroscopy or spectral energy
    
    nsample= int(fm[31].split()[1])
    
    tau_particle=np.zeros((Nkpoints,Nbranch)) 
    tau_wave=np.zeros((Nkpoints,Nbranch))    

if calculation_mode=='1':
    
    print(' 1. Only times calculations.')
    
elif calculation_mode=='2':
    
    print(' 1. Times + thermal conductivity calculations.')
    
elif calculation_mode=='3':
    
    print(' 1. Thermal conductivity calculations.')
    filel=fm[32].split()[1]
    filec=fm[33].split()[1]
    
if calculation_mode!='1':
    
    volume=float(fm[34].split()[1])

    if len(fm[35].split())!=2 and len(fm[35].split())!=4:
        
        print('ERROR: The temperature settings are inaccurate...')
        sys.exit()
        
    elif len(fm[35].split())==2:
        
        Temperature=np.zeros((1));Temperature[0]=float(fm[35].split()[1])
        print('    The calculating temperature is %sK.'%(Temperature[0]))
        
    else :
        
        T1=float(fm[35].split()[1]);T2=float(fm[35].split()[2]);T3=float(fm[35].split()[3])
        Temperature=np.arange(T1, T2, T3) 
        print('    The calculating temperature region is [%.1f,%.1f,%.1f]K.'%(T1, T2, T3))
    
    isotope=fm[36].split()[1]
    if isotope=='True':
        
        print('    The isotope scattering is included but from file BTE.w_isotopic.')
        
    size=fm[37].split()[1]
    if size=='True':
        
        if len(fm[37].split())<3:
            print('ERROR: The size settings are inaccurate...')
            sys.exit()
            
        if len(fm[37].split())==3:
            T_size=float(fm[37].split()[2])
            print('   The length dependence will be calculated at T=%.1fK.'%(T_size))
     
if times_mode=='1':
    
    print(' 2. Conventional lifetimes and coherence times will be calculated')
    print('      respcectively from SED and Coherence distribution. ')
    
elif times_mode=='2':
    
    print(' 2. Corrected lifetimes and coherence times will be calculated')
    print('      from phonon decay fitting. ')
    
elif times_mode=='3':
    
    print(' 2. Conventional lifetimes and coherence times will be calculated')
    print('      from spectroscopy or spectral energy. ')

if nsample==1:
    print(' 3. There is %d sample prepared. '%(nsample))
else:
    print(' 3. There are %d samples prepared. '%(nsample))

filename='frequency.dat'
file_have=os.path.exists(filename)

if file_have==False:
    
    print('ERROR: The file of frequency.dat does not exist...')
    sys.exit()
    
else:
    
    f=open(filename,'r')
    fm=f.readlines()
    f.close()
    
    omega=np.zeros((Nkpoints,Nbranch)) 
    
    for i in range(Nkpoints):
        for j in range(Nbranch):
            omega[i,j]=float(fm[i].split()[1+j])

date_info=time.strftime('-%d-%H-%M',time.localtime(time.time()))

#---------------------------------- calculating part ----------------------------------

print(' ')
print('------------------------------------------------------------------------')

if calculation_mode=='3':
    
    print('######## II. Times are reading from files...')
   
    f=open(filel,'r')
    fl=f.readlines()
    f.close()
    
    f=open(filec,'r')
    fc=f.readlines()
    f.close()
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):
            tau_particle[iv,ib]=float(fl[1+ib+iv*Nbranch].split()[1])
            tau_wave[iv,ib]=float(fc[1+ib+iv*Nbranch].split()[1])
    
else: 
    
    print('######## II. Times are calculating...')
    print('------------------------------------------------------------------------')
    print(' ')

#-------------------------------------normal mode analysis------------------------------

    if times_mode=='1':
        
        print(' 1. The fitting of conventional lifetime starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        tau_particle=fitting_models.conventional_taul(date_info, nsample, omega)
                
        print(' 2. The calculation of coherence time starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        tau_wave=fitting_models.density_tauC(date_info, nsample, omega)
         
        print('  Output files i. tau_P%s.dat; ii. tau_C%s.dat'%(date_info,date_info))
    
        print(' ')


#---------------------------------------density-decay------------------------------------
 
    if times_mode=='2':
        
        print(' 1. The calculation of coherence time starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        tau_wave=fitting_models.density_tauC(date_info, nsample, omega) 
        
        print(' 2. The calculation of corrected lifetimes starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        tau_particle=fitting_models.decay_taul(date_info, nsample, omega,tau_wave)
        
        print('  Output files i. tau_P%s.dat; ii. tau_C%s.dat'%(date_info,date_info))
        print(' ')
        
#---------------------------------------spectroscopy_fitting------------------------------------
    
    if times_mode=='3':
        
        print(' 1. The spectroscopy fitting starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        tau_particle,tau_wave=fitting_models.spectroscopy_fitting(date_info, nsample, omega)
       
        print('  Output files i. tau_P%s.dat; ii. tau_C%s.dat'%(date_info,date_info))
        print(' ')
        
    print(' The calculations of different times finish!')

#----------------------------------kappa calculation-----------------------------------------------------------------

if calculation_mode=='2' or calculation_mode=='3':

    print('------------------------------------------------------------------------')
    print('######## III. Thermal conductivities are calculating...')
    print('------------------------------------------------------------------------')
    print('')
        
    filename='Groupvelocity.dat'
    file_have=os.path.exists(filename)

    if file_have==False:
        
        print('ERROR: The file of Groupvelocity.dat does not exist...')
        sys.exit()
        
    else:
        
        f=open(filename,'r')
        inn=f.readlines()
        f.close()                
           
        velocity=np.zeros((Nkpoints*Nbranch,3))
               
        for i in range(Nkpoints*Nbranch):
            for j in range(3):
                velocity[i,j]=float(inn[i].split()[j])
     
    filename='Kpoints.dat'
    file_have=os.path.exists(filename)

    if file_have==False:
        
        print('ERROR: The file of Kpoints.dat does not exist...')
        sys.exit()
        
    else: 
         
        f=open(filename,'r')
        inn=f.readlines()
        f.close()                
           
        weight=np.zeros((Nkpoints))               
        for i in range(Nkpoints):
            weight[i]=float(inn[i+1].split()[3])
 
    #----------------------------------mode kappa-----------------------------------
 
    print(' 1. The calculation of mode thermal conductivity starts:')
    print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    frequency,MFP,mode_kappa_classical_particle,mode_kappa_quantum_particle,mode_kappa_classical_wave,mode_kappa_quantum_wave=\
        kappa_omega.kappa_w(volume, Temperature, omega, tau_particle, tau_wave, velocity, weight, isotope)
    
    print('  Output files i. mode_kappa*dat.')
    print(' ')
    
    #----------------------------------accumulative kappa-----------------------------------
    
    print(' 2. The calculation of frequency accumulative thermal conductivity starts:')
    print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    kappa_accum.kappa_modal(Temperature,frequency,MFP,mode_kappa_classical_particle,mode_kappa_quantum_particle,mode_kappa_classical_wave,mode_kappa_quantum_wave,weight)
    
    print('  Output files i. modal_kappa*dat, ii. accum_kappa*dat.')
    print(' ')
    
    #----------------------------------MFP kappa-----------------------------------
    
    print(' 3. The calculation of MFP accumulative thermal conductivity starts:')
    print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    kappa_accum.kappa_mfp(Temperature,len(frequency),MFP,mode_kappa_classical_particle,mode_kappa_quantum_particle,mode_kappa_classical_wave,mode_kappa_quantum_wave,weight)
        
    print('  Output files i. MFP_kappa*dat, ii. MFP_accum_kappa*dat.')
    print(' ')
    
    #----------------------------------size kappa-----------------------------------

    if size=='True':
                   
        print(' 4. The calculation of length dependent thermal conductivity starts:')
        print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        kappa_L.kappa_length(volume,T_size,tau_particle,tau_wave,velocity,omega,weight,isotope)
        
        print('  Output file i. Length_kappa.dat.')
        print(' ')
        
    print(' The calculation of thermal conductivities finish!')
          
print('------------------------------------------------------------------------')
print('                 Thank you for using WPPT package!')
