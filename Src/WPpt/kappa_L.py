#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 21:27:49 2023

@author: zhongweizhang
"""
import numpy as np
import os
import sys

global kb
kb=1.38064852e-23

global hbar
hbar=6.62607015e-22         
  
def cv_classical(T,w):
    
    heat=0.0
    heat=heat+kb
                    
    return heat

def cv_quantum(T,w):
    
    heat=0.0
    heat=heat+kb*(hbar*w/kb/T)**2*np.exp(hbar*w/kb/T)\
                    /(np.exp(hbar*w/kb/T)-1)**2
    
    return heat

def kappa_length(volume,T_size,tau_particle,tau_wave,velocity,omega,weight,isotope):
    
    Ntotal=50;Length=np.zeros((Ntotal))
    Length[:10]=np.linspace(0.1,1, 10, endpoint=False)
    Length[10:20]=np.linspace(1,10, 10, endpoint=False)
    Length[20:30]=np.linspace(10,100, 10, endpoint=False)
    Length[30:40]=np.linspace(100,1000, 10, endpoint=False)
    Length[40:50]=np.linspace(1000,10000, 10, endpoint=False)
         
    Kappa_L=np.zeros((Ntotal,4))
    
    Nkpoints,Nbranch=np.shape(tau_particle)
    if isotope=='True':
        
        filename='BTE.w_isotopic'
        file_have=os.path.exists(filename)
    
        if file_have==False:
            print('ERROR: The file of BTE.w_isotopic does not exist...')
            sys.exit()
            
        else:        
            f=open(filename,'r')
            inn1=f.readlines()
            f.close()                
        
            n=len(inn1)
            tau_isotope=np.zeros((n,2))
        
            for i in range(n):
                for j in range(2):
                    tau_isotope[i,j]=float(inn1[i].split()[j])
            
            if n!=Nkpoints*Nbranch:
                print('ERROR: Number of modes in BTE.w_isotopic are inconsistent with nkpoitns*nbranch...')
                sys.exit()
    
    filename='Length_kappa%dK.dat'%(int(T_size))
    f1=open(filename,'w')
    f1.write('L/A  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
    
    for iL in range(Ntotal):
        k=-1
        for j in range(Nbranch):
            for i in range(Nkpoints):
                
                k=k+1;times=tau_particle[i,j]
                ave_vg=np.sqrt((velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0)
                if isotope=='True':                
                    times=1.0/(1.0/times+tau_isotope(k,1)+ave_vg/Length[iL])
                else:
                    times=1.0/(1.0/times+ave_vg/Length[iL])
                    
                Kappa_L[iL,0]=Kappa_L[iL,0]+1.0e-8*cv_classical(T_size,omega[i,j])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*times*weight[i]
                Kappa_L[iL,1]=Kappa_L[iL,1]+1.0e-8*cv_quantum(T_size,omega[i,j])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*times*weight[i]
                
                times=np.sqrt(np.pi/4.0/np.log(2))*tau_wave[i,j]*np.exp(tau_wave[i,j]**2/128.0/np.log(2)/tau_particle[i,j]**2)
                if tau_particle[i,j]==0.0:
                    times=np.sqrt(np.pi/4.0/np.log(2))*tau_wave[i,j]*np.exp(1.0/128.0/np.log(2))
               
                if isotope=='True':                
                    times=1.0/(1.0/times+tau_isotope(k,1)+ave_vg/Length[iL])
                else:
                    times=1.0/(1.0/times+ave_vg/Length[iL]) 
                 
                Kappa_L[iL,2]=Kappa_L[iL,2]+1.0e-8*cv_classical(T_size,omega[i,j])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*weight[i]*times
                Kappa_L[iL,3]=Kappa_L[iL,3]+1.0e-8*cv_quantum(T_size,omega[i,j])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*weight[i]*times
        
        Kappa_L[iL,:]=Kappa_L[iL,:]/np.sum(weight)
        
        f1.write('%f  %f  %f  %f  %f'%(Length[iL],Kappa_L[iL,0],Kappa_L[iL,1],Kappa_L[iL,2],Kappa_L[iL,3])+'\n')
               
    f1.close()        