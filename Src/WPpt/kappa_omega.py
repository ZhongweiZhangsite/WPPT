#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 22:18:15 2023

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

def kappa_w(volume,Temperature,omega,tau_particle,tau_wave,velocity,weight,isotope):
    
    Nkpoints,Nbranch=np.shape(omega)
    frequency=np.zeros((Nkpoints*Nbranch))
    mode_kappa_classical_particle=np.zeros((Nkpoints*Nbranch,len(Temperature)))
    mode_kappa_quantum_particle=np.zeros((Nkpoints*Nbranch,len(Temperature)))
    
    mode_kappa_classical_wave=np.zeros((Nkpoints*Nbranch,len(Temperature)))
    mode_kappa_quantum_wave=np.zeros((Nkpoints*Nbranch,len(Temperature)))
    
    MFP=np.zeros((Nkpoints*Nbranch,2))
    
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
                
    for T in range(len(Temperature)):
        
        filename='mode_kappa%dK.dat'%(int(Temperature[T]))
        f1=open(filename,'w')
        f1.write('omega/THz  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
        
        k=-1    
        for j in range(Nbranch):
            for i in range(Nkpoints):
               
                k=k+1;frequency[k]=omega[i,j]
                times=tau_particle[i,j]
                if isotope=='True':                
                    times=1.0/(1.0/times+tau_isotope(k,1))
                mode_kappa_classical_particle[k,T]=1.0e-8*cv_classical(Temperature[T],frequency[k])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*times*weight[i]
                mode_kappa_quantum_particle[k,T]=1.0e-8*cv_quantum(Temperature[T],frequency[k])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*times*weight[i]
                
                MFP[k,0]=times*np.sqrt((velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0)
                
                times=np.sqrt(np.pi/4.0/np.log(2))*tau_wave[i,j]*np.exp(tau_wave[i,j]**2/128.0/np.log(2)/tau_particle[i,j]**2)
                if tau_particle[i,j]==0.0:
                    times=np.sqrt(np.pi/4.0/np.log(2))*tau_wave[i,j]*np.exp(1.0/128.0/np.log(2))
               
                if isotope=='True':                
                    times=1.0/(1.0/times+tau_isotope(k,1))
                    
                MFP[k,1]=times*np.sqrt((velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0)
                
                mode_kappa_classical_wave[k,T]=1.0e-8*cv_classical(Temperature[T],frequency[k])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*weight[i]*times
                mode_kappa_quantum_wave[k,T]=1.0e-8*cv_quantum(Temperature[T],frequency[k])/volume*(velocity[k,0]**2+velocity[k,1]**2+velocity[k,2]**2)/3.0*weight[i]*times
                
                f1.write('%f  %f  %f  %f  %f'%(frequency[k],mode_kappa_classical_particle[k,T],mode_kappa_quantum_particle[k,T],mode_kappa_classical_wave[k,T],mode_kappa_quantum_wave[k,T])+'\n')
        f1.close()

    return frequency,MFP,mode_kappa_classical_particle,mode_kappa_quantum_particle,mode_kappa_classical_wave,mode_kappa_quantum_wave
