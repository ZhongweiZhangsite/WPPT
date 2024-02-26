#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 21:46:37 2023

@author: zhongweizhang
"""

import numpy as np

def kappa_mfp(Temperature,NN,MFP,mode_k_clas_part,mode_k_quan_part,mode_k_clas_wav,mode_k_quan_wav,weight):
    
    Ntotal=30;MFP_1=np.linspace(0,np.max(MFP)*1.1, Ntotal, endpoint=False) 

    accum_k_clas_part=np.zeros((Ntotal,len(Temperature)))
    accum_k_quan_part=np.zeros((Ntotal,len(Temperature)))

    accum_k_clas_wav=np.zeros((Ntotal,len(Temperature)))
    accum_k_quan_wav=np.zeros((Ntotal,len(Temperature)))

    modal_clas_part=np.zeros((Ntotal,len(Temperature)))
    modal_quan_part=np.zeros((Ntotal,len(Temperature)))

    modal_clas_wav=np.zeros((Ntotal,len(Temperature)))
    modal_quan_wav=np.zeros((Ntotal,len(Temperature)))
    
    for T in range(len(Temperature)):
            
        accum_k_clas_part[0,T]=0.0;accum_k_quan_part[0,T]=0.0
        accum_k_clas_wav[0,T]=0.0;accum_k_quan_wav[0,T]=0.0
           
        modal_clas_part[0,T]=0.0;modal_quan_part[0,T]=0.0
        modal_clas_wav[0,T]=0.0;modal_quan_wav[0,T]=0.0
             
        for j in range(1,Ntotal):
    
            for i in range(NN):
                   
                if MFP_1[j-1]<MFP[i,0] and MFP_1[j]>MFP[i,0] :
    
                    accum_k_clas_part[j,T]=accum_k_clas_part[j,T]+mode_k_clas_part[i,T]    
                    accum_k_quan_part[j,T]=accum_k_quan_part[j,T]+mode_k_quan_part[i,T]
                        
                    modal_clas_part[j,T]=modal_clas_part[j,T]+mode_k_clas_part[i,T]    
                    modal_quan_part[j,T]=modal_quan_part[j,T]+mode_k_quan_part[i,T]
                       
                         
                if MFP_1[j-1]<MFP[i,1] and MFP_1[j]>MFP[i,1]:
                         
                    accum_k_clas_wav[j,T]=accum_k_clas_wav[j,T]+mode_k_clas_wav[i,T]    
                    accum_k_quan_wav[j,T]=accum_k_quan_wav[j,T]+mode_k_quan_wav[i,T]
                     
                    modal_clas_wav[j,T]=modal_clas_wav[j,T]+mode_k_clas_wav[i,T]    
                    modal_quan_wav[j,T]=modal_quan_wav[j,T]+mode_k_quan_wav[i,T]
               
                if 1==1:
                    accum_k_quan_part[j,T]=accum_k_quan_part[j,T]+accum_k_quan_part[j-1,T]
                    accum_k_clas_part[j,T]=accum_k_clas_part[j,T]+accum_k_clas_part[j-1,T]
                   
                    accum_k_quan_wav[j,T]=accum_k_quan_wav[j,T]+accum_k_quan_wav[j-1,T]
                    accum_k_clas_wav[j,T]=accum_k_clas_wav[j,T]+accum_k_clas_wav[j-1,T]
        
        accum_k_clas_part[:,T]=accum_k_clas_part[:,T]/np.sum(weight)
        accum_k_quan_part[:,T]=accum_k_quan_part[:,T]/np.sum(weight)
       
        accum_k_clas_wav[:,T]=accum_k_clas_wav[:,T]/np.sum(weight)
        accum_k_quan_wav[:,T]=accum_k_quan_wav[:,T]/np.sum(weight)
    
        modal_clas_part[:,T]=modal_clas_part[:,T]/np.sum(weight)
        modal_quan_part[:,T]=modal_quan_part[:,T]/np.sum(weight)
       
        modal_clas_wav[:,T]=modal_clas_wav[:,T]/np.sum(weight)
        modal_quan_wav[:,T]=modal_quan_wav[:,T]/np.sum(weight)
         
        filename='MFP_kappa%dK.dat'%(int(Temperature[T]))
        f1=open(filename,'w')
        f1.write('omega/THz  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
            
        filename='MFP_accum_kappa%dK.dat'%(int(Temperature[T]))
        f2=open(filename,'w')
        f2.write('MFP/A  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
             
        f1.write('%f  %f  %f  %f  %f'%(MFP_1[0],0.0,0.0,0.0,0.0)+'\n')
        f2.write('%f  %f  %f  %f  %f'%(MFP_1[0],0.0,0.0,0.0,0.0)+'\n')
        for j in range(1,Ntotal):
            f1.write('%f  %f  %f  %f  %f'%(MFP_1[j],modal_clas_part[j,T],modal_quan_part[j,T],modal_clas_wav[j,T],modal_quan_wav[j,T])+'\n')
            f2.write('%f  %f  %f  %f  %f'%(MFP_1[j],accum_k_clas_part[j,T],accum_k_quan_part[j,T],accum_k_clas_wav[j,T],accum_k_quan_wav[j,T])+'\n')
     
        f1.close() 
        f2.close()
 
    
def kappa_modal(Temperature,frequency,MFP,mode_k_clas_part,mode_k_quan_part,mode_k_clas_wav,mode_k_quan_wav,weight):
    
    Ntotal=30;omega_1=np.linspace(0,np.max(frequency)*1.1, Ntotal, endpoint=False)
    
    accumulative_kappa_classical_particle=np.zeros((Ntotal,len(Temperature)))
    accumulative_kappa_quantum_particle=np.zeros((Ntotal,len(Temperature)))
    
    accumulative_kappa_classical_wave=np.zeros((Ntotal,len(Temperature)))
    accumulative_kappa_quantum_wave=np.zeros((Ntotal,len(Temperature)))
    
    modal_classical_particle=np.zeros((Ntotal,len(Temperature)))
    modal_quantum_particle=np.zeros((Ntotal,len(Temperature)))
    
    modal_classical_wave=np.zeros((Ntotal,len(Temperature)))
    modal_quantum_wave=np.zeros((Ntotal,len(Temperature)))
    
    for T in range(len(Temperature)):
        
        accumulative_kappa_classical_particle[0,T]=0.0;accumulative_kappa_quantum_particle[0,T]=0.0
        accumulative_kappa_classical_wave[0,T]=0.0;accumulative_kappa_quantum_wave[0,T]=0.0
       
        modal_classical_particle[0,T]=0.0;modal_quantum_particle[0,T]=0.0
        modal_classical_wave[0,T]=0.0;modal_quantum_wave[0,T]=0.0 
       
        for j in range(1,Ntotal):
             k1=0
             for i in range(len(frequency)):
               
                 if omega_1[j-1]<frequency[i] and omega_1[j]>frequency[i] :
                     k1=k1+1
                     accumulative_kappa_classical_particle[j,T]=accumulative_kappa_classical_particle[j,T]+mode_k_clas_part[i,T]    
                     accumulative_kappa_quantum_particle[j,T]=accumulative_kappa_quantum_particle[j,T]+mode_k_quan_part[i,T]
                    
                     accumulative_kappa_classical_wave[j,T]=accumulative_kappa_classical_wave[j,T]+mode_k_clas_wav[i,T]    
                     accumulative_kappa_quantum_wave[j,T]=accumulative_kappa_quantum_wave[j,T]+mode_k_quan_wav[i,T]
                   
                     modal_classical_particle[j,T]=modal_classical_particle[j,T]+mode_k_clas_part[i,T]    
                     modal_quantum_particle[j,T]=modal_quantum_particle[j,T]+mode_k_quan_part[i,T]
                   
                     modal_classical_wave[j,T]=modal_classical_wave[j,T]+mode_k_clas_wav[i,T]    
                     modal_quantum_wave[j,T]=modal_quantum_wave[j,T]+mode_k_quan_wav[i,T]
           
             if 1==1:
                accumulative_kappa_quantum_particle[j,T]=accumulative_kappa_quantum_particle[j,T]+accumulative_kappa_quantum_particle[j-1,T]
                accumulative_kappa_classical_particle[j,T]=accumulative_kappa_classical_particle[j,T]+accumulative_kappa_classical_particle[j-1,T]
               
                accumulative_kappa_quantum_wave[j,T]=accumulative_kappa_quantum_wave[j,T]+accumulative_kappa_quantum_wave[j-1,T]
                accumulative_kappa_classical_wave[j,T]=accumulative_kappa_classical_wave[j,T]+accumulative_kappa_classical_wave[j-1,T]
         
    
        accumulative_kappa_classical_particle[:,T]=accumulative_kappa_classical_particle[:,T]/np.sum(weight)
        accumulative_kappa_quantum_particle[:,T]=accumulative_kappa_quantum_particle[:,T]/np.sum(weight)
       
        accumulative_kappa_classical_wave[:,T]=accumulative_kappa_classical_wave[:,T]/np.sum(weight)
        accumulative_kappa_quantum_wave[:,T]=accumulative_kappa_quantum_wave[:,T]/np.sum(weight)
    
        modal_classical_particle[:,T]=modal_classical_particle[:,T]/np.sum(weight)
        modal_quantum_particle[:,T]=modal_quantum_particle[:,T]/np.sum(weight)
       
        modal_classical_wave[:,T]=modal_classical_wave[:,T]/np.sum(weight)
        modal_quantum_wave[:,T]=modal_quantum_wave[:,T]/np.sum(weight)
        
        filename='modal_kappa%dK.dat'%(int(Temperature[T]))
        f1=open(filename,'w')
        f1.write('omega/THz  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
     
        filename='accum_kappa%dK.dat'%(int(Temperature[T]))
        f2=open(filename,'w')
        f2.write('omega/THz  p-k_clas  p-k_quan  w-k_clas  w-k_quan'+'\n')
      
        f1.write('%f  %f  %f  %f  %f'%(omega_1[0],0.0,0.0,0.0,0.0)+'\n')
        f2.write('%f  %f  %f  %f  %f'%(omega_1[0],0.0,0.0,0.0,0.0)+'\n')
        for j in range(1,Ntotal):
            f1.write('%f  %f  %f  %f  %f'%(omega_1[j],modal_classical_particle[j,T],modal_quantum_particle[j,T],modal_classical_wave[j,T],modal_quantum_wave[j,T])+'\n')
            f2.write('%f  %f  %f  %f  %f'%(omega_1[j],accumulative_kappa_classical_particle[j,T],accumulative_kappa_quantum_particle[j,T],accumulative_kappa_classical_wave[j,T],accumulative_kappa_quantum_wave[j,T])+'\n')
        
        f1.close() 
        f2.close()    
     