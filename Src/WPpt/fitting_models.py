#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 22:52:12 2023

@author: zhongweizhang
"""
import numpy as np
import os
import sys
from scipy.optimize import curve_fit
np.seterr(divide='ignore', invalid='ignore')

def lorentz(w,I0,w0,gamma):

    return I0/((w-w0)**2+gamma)

def decay(t,taul):
    
    global taup

    return np.exp(-t/2.0/taul)*np.exp(-4.0*np.log(2)*t**2/taup**2)

def spectroscopy(x,A,gamma,sigma):

    return A*np.exp(-(x-omega0)**2*gamma)*np.cos(gamma*sigma*(x-omega0))


def spectroscopy_fitting(date_info,nsample,omega):
    
    global omega0    
    
    filename5='tau_P%s.dat'%(date_info)
    fP=open(filename5,'w')
    fP.write('The conventional lifetimes from lorentzian fitting'+'\n') 
    
    filename6='tau_C%s.dat'%(date_info)
    fC=open(filename6,'w')
    fC.write('The coherence times from spectroscopy fitting'+'\n')
           
    filename='NMAw-1/NMAw_1_1.dat'         
    file_have=os.path.exists(filename)
    if file_have==False:
        print('The files of NMAw* do not exist...')
        sys.exit()
        
    else:
        
        f=open(filename,'r')
        am=f.readlines()
        f.close()
        N_total=int(len(am))-3
    
    Nkpoints,Nbranch=np.shape(omega)
    
    tau_particle=np.zeros((Nkpoints,Nbranch)) 
    tau_wave=np.zeros((Nkpoints,Nbranch))    
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='NMAw-1/NMAw_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have==False:
               print('ERROR: The file of %s does not exist...'%(filename))
               sys.exit()
               
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
            
            ffy=np.zeros((N_total))
            fre=np.zeros((N_total))
            
            for i in range(N_total):                    
                ffy[i]=float(am[i+3].split()[1])
                fre[i]=float(am[i+3].split()[0])
                
            for isample in range(1,nsample):
                
                filename='NMAw-%d/NMAw_%d_%d.dat'%(isample+1,iv+1,ib+1)            
                file_have=os.path.exists(filename)
                if file_have==False:
                    print('ERROR: The file of %s does not exist...'%(filename))
                    sys.exit()
                    
                else:
                        
                    f=open(filename,'r')
                    am=f.readlines()
                    f.close() 
                
                for i in range(N_total):                    
                    ffy[i]=ffy[i]+float(am[i+3].split()[1])
                    fre[i]=fre[i]+float(am[i+3].split()[0])
                    
            ffy=ffy/float(nsample)
            fre=fre/float(nsample)
                        
            shift=0.1
            df=fre[1]-fre[0]
                
            p1=np.max([0,int((fre[int(np.argmax(ffy[:int(N_total/2)]))]-shift)/df)])
            p2=int((fre[int(np.argmax(ffy[:int(N_total/2)]))]+shift)/df)
                
            fx=fre[p1:p2]
            fy=ffy[p1:p2]
                
            Aguess=np.max(fy)
            Wguess=fx[int(np.argmax(fy))]
    
            j=0
            for i in range(len(fx)):
                if fy[i]>Aguess*0.2 :
                    j=j+1
    
            Gguess=(j*(fx[2]-fx[1]))**2
                     
            popt, pcov = curve_fit(lorentz,fx,fy,maxfev=500000,\
            bounds=([Aguess*0.3,Wguess*0.5,Gguess*0.3],[Aguess*2,Wguess*1.5,Gguess*1.5]),\
            p0=[np.max(fy),fx[int(np.argmax(fy))],Gguess])
            tau_particle0=1/2/np.sqrt(popt[2])
            tau_particle[iv,ib]=tau_particle0
            
            omega0=popt[1] 
            Gguess=(tau_particle0)**2*np.pi**2/2/np.log(2)
            Sguess=np.pi/(tau_particle0)/4.0/np.pi
            #Bguess=Aguess
            popt, pcov = curve_fit(spectroscopy,fx,fy,bounds=([Aguess*0.3,Gguess*0.7,Sguess*0.7],[Aguess*2,Gguess*2,Sguess*2]),maxfev=500000)
                
            popt[1]=popt[1]/(4.0*np.pi**2/8/np.log(2))
            tau_p=np.sqrt(popt[1])*np.sqrt(2) 
            tau_particle[iv,ib]=1/4/np.pi/popt[2]
            
            tau_wave[iv,ib]=tau_p
            
            fP.write('%f %f '%(omega[iv,ib],tau_particle[iv,ib])+'\n')
            fC.write('%f %f '%(omega[iv,ib],tau_wave[iv,ib])+'\n')
            
    fC.close()
    fP.close()
  
    return tau_particle,tau_wave


def density_tauC(date_info,nsample,omega):
        
    filename3='tau_C%s.dat'%(date_info)
    fC=open(filename3,'w')
    fC.write('The coherence times from density distribution'+'\n')
    
    filename='TauC_Dens-1/coh_1_1.dat'
    file_have=os.path.exists(filename)
    if file_have==False:
        print('The files of TauC_Dens-1/* do not exist...')
        sys.exit()
        
    else:
        
        f=open(filename,'r')
        am=f.readlines()
        f.close()
        N_total=int(len(am))-1
    
    Nkpoints,Nbranch=np.shape(omega)
    tau_wave=np.zeros((Nkpoints,Nbranch))   
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='TauC_Dens-1/coh_%d_%d.dat'%(iv+1,ib+1)
            file_have=os.path.exists(filename)
            if file_have==False:
                print('ERROR: The file of %s does not exist...'%(filename))
                sys.exit()
                
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
            
            Coh_density=np.zeros((N_total,2)) 
            
            for i in range(N_total):                    
                Coh_density[i,0]=float(am[i+1].split()[0])
                Coh_density[i,1]=float(am[i+1].split()[1])
                
            for isample in range(1,nsample):
    
                filename='TauC_Dens-%d/coh_%d_%d.dat'%(isample+1,iv+1,ib+1)
                file_have=os.path.exists(filename)
                
                if file_have==False:
                    print('ERROR: The file of %s does not exist...'%(filename))
                    sys.exit()
                    
                else:
                        
                    f=open(filename,'r')
                    am=f.readlines()
                    f.close()
                  
                    for i in range(N_total):                    
                        Coh_density[i,0]=Coh_density[i,0]+float(am[i+1].split()[0])
                        Coh_density[i,1]=Coh_density[i,1]+float(am[i+1].split()[1])
                         
            Coh_density=Coh_density/float(nsample)
                
            tau_wave[iv,ib]=np.sum(Coh_density[:,0]*Coh_density[:,1])
                
            fC.write('%f  %f '%(omega[iv,ib],tau_wave[iv,ib])+'\n')
    
    fC.close()
    
    return tau_wave


def decay_taul(date_info,nsample,omega,tauC):
     
    global taup
    filename4='tau_P%s.dat'%(date_info)
    fP=open(filename4,'w')
    fP.write('The corrected lifetimes from phonon decay'+'\n')
    
    filename='Decay-1/decay_1_1.dat'         
    file_have=os.path.exists(filename)
    if file_have==False:
        print('The files of Decay-1/* do not exist...')
        sys.exit()
        
    else:
       
        f=open(filename,'r')
        am=f.readlines()
        f.close()
        N_total=int(len(am))-1
    
    Nkpoints,Nbranch=np.shape(omega)
    tau_particle=np.zeros((Nkpoints,Nbranch))   
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='Decay-1/decay_1_1.dat'
            file_have=os.path.exists(filename)
            if file_have==False:
                print('ERROR: The file of %s does not exist...'%(filename))
                sys.exit()
                
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
                N_total=int(len(am))-1
                
                Decay=np.zeros((N_total,2)) 
                
                for i in range(N_total):                    
                    Decay[i,0]=float(am[i+1].split()[0])
                    Decay[i,1]=float(am[i+1].split()[1])
                
            for isample in range(1,nsample):
               
               filename='Decay-%d/decay_%d_%d.dat'%(isample+1,iv+1,ib+1)            
               file_have=os.path.exists(filename)
               if file_have==False:
                   print('ERROR: The file of %s does not exist...'%(filename))
                   sys.exit()
                   
               else:
                       
                   f=open(filename,'r')
                   am=f.readlines()
                   f.close() 
               
               for i in range(N_total):                    
                   Decay[i,0]=Decay[i,0]+float(am[i+1].split()[0])
                   Decay[i,1]=Decay[i,1]+float(am[i+1].split()[1])
                   
            Decay=Decay/float(nsample)
           
            taup=tauC[iv,ib]                
            popt, pcov = curve_fit(decay,Decay[:,0],Decay[:,1],maxfev=500000,p0=[taup])
                
            tau_particle[iv,ib]=abs(popt[0])
                
            fP.write('%f  %f '%(omega[iv,ib],tau_particle[iv,ib])+'\n')
    
    fP.close() 
    
    return tau_particle


def conventional_taul(date_info,nsample,omega):
        
    filename1='tau_P%s.dat'%(date_info)
    fNMA=open(filename1,'w')
    
    fNMA.write('The conventional lifetimes from lorentzian fitting'+'\n')
    
    filename='NMAw-1/NMAw_1_1.dat'
    file_have=os.path.exists(filename)
    if file_have==False:
        print('ERROR: The files of NMAw-1/* do not exist...')
        sys.exit()
        
    else:
        
        f=open(filename,'r')
        am=f.readlines()
        f.close()
        N_total=int(len(am))-3
    
    Nkpoints,Nbranch=np.shape(omega)
    tau_particle=np.zeros((Nkpoints,Nbranch))   
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='NMAw-1/NMAw_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have==False:
                print('ERROR: The file of %s does not exist...'%(filename))
                sys.exit()
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
            
            ffy=np.zeros((N_total))
            fre=np.zeros((N_total))
            
            for i in range(N_total):                    
                ffy[i]=float(am[i+3].split()[1])
                fre[i]=float(am[i+3].split()[0])
                
            for isample in range(1,nsample):
                
                filename='NMAw-%d/NMAw_%d_%d.dat'%(isample+1,iv+1,ib+1)            
                file_have=os.path.exists(filename)
                if file_have==False:
                    print('ERROR: The file of %s does not exist...'%(filename))
                    sys.exit()
                else:
                        
                    f=open(filename,'r')
                    am=f.readlines()
                    f.close() 
                
                for i in range(N_total):                    
                    ffy[i]=ffy[i]+float(am[i+3].split()[1])
                    fre[i]=fre[i]+float(am[i+3].split()[0])
                    
            ffy=ffy/float(nsample)
            fre=fre/float(nsample)
            
            shift=0.1
            df=fre[1]-fre[0]
                
            p1=np.max([0,int((fre[int(np.argmax(ffy[:int(N_total/2)]))]-shift)/df)])
            p2=int((fre[int(np.argmax(ffy[:int(N_total/2)]))]+shift)/df)
                
            fx=fre[p1:p2]
            fy=ffy[p1:p2]
                    
            Aguess=np.max(fy)
            Wguess=fx[int(np.argmax(fy))]
    
            j=0
            for i in range(len(fx)):
                if fy[i]>Aguess*0.2 :
                    j=j+1

            Gguess=(j*(fx[2]-fx[1]))**2
            Bguess=Aguess
            popt, pcov = curve_fit(lorentz,fx,fy,maxfev=500000,\
            bounds=([Aguess*0.3,Wguess*0.5,Gguess*0.3],[Aguess*2,Wguess*1.5,Gguess*1.5]),\
            p0=[np.max(fy),fx[int(np.argmax(fy))],Gguess])
            
            tau_particle[iv,ib]=1/2/np.sqrt(popt[2])
         
            fNMA.write('%f  %f  %f'%(omega[iv,ib],popt[1],tau_particle[iv,ib])+'\n')
    
    fNMA.close()
    
    return tau_particle
