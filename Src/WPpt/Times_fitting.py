import numpy as np
from scipy.optimize import curve_fit

import os
import sys

def lorentz(w,I0,w0,gamma):

    return I0/((w-w0)**2+gamma)

def decay(t,taul):
    
    global taup

    return np.exp(-t/2.0/taul)*np.exp(-4.0*np.log(2)*t**2/taup**2)

def proposed0(x,A,gamma,sigma):

    return A*np.exp(-(x-omega0)**2*gamma)*np.cos(gamma*sigma*(x-omega0))

print('-----------This python code used for modal fitting------------')
print('---------------This is a part of WPPT package-----------------')
print(' ')
print('----------------The homepage for WPPT package-----------------')
print('----------https://github.com/ZhongweiZhangsite/WPPT-----------')
print(' ')
print('The fitted properties can be further applied to WPpt module.')
print('The command is reading from CONTROL.WPPT ..')
print(' ')

filename='CONTROL.WPPT'

file_have=os.path.exists(filename)
if file_have=='False':
    print('The target file CONTROL.WPPT does not exist...')
else:
                    
    f=open(filename,'r')
    fm=f.readlines()
    f.close()
    
    Nbranch=int(fm[6].split()[1])*3
    Nkpoints=int(fm[8].split()[1])
    
    tau_NMA= fm[31].split()[1]
    tau_c= fm[32].split()[1]
    tau_p= fm[33].split()[1]
    
    tau_SSC= fm[37].split()[1]
    
print('1. The coventional lifetime from NMA: %s'%(tau_NMA))
print('2. The coherence time from PWave:     %s'%(tau_c))
print('3. The corrected lifetime from PWave: %s'%(tau_p))
print('4. The spectroscopy fitting: %s'%(tau_SSC))
print(' ')

Nbranch=int(fm[6].split()[1])*3
Nkpoints=int(fm[8].split()[1])

filename='frequency.dat'
file_have=os.path.exists(filename)

if file_have==False:
    print('The target file frequency.dat does not exist...')
    sys.exit()
else:
        
    f=open(filename,'r')
    fm=f.readlines()
    f.close()
    
    omega=np.zeros((Nkpoints,Nbranch))
    
    for i in range(Nkpoints):
        for j in range(Nbranch):
            omega[i,j]=float(fm[i].split()[1+j])


if tau_NMA=='True' :
    
    print('1. The fitting of coventional lifetime start:')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    filename1='tau_NMA.dat'
    fNMA=open(filename1,'w')
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='NMAw_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have=='False':
                print('The target files NMAw does not exist...')
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
                N_total=int(len(am))-3
                
                ffy=np.zeros((N_total))
                fre=np.zeros((N_total))
                
                for i in range(N_total):                    
                    ffy[i]=float(am[i+3].split()[1])
                    fre[i]=float(am[i+3].split()[0])
                            
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
                     
                popt, pcov = curve_fit(decay,fx,fy,maxfev=500000,bounds=([Aguess*0.3,Wguess*0.8,Gguess*0.3,Bguess*0.3],[Aguess*2,Wguess*1.2,Gguess*1.5,Bguess*2]),p0=[np.max(fy),fx[int(np.argmax(fy))],50.0])
                tau_pNMA=1/2/np.sqrt(popt[2])
                
                fNMA.write('%f  %f  %f'%(omega[iv,ib],popt[1],tau_pNMA)+'\n')
    
    fNMA.close()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
print(' ')

if tau_c=='True' :  
    
    print('2. The calculation of coherence time start:')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    filename2='tau_C.dat'
    fC=open(filename2,'w')
    
    tau_coh=np.zeros((Nkpoints,Nbranch)) 
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='coh_density_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have=='False':
                print('The target file does not exist...')
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
                N_total=int(len(am))-1
                
                Coh_density=np.zeros((N_total,2)) 
                
                for i in range(N_total):                    
                    Coh_density[i,0]=float(am[i+1].split()[0])
                    Coh_density[i,1]=float(am[i+1].split()[1])
                            
                tau_coh[iv,ib]=np.sum(Coh_density[:,0]*Coh_density[:,1])
                
                fC.write('%f  %f '%(omega[iv,ib],tau_coh[iv,ib])+'\n')
    
    fC.close()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
print(' ')

if tau_p=='True' :
    
    print('3. The fitting of lifetime time start:')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    global taup
    filename3='tau_P.dat'
    fP=open(filename3,'w')
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):            
            
            filename='Decay_coh_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have=='False':
                print('The target file does not exist...')
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
                N_total=int(len(am))-3
                
                Decay=np.zeros((N_total,2)) 
                
                for i in range(N_total):                    
                    Decay[i,0]=float(am[i+3].split()[0])
                    Decay[i,1]=float(am[i+3].split()[1])
                            
                taup=tau_coh[iv,ib]
                
                popt, pcov = curve_fit(decay,Decay[:,0],Decay[:,1],maxfev=500000,p0=[taup])
                
                tau_P=abs(popt[0])
                
                fP.write('%f  %f '%(omega[iv,ib],tau_P)+'\n')
    
    fP.close()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')

    print(' ')

if tau_SSC=='True' :
    
    print('4. The spectroscopy fitting start:')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    global omega0
    
    filename4='SSC_tau.dat'
    f_SCC=open(filename4,'w')
    
    for iv in range(Nkpoints):
        for ib in range(Nbranch):  
            
            filename='NMAw_%d_%d.dat'%(iv+1,ib+1)            
            file_have=os.path.exists(filename)
            if file_have=='False':
                print('The target files NMAw does not exist...')
            else:
                    
                f=open(filename,'r')
                am=f.readlines()
                f.close()
                N_total=int(len(am))-3
                
                ffy=np.zeros((N_total))
                fre=np.zeros((N_total))
                
                for i in range(N_total):                    
                    ffy[i]=float(am[i+3].split()[1])
                    fre[i]=float(am[i+3].split()[0])
                            
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
                     
                popt, pcov = curve_fit(lorentz,fx,fy,maxfev=500000,bounds=([Aguess*0.3,Wguess*0.8,Gguess*0.3,Bguess*0.3],[Aguess*2,Wguess*1.2,Gguess*1.5,Bguess*2]),p0=[np.max(fy),fx[int(np.argmax(fy))],50.0])
                tau_l0=1/2/np.sqrt(popt[2])
                
                omega0=popt[1] 
                Gguess=(tau_l0)**2*np.pi**2/2/np.log(2)
                Sguess=np.pi/(tau_l0)/4.0/np.pi
                #Bguess=Aguess
                popt, pcov = curve_fit(proposed0,fx,fy,bounds=([Aguess*0.3,Gguess*0.7,Sguess*0.7],[Aguess*2,Gguess*2,Sguess*2]),maxfev=500000)
                
                popt[1]=popt[1]/(4.0*np.pi**2/8/np.log(2))
                tau_p=np.sqrt(popt[1])*np.sqrt(2) 
                tau_l=1/4/np.pi/popt[2]
                
                f_SCC.write('%f %f %f %f %f'%(omega0,tau_l0,tau_l,tau_p,tau_p/tau_l)+'\n')
    
    f_SCC.close()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%')


print('The calculation and fitting for different times properties finish!')
print(' ')
print('Thank you for your using!')