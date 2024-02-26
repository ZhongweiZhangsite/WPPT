 !GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.
 
module temporal

    contains
    subroutine temporal_coh

use wl_trans_temp
use cor_transform
use cudafor

implicit none
integer:: i,j,k,Ln,num_time,width_number,ii,ierr
integer:: Min_time,Max_time,Step_time,displace,Step_displace
integer:: n_f,n_k,num_cell,na,nb,nc,n,ik,ib,wvector,branch,f(1)
integer:: i_sample,n_sample,n_interval,n1,n2,cut_off
double precision:: pos(2),width,dt
double precision:: time0,time1
double precision:: omega_target,tau_up,tau_down,sc_factor1,Size_WP
double precision:: Min_width,Max_width,Step_width
double precision,allocatable:: wavelet_results_node(:,:),coh_density(:),width_coh(:)
double precision,allocatable:: omega0(:),time(:),omega(:),fft_node(:,:)
double precision,allocatable:: cor_each(:,:),cor_ave(:),wave_in(:)
double complex,allocatable::modal_velocity(:)

double precision,device:: omega_targetd,dtd,Size_WPd
integer,device:: Step_displaced
double precision,device,allocatable:: width_cohd(:),timed(:),wave_ind(:),corrd(:)
double precision,device,allocatable:: wavelet_results_noded(:,:),omegad(:)
double complex,device,allocatable:: modal_velocityd(:)

real, parameter :: pi=3.1415926,mass=1.66054e-27,A2m=1.0e-10,ps2s=1.0e-12,ev2J=1.60218e-19,hbar=6.626e-34

character:: char0,TCD_WC,Omega_OR,decay_WC,each_WC
character(len=10) :: file_id
character(len=50) :: file_name

CHARACTER(256) :: folder_path
LOGICAL :: folder_exists

integer :: plan

integer::tPB=32
type(dim3) :: grid, tBlock

!-----------------------parameters setting------------------------

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,Ln
read(1,*)
read(1,*) char0,num_cell
read(1,*) char0,n_k
do i=1,3
read(1,*)
enddo
n_f=num_cell*3
read(1,*) char0,dt
do i=1,7
    read(1,*)
enddo
read(1,*) char0,Step_displace
read(1,*) char0,tau_up
read(1,*) char0,tau_down
read(1,*) char0,width_number
read(1,*) char0,Size_WP
read(1,*) char0,n_sample
read(1,*) char0,Omega_OR
read(1,*) char0,TCD_WC
read(1,*) char0,decay_WC
read(1,*) char0,each_WC
close(1)

dtd=dt
Size_WPd=Size_WP
Step_displaced=Step_displace

!-----------------------------------Initialization-----------------------------------------

allocate(time(Ln))
allocate(timed(Ln))

folder_path = './NMAt'
INQUIRE(FILE=folder_path, EXIST=folder_exists)

IF(folder_exists) THEN
    write(file_id, '(i0,A1,i0)') 1,'_',1
    file_name = 'NMAt/NMAt_' // trim(adjustl(file_id)) // '.dat'
     
    Open(1,file=trim(file_name))
    read(1,*)
    read(1,*)
    read(1,*)
    do i=1,Ln
        read(1,*) time(i)
    enddo
    close(1)
ELSE
    write(*,*) "ERROR: The files in NMAt/* do not exist..."
    stop
ENDIF
 
num_time=int(Ln/Step_displace)

!----------------------parallel message sending------------------

write(*,*) "######## I. Successfully read pre-prepared files!"
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
write(*,'(A24,I9,A2,f15.3,A5)') 'The nstep, totaltime:',Ln,',',maxval(time),'ps.'
write(*,'(A25,I5,A21,I5,A2)') 'Number of wavevectors:',n_k,'; Number of branches:',n_f,'.'
write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'

write(*,*) '######## II. Setup for the phonon wavelet transform...'
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '

!-------------------------------------frequency/vector----------------------------------

Min_width=tau_down !/dt/(2.0*sqrt(2*log(2.0)))      !1.0/maxval(vector_target(:2))*sc_factor1
Max_width=tau_up !/dt/(2.0*sqrt(2*log(2.0)))        !1.0/maxval(vector_target(:2))*tau_up  ! considering the maximum coherence time
Step_width=(Max_width-Min_width)/float(width_number)

cut_off=int(Max_width*Size_WP/dt/Step_displace)         ! cut-off at the begin and end

allocate(width_coh(width_number))
allocate(width_cohd(width_number))

do i=1,width_number
width_coh(i)=Min_width+float(i-1)*Step_width
enddo

width_cohd=width_coh

if(Max_width>maxval(time)*0.9) then
write(*,*) 'ERROR: Maximun coherence time is not reasonable...'
stop
endif

write(*,'(A30,f15.3,A5)') ' The minimun coherence time:',tau_down, 'ps.'
write(*,'(A30,f15.3,A5)') ' The maximun coherence time:',tau_up, 'ps.'
write(*,'(A28,I9,A2)') 'The coherence time steps:',width_number,'.'
write(*,'(A39,I9,A2)') ' The time evolution steps in Wavelet:',num_time,'.'
write(*,'(A23,I9,A2)') ' The cut-off steps:  ',cut_off,'.'

write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'
write(*,*) "######## III. GPU parallel calculation is beginning..."

if(Max_width*2.0*Size_WP>maxval(time)) then
write(*,*) 'ERROR: The coherence time is too long...'
stop
endif

write(*,*) ' '

allocate(wavelet_results_node(width_number,num_time))
allocate(wavelet_results_noded(width_number,num_time))
wavelet_results_node=0.0d0

allocate(modal_velocity(Ln))
allocate(modal_velocityd(Ln))
allocate(coh_density(width_number))

folder_path = './TauC_Dens'
INQUIRE(FILE=folder_path, EXIST=folder_exists)

IF(folder_exists) THEN
    call system("rm -fr TauC_Dens/*") 
ELSE
    call system("mkdir TauC_Dens")
ENDIF

if(TCD_WC=='T') then
    folder_path = './T_C_D'
    INQUIRE(FILE=folder_path, EXIST=folder_exists)

    IF(folder_exists) THEN
        call system("rm -fr T_C_D/*")
    ELSE
        call system("mkdir T_C_D")
    ENDIF
else
    folder_path = './T_C_D'
    INQUIRE(FILE=folder_path, EXIST=folder_exists)

    IF(folder_exists) THEN
        call system("rm -fr T_C_D")
    ENDIF
  
endif

if(decay_WC=='T') then

    folder_path = './Decay'
    INQUIRE(FILE=folder_path, EXIST=folder_exists)

    IF(folder_exists) THEN
        call system("rm -fr Decay/*")
    ELSE
        call system("mkdir Decay")
    ENDIF

    n_interval=int((num_time-2.0*cut_off)/5.0)
    allocate(cor_each(width_number,n_interval))
    allocate(cor_ave(n_interval))
    allocate(wave_in(n_interval))
    allocate(wave_ind(n_interval))
    allocate(corrd(n_interval))

    if(each_WC=='T') then

        folder_path = './Decay_coh'
        INQUIRE(FILE=folder_path, EXIST=folder_exists)
        IF(folder_exists) THEN
            call system("rm -fr Decay_coh/*")
        ELSE
            call system("mkdir Decay_coh")
        ENDIF

    endif

else

    if(each_WC=='T') then
        write(*,'(A40)') ' The settings for decay are inconsistent.'
    endif

endif


if(Omega_OR=='F') then
    allocate(omega(Ln))
    allocate(omegad(Ln))
    do i=1,Ln
        omega(i)=1.0/dt/float(Ln)*float(i)
    enddo
    omegad=omega

    allocate(fft_node(Ln,2))
  !  allocate(fft_noded(Ln))
         
endif

write(*,'(A43)') '                      wavevector-branch:'

!-------------------------------------calculating----------------------------------

do ik=1,n_k
    do ib=1,n_f

        write(*,'(A25,2I7)') '   In calculating mode',ik,ib
        
        write(file_id, '(i0,A1,i0)') ik,'_',ib
        file_name = 'NMAt/NMAt_' // trim(adjustl(file_id)) // '.dat'
         
        Open(1,file=trim(file_name))
        read(1,*)
        read(1,*) char0,wvector,char0,branch
        read(1,*)

        do i=1,Ln
            read(1,*) time(i),pos(:)
            modal_velocity(i)=cmplx(pos(1),pos(2))
        enddo
        close(1)
        
        timed=time
       
        modal_velocityd=modal_velocity

        if(Omega_OR=='T') then

            open(1,file='frequency.dat')
            do i=1,(wvector-1)
            read(1,*)
            enddo
            allocate(omega0(branch+1))
            read(1,*) omega0
            omega_target=omega0(branch+1)
            close(1)
            deallocate(omega0)
 
        else
             
            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'NMAw/NMAw_' // trim(adjustl(file_id)) // '.dat'
            Open(1,file=trim(file_name))
            read(1,*)
            read(1,*)
            read(1,*)
            ii=0
            do while(.true.)
                read(1,*,iostat=ierr)  pos(:)
                if(ierr/=0)exit
                    ii=ii+1
                    fft_node(ii,:)=pos(:)
            enddo
            close(1)

            f(1)=maxloc(fft_node(:,2),1)
            omega_target=fft_node(int(f(1)),1)
 
        endif

!-----------------------------wavelet transform-------------------------------------------

        omega_targetd=omega_target

        tBlock = dim3(tPB,1,1)
        grid = dim3(ceiling(real(num_time)/tBlock%x),1,1)

        call Cal_temp_1D<<<grid, tBlock>>>(modal_velocityd,width_cohd,omega_targetd,&
                Size_WPd,timed,dtd,Step_displaced,wavelet_results_noded)

        wavelet_results_node=wavelet_results_noded
        wavelet_results_node=wavelet_results_node*0.5*mass*(A2m/ps2s)**2/(hbar/ps2s*omega_target)
        wavelet_results_node=wavelet_results_node*2.0

!-----------------------------output transformc-------------------------------------------

        if(TCD_WC=='T') then

            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'T_C_D/time_tau_' // trim(adjustl(file_id)) // '.dat'
            
            open(1,file=trim(file_name))
            write(1,*) 't_evo   tau_c   Ph_number'

            do i=cut_off,num_time-cut_off
                do j=1,width_number
                    write(1,'(3E15.7)') time(i*Step_displace),width_coh(j),wavelet_results_node(j,i)
                enddo
            enddo
            close(1)

        endif

        coh_density=0.0d0

        do i=cut_off,num_time-cut_off
            coh_density=coh_density+wavelet_results_node(:,i)
        enddo
        
        coh_density=coh_density/sum(coh_density)

        write(file_id, '(i0,A1,i0)') ik,'_',ib
        file_name = 'TauC_Dens/coh_' // trim(adjustl(file_id)) // '.dat'

        open(1,file=trim(file_name))
        write(1,*) 'tau_c   SAPND'
        do j=1,width_number
            write(1,*)  width_coh(j),coh_density(j)
        enddo
        close(1)

!-----------------------------------------phonon decay------------------------------

        if(decay_WC=='T') then
             
            cor_each=0.0

            tBlock = dim3(tPB,1,1)
            grid = dim3(ceiling(real(n_interval)/tBlock%x),1,1)

            do j=1,width_number

                    do i_sample=1,n_sample

                        n1=(i_sample-1)*int((num_time-2.0*cut_off-n_interval)/(n_sample-1))+1+cut_off
                        n2=(i_sample-1)*int((num_time-2.0*cut_off-n_interval)/(n_sample-1))+n_interval+cut_off
        
                        wave_in(:)=wavelet_results_node(j,n1:n2);wave_in(:)=wave_in(:)-sum(wave_in(:))/float(n_interval)
                        wave_ind=wave_in
  
                        call Cal_cor1D<<<grid, tBlock>>>(wave_ind,corrd)
                        cor_each(j,:)=cor_each(j,:)+corrd
                    
                    enddo
             
                    cor_each(j,:)=cor_each(j,:)/cor_each(j,1)

            enddo
              
            if(each_WC=='T') then

                write(file_id, '(i0,A1,i0)') ik,'_',ib
                file_name = 'Decay_coh/decay_coh_' // trim(adjustl(file_id)) // '.dat'
                open(1,file=trim(file_name))
                write(1,*) 'evo_t   decay',width_number,n_interval

                do j=1,width_number
                    write(1,*) width_coh(j)
                    do i=1,n_interval
                        write(1,*)  time(i*Step_displace),cor_each(j,i)
                    enddo
                enddo
                close(1)

            endif

            cor_ave=0.0
            do j=1,width_number
                cor_ave=cor_ave+coh_density(j)*cor_each(j,:)
            enddo

            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'Decay/decay_' // trim(adjustl(file_id)) // '.dat'
            open(1,file=trim(file_name))
            write(1,*) 'evo_t decay'

            do i=1,n_interval
                write(1,*)  time(i*Step_displace),cor_ave(i)
            enddo
            close(1)

        endif

    enddo
enddo

!------------------------------------------------------------------------
 write(*,*) ' '
 write(*,*) "  Coherence calculation is successfully done!"
 write(*,*) '------------------------------------------------------------------------'
 write(*,*) "######## IV. The output files:"
 write(*,'(A23)',advance='no') ' 1. TauC_Dens/coh_*.dat'
 n=1
 if(TCD_WC=='T') then
     n=n+1
     write(*,'(A2,I2,A22)',advance='no') ' ',n,'. T_C_D/time_tau_*.dat'
 endif

 if(decay_WC=='T') then
    n=n+1
    write(*,'(A2,I2,A19)',advance='no') ' ',n,'. Decay/decay_*.dat'
 endif

 if(each_WC=='T') then
    n=n+1
    write(*,'(A2,I2,A27)',advance='no') ' ',n,'. Decay_coh/decay_coh_*.dat'
 endif

 write(*,*) '.'

end subroutine temporal_coh
end module temporal
