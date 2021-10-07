!GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the phonon coherence from a parallel running code

program main

use mpi
include 'fftw3.f'

integer i,id,ii,jj,m,xx,num,Ln,num_time,width_number,f(1)
integer Min_time,Max_time,Step_time,displace,Step_displace
integer wvector,branch,sc_factor2
double precision,allocatable:: time_number_node(:),time_id(:)
double precision pos(2),omega_target,wavelet_results1,width,dt
double precision tau_up,sc_factor1,Size_WP
double precision Min_width,Max_width,Step_width,Min_displace,Max_displace,real_time,process,wavelet_results2
double complex::  temp1,temp2
double precision,allocatable::time(:),width_coh(:),coherence_dist(:,:),coherence_dist0(:,:),wavelet_results_node(:,:),coh_density(:)
double precision,allocatable:: omega_range(:),coh_spectral(:,:),w_t_coh(:,:,:),omega0(:)
double complex,allocatable::modal_velocity(:),q0(:),wavelet_results(:),motherwavelet(:)
integer :: plan
real, parameter :: pi=3.1415926,mass=1.993e-26,Atom=1.0e-10,pstos=1.0e-12,hbar=6.62607015e-22
character tmp

integer(kind=MPI_OFFSET_KIND)    offset
integer  fh,nbins,nlest,recvcount
integer,allocatable:: sendcount(:),displs(:)
character*16 infile1,infile2,infile3
character char0
character Omega_cw,Spectral_cw
character(len=10) :: file_id
character(len=50) :: file_name
integer :: tx
character( Len = 85 ) :: cStr
integer ierr,my_id,num_procs,root_proc,filetype

!----------------------------------------------------------------

call CPU_TIME(time0)

root_proc=0
    
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

if(my_id==root_proc)then

    write(*,*) '------------------------------------------------------------'
    write(*,*) '----***----******---***---*******---*******--********-------'
    write(*,*) '-----***---******---***--***---**--***---**-----**----------'
    write(*,*) '-----***--***-***--***---********--********-----**----------'
    write(*,*) '------***-***--***-***---***-------***----------**----------'
    write(*,*) '------*******--******----***-------***----------**----------'
    write(*,*) '------*******--******----**--------**-----------**---V.1.0--'
    write(*,*) '------------------------------------------------------------'
    write(*,*) '-------------The PWave module in WPPT package---------------'
    write(*,*) '---------------------By Zhongwei Zhang----------------------'
    write(*,*) '---------------Email: zhongwei@tongji.edu.cn----------------'
    write(*,*) '---------------------Date: Oct. 7 2021----------------------'
    write(*,*) '------------------------------------------------------------'

endif

CALL getarg(1,infile1)

!-----------------------parameters setting------------------------

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,Ln
read(1,*) char0,dt
read(1,*)
do i=1,15
read(1,*)
enddo
read(1,*) char0,Step_displace
read(1,*) char0,Omega_cw
read(1,*) char0,tau_up
read(1,*) char0,sc_factor1
read(1,*) char0,sc_factor2
read(1,*) char0,Size_WP
read(1,*) char0,Spectral_cw
close(1)

!-----------------------------------Initialization-----------------------------------------

allocate(modal_velocity(Ln))
allocate(time(Ln))

Open(1,file=infile1)
read(1,*)
read(1,*) char0,wvector,char0,branch
read(1,*)
do i=1,Ln
	read(1,*) time(i),pos(:)
	modal_velocity(i)=cmplx(pos(1),pos(2))
enddo
close(1)
 
!----------------------parallel message sending------------------

num_time=int(Ln/Step_displace)

allocate(sendcount(num_procs))
allocate(displs(num_procs))

nbins=num_time/num_procs
nlest=mod(num_time,num_procs)

do i=1,num_procs

    if(i>nlest) then
        sendcount(i)=nbins
    else
        sendcount(i)=nbins+1
    endif

enddo
 
displs(1)=0
do i=2,num_procs
     displs(i)=displs(i-1)+sendcount(i-1)
     displs(i)=displs(i)
enddo

do i=1,num_procs
    
    if(my_id==(i-1)) then
        allocate(time_number_node(sendcount(i)))
        recvcount=sendcount(i)
    endif
enddo

!-----------------------------------------------------------------

if (my_id==root_proc) then
    
    allocate(time_id(num_time))

    do i=1,num_time
        time_id(i)=float(i*Step_displace)
    enddo

    write(*,*) "######## 1. Successfully read pre-prepared files!"
    write(*,*) "------------------------------------------------------------"
    write(*,*) ' '
    write(*,'(A21,I9,A3,f15.3,A4)') ' The nstep, totaltime:',Ln,',',time(Ln),'ps.'
    write(*,'(A14,f15.3,A4)') ' The timestep:',dt,'ps.'
    write(*,*) ' '
    write(*,*) "-----------------------------------------------------------"
    
endif

if (my_id==root_proc) then
    write(*,*) '######## 2. Determining the frequency...'
    write(*,*) '------------------------------------------------------------'
    write(*,*) ' '
    if(Omega_cw=='T') then
        write(*,*) 'The frequency is read from file frequency.dat.'
    else
        write(*,*) 'The frequency is deterimed by FFTW.'
    endif
endif

!-------------------------------------frequency----------------------------------

if(Omega_cw=='T') then

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

    allocate(q0(Ln))
    q0=cmplx(0.0d0,0.0d0)
    call dfftw_plan_dft_1d(plan,Ln,modal_velocity,q0,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,modal_velocity,q0)
    call dfftw_destroy_plan(plan)

    q0(:int(Ln/2))=real(q0(:int(Ln/2)))**2+imag(q0(:int(Ln/2)))**2
    f(1)=maxloc(real(q0(:int(Ln/2))),1)
    omega_target=float(f(1))/dt/Ln
    deallocate(q0)

endif

if (my_id==root_proc) then
    write(*,'(A23,f15.3,A5)') ' The center frequency is',omega_target, 'THz.'
    write(*,*) ' '
    write(*,*) '------------------------------------------------------------'
endif

if(omega_target<0.0 .or. omega_target==0.0) then
if (my_id==root_proc) then
write(*,*) 'The center frequency is not reasonable. Stop!'
endif
stop
endif

call MPI_SCATTERV(time_id(:),sendcount,displs,MPI_DOUBLE_PRECISION,time_number_node(:),recvcount,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)

!-------------------------------------calculations-------------------------------------------

if (my_id==root_proc) then
    write(*,*) '######## 3. Setup for the phonon wavelet transform..'
    write(*,*) '------------------------------------------------------------'
    write(*,*) ' '
endif

Min_width=1.0/omega_target*sc_factor1
Max_width=time(Ln)*tau_up                     ! considering the maximum coherence time
width_number=sc_factor2
Step_width=(Max_width-Min_width)/float(width_number)
allocate(width_coh(width_number))
do i=1,width_number
width_coh(i)=Min_width+float(i-1)*Step_width
enddo

if (my_id==root_proc) then
    write(*,'(A27,f15.3,A4)') ' The minimun coherence time:',Min_width, 'ps.'
    write(*,'(A27,f15.3,A4)') ' The maximun coherence time:',Max_width, 'ps.'
    write(*,'(A36,I9,A2)') ' The time evolution steps in Wavelet:',num_time,'.'
    write(*,'(A26,I9,A2)') 'The coherence time steps:',width_number,'.'
    if(Spectral_cw=='F') then
        write(*,'(A40)') ' The calculation of spec_coh_t_N: Flase.'
    else
    write(*,*) 'The calculation of spec_coh_t_N: Ture.'
    endif
    write(*,*) ' '
    write(*,*) '------------------------------------------------------------'
    write(*,*) "######## 4. The calculation is begin.."

endif

if(Max_width*2.0*Size_WP>time(Ln)) then
if (my_id==root_proc) then
write(*,*) 'The coherence time is too large. Stop!'
endif
stop
endif

!------------------------------------------------------------------

call MPI_Barrier(MPI_COMM_WORLD,IERROR)

if (my_id==root_proc) then
    write(*,*) '------------------------------------------------------------'
    write(*,*) 'mpirun calculating...'
    write(*,*) '------------------------------------------------------------'
    write(*,*) ' '
endif

m=int(Max_width*2.0*Size_WP/dt)+1
allocate(motherwavelet(m))

allocate(wavelet_results_node(width_number,sendcount(my_id+1)))
allocate(coherence_dist(width_number,num_time))

if(Spectral_cw=='F') then

    do i=1,sendcount(my_id+1)
        
        if(my_id==root_proc) then
            write(*,'(A20,2I7)') 'In calculating step',i,sendcount(my_id+1)
        endif

        !write(*,*) '--------',0.0,'%  ----------'
        displace=int(time_number_node(i))

        do j=1,width_number
     
            width=width_coh(j)
     
            Min_time=int(maxval((/time(displace)-Size_WP*width,0.0d0/))/dt)+1
            Max_time=int(minval((/time(displace)+Size_WP*width,time(Ln)/))/dt)+1
     
            Step_time=1    !int(interval)
            ntime=Max_time-Min_time+1
            
            motherwavelet=cmplx(0.0d0,0.0d0);wavelet_results=cmplx(0.0d0,0.0d0)

            motherwavelet(:ntime)=exp(cmplx(0.0,2.0*pi*omega_target*(time(Min_time:Max_time)-time(displace))))*exp(-(time(Min_time:Max_time)-time(displace))**2/width**2/2.0)

            motherwavelet=motherwavelet/sqrt(sum(real(motherwavelet(:ntime))**2)+sum(imag(motherwavelet(:ntime))**2)+1.0e-20)
            wavelet_results=motherwavelet(:ntime)*modal_velocity(Min_time:Max_time)

            wavelet_results2=real(sum(wavelet_results))**2+imag(sum(wavelet_results))**2
            wavelet_results2=wavelet_results2/sqrt(width)

            wavelet_results_node(j,i)=wavelet_results2

            enddo

    enddo

    call MPI_Barrier(MPI_COMM_WORLD,IERROR)
    !-------------------------------------------------------------------------

    call MPI_GATHERV(wavelet_results_node,width_number*recvcount,MPI_DOUBLE_PRECISION,coherence_dist,width_number*sendcount,width_number*displs,MPI_DOUBLE_PRECISION,&
                        root_proc,MPI_COMM_WORLD,ierr)


    deallocate(motherwavelet)
    deallocate(wavelet_results_node)

else

    nomega=10
    allocate(omega_range(nomega))
    do i=1,nomega
        omega_range(i)=omega_target-nomega/2.0*0.1+0.1*i
    enddo
     
    allocate(coh_spectral(nomega,num_time))
    allocate(w_t_coh(nomega,width_number,num_time))
    allocate(coherence_dist0(width_number,num_time))
    
    do iw=1,nomega
    
        if(my_id==root_proc) then
            write(*,'(A20,2I7)') 'In calculating omega',iw,nomega
        endif

        do i=1,sendcount(my_id+1)
            
            
            !write(*,*) '--------',0.0,'%  ----------'
            displace=int(time_number_node(i))

            do j=1,width_number
         
                width=width_coh(j)
         
                Min_time=int(maxval((/time(displace)-Size_WP*width,0.0d0/))/dt)+1
                Max_time=int(minval((/time(displace)+Size_WP*width,time(Ln)/))/dt)+1
         
                Step_time=1    !int(interval)
                ntime=Max_time-Min_time+1
                
                motherwavelet=cmplx(0.0d0,0.0d0);wavelet_results=cmplx(0.0d0,0.0d0)

                motherwavelet(:ntime)=exp(cmplx(0.0,2.0*pi*omega_range(iw)*(time(Min_time:Max_time)-time(displace))))*exp(-(time(Min_time:Max_time)-time(displace))**2/width**2/2.0)

                motherwavelet=motherwavelet/sqrt(sum(real(motherwavelet(:ntime))**2)+sum(imag(motherwavelet(:ntime))**2)+1.0e-20)
                wavelet_results=motherwavelet(:ntime)*modal_velocity(Min_time:Max_time)

                wavelet_results2=real(sum(wavelet_results))**2+imag(sum(wavelet_results))**2
                wavelet_results2=wavelet_results2/sqrt(width)

                wavelet_results_node(j,i)=wavelet_results2

                enddo

        enddo

        call MPI_Barrier(MPI_COMM_WORLD,IERROR)
        !-------------------------------------------------------------------------

        call MPI_GATHERV(wavelet_results_node,width_number*recvcount,MPI_DOUBLE_PRECISION,coherence_dist0,width_number*sendcount,width_number*displs,MPI_DOUBLE_PRECISION,&
                            root_proc,MPI_COMM_WORLD,ierr)
        
        if(iw==int(nomega/2.0)) then
            coherence_dist=coherence_dist0
        endif
        
        coh_spectral(iw,:)=sum(coherence_dist0,dim=2)
        w_t_coh(iw,:,:)=coherence_dist0

    enddo

    deallocate(motherwavelet)
    deallocate(wavelet_results_node)
    deallocate(coherence_dist0)

endif

if(my_id==root_proc) then

    width_coh=width_coh*(2.0*sqrt(2*log(2.0)))
    coherence_dist=coherence_dist/hbar/omega_target

    open(1,file='coherence_time.dat')
    write(1,*) 't_evo    tau_c    Ph_number'

    do i=1,num_time
        do j=1,width_number
            write(1,*) time(int(time_id(i))),width_coh(j),coherence_dist(j,i)
        enddo
    enddo
    close(1)

    allocate(coh_density(width_number))
    coh_density=0.0d0
    do i=1,num_time
        coh_density=coh_density+coherence_dist(:,i)
    enddo
    coh_density=coh_density/sum(coh_density)

    open(1,file='coh_density.dat')
    write(1,*) 'tau_c   TAPND'
    do j=1,width_number
        write(1,*)  width_coh(j),coh_density(j)
    enddo
    close(1)

    if(Spectral_cw=='T') then

        open(1,file='Spectral_evo.dat')
        write(1,*) 'omega   t_evo   Ph_number'
        do i=1,nomega
            do j=1,num_time
                write(1,*)  omega_range(i),time(int(time_id(j))),coh_spectral(i,j)
            enddo
        enddo
        close(1)
     
        open(1,file='w_t_coh.dat')
        write(1,*) 'omega   t_evo   tau_c   Ph_number'
        do iw=1,nomega
            do i=1,num_time
                do j=1,width_number
                    write(1,*) time(int(time_id(i))),width_coh(j),w_t_coh(iw,j,i)
                enddo
            enddo
            write(1,*)
        enddo
        close(1)

    endif

!------------------------------------------------------------------------
write(*,*) ' '
write(*,*) "PWave calculation is successfully done!"
write(*,*) '------------------------------------------------------------'
write(*,*) "######## 5. The output files:"
write(*,*) '------------------------------------------------------------'
write(*,*) ' '
write(*,'(A43)',advance='no') '1. coherence_time.dat,  2. coh_density.dat'
n=2
if(Spectral_cw=='T') then
n=n+1
write(*,'(I2,A19)',advance='no') n,'. Spectral_evo.dat'
n=n+1
write(*,'(I2,A11)',advance='no') n,'. w_t_coh.dat'
endif

write(*,*) '.'
write(*,*) ' '
write(*,*) '------------------------------------------------------------'
write(*,*) "Thank you for your using!"
write(*,*) ' '

call CPU_TIME(time1)
write(*,*) 'The used times:',time1-time0,'seconds.'


endif

call MPI_FINALIZE(ierr)

end program
