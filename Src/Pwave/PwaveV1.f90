program main

use mpi
include '/home/zhongwei/software//fftw/include/fftw3.f'

integer i,id,ii,jj,m,xx,num,Ln,num_time,width_number,f(1)
integer Min_time,Max_time,Step_time,displace,Step_displace
double precision,allocatable:: time_number_node(:),time_id(:)
double precision pos(2),omega_target,wavelet_results1,width,dt
double precision tau_up,sc_factor1,sc_factor2,Size_WP
double precision Min_width,Max_width,Step_width,Min_displace,Max_displace,real_time,process,wavelet_results2
double complex::  temp1,temp2
double precision,allocatable::time(:),width_coh(:),coherence_dist(:,:),wavelet_results_node(:,:),coh_density(:)
double complex,allocatable::modal_velocity(:),q0(:),wavelet_results(:),motherwavelet(:)
integer :: plan
real, parameter :: pi=3.1415926,mass=1.993e-26,Atom=1.0e-10,pstos=1.0e-12,hbar=6.62607015e-22
character tmp

integer(kind=MPI_OFFSET_KIND)    offset
integer  fh,nbins,nlest,recvcount
integer,allocatable:: sendcount(:),displs(:)
character*16 infile1,infile2,infile3
character char0
character Omega_cw
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
close(1)

!-----------------------------------Initialization-----------------------------------------

allocate(modal_velocity(Ln))
allocate(time(Ln))

Open(1,file=infile1)
read(1,*)
read(1,*) char0,wvector,char0,branch
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
    write(*,*) 'The nstep,timestep:',Ln,',',dt,'ps'
    write(*,*) 'The calculated step:',num_time
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
    do i=1,branch
    read(1,'(f)',advance='no') omega_target
    enddo
    read(1,*) omega_target
    close(1)

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
    write(*,*) 'The center frequency is',omega_target, 'THz.'
    write(*,*) ' '
    write(*,*) '------------------------------------------------------------'
endif

if(omega_target<0.0 .or. omega_target==0.0) then
write(*,*) 'The center frequency is not reasonable. Stop.'
stop
endif

call MPI_SCATTERV(time_id(:),sendcount,displs,MPI_DOUBLE_PRECISION,time_number_node(:),recvcount,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)

!-------------------------------------calculations-------------------------------------------

if (my_id==root_proc) then
    write(*,*) '######## 3. Setup for the phonon wavelet transform...'
    write(*,*) '------------------------------------------------------------'
    write(*,*) ' '
endif

Min_width=1.0/omega_target*sc_factor2
Max_width=time(Ln)*tau_up                     ! considering the maximum coherence time
width_number=sc_factor1
Step_width=(Max_width-Min_width)/float(width_number)
allocate(width_coh(width_number))
do i=1,width_number
width_coh(i)=Min_width+float(i-1)*Step_width
enddo

if (my_id==root_proc) then
    write(*,*) 'The minimun coherence time:',Min_width, 'ps.'
    write(*,*) 'The maximun coherence time:',Max_width, 'ps.'
    write(*,*) 'The time evolution steps in Wavelet:',num_time,'.'
    write(*,*) 'The coherence time steps:',width_number,'.'
    write(*,*) ' '
    write(*,*) '------------------------------------------------------------'
    write(*,*) "######## 4. The calculation is begin..."

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

if(my_id==root_proc) then

width_coh=width_coh*(2.0*sqrt(2*log(2.0)))
coherence_dist=coherence_dist/hbar/omega_target

open(1,file='coherence_time.dat')
write(1,*) 't_evo  tau_c  Ph_number'

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
write(1,*) 'tau_c TAPND'
do j=1,width_number
    write(1,*)  width_coh(j),coh_density(j)
enddo
close(1)

!------------------------------------------------------------------------
write(*,*) ' '
write(*,*) '------------------------------------------------------------'
write(*,*) "PWave calculation is successfully done!"
write(*,*) "######## 5. The output files:"
write(*,*) '------------------------------------------------------------'
write(*,*) ' '
write(*,'(A17)',advance='no') '1. coherence_time.dat,  2. coh_density.dat'
n=2
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
