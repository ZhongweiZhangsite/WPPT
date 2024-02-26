!GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the modal velocity from a parallel running code

program main

use mpi
include 'fftw3.f'

integer:: i,j,k,m,n,i_sample,n_sample,x
integer:: Ln,Ln2,n_omega,n_calculation,n_interval,i_step,n1,n2
integer:: num,num_cell,j_num,num1,num_cell1,dimen
integer:: i_v,n_k,n_k1,n_f
integer:: nbins,nlest,recvcount,npieces,ipiece
double precision,allocatable:: vector(:),q_position(:),q_mode(:,:),wavevector(:,:)
double precision,allocatable:: cell_position(:,:),ve_x(:,:),ve_y(:,:),ve_z(:,:)
double precision,allocatable:: sed_omega(:),sed(:),omega(:),DOS(:)
double precision:: eigen(2),para,dt,time0,time1
double complex::  temp1,temp2
double complex,allocatable:: eigenvector(:,:,:,:),q(:,:),q0(:),q1(:)
real, parameter :: pi=3.1415926,mass=1.66054e-27,Atom=1.0e-10,pstos=1.0e-12,ev2J=1.60218e-19
character::  calcul_mode
integer	status(MPI_STATUS_SIZE)
integer(kind=MPI_OFFSET_KIND)	offset
integer  fh,temporal_coh,spatial_coh,endprog
integer,allocatable:: sendcount(:),displs(:)
character*16 infile1,infile2,infile3
character(len=10) :: file_id
character(len=50) :: file_name
character char0
character*1 NMAw_cw,SED_cw,DOS_cw         ! The control characters to calculate (c) and write (w) properties

CHARACTER(256) :: folder_path
LOGICAL :: folder_exists

integer ierr,myid,numprocs,root_proc,filetype
integer :: plan

!----------------------------------------------------------------

integer :: tx
character( Len = 85 ) :: cStr

call CPU_TIME(time0)
    
root_proc=0
    
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
 
if(myid==root_proc)then
    
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) '----------***----******---***---*******---*******--********-------------'
    write(*,*) '-----------***---******---***--***---**--***---**-----**----------------'
    write(*,*) '-----------***--***-***--***---********--********-----**----------------'
    write(*,*) '------------***-***--***-***---***-------***----------**----------------'
    write(*,*) '------------*******--******----***-------***----------**----------------'
    write(*,*) '------------*******--******----**--------**-----------**---V.3.1--------'
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) '----------The NMA module for modal calculations in WPPT package---------'
    write(*,*) '---------------------------By Zhongwei Zhang----------------------------'
    write(*,*) '---------------------Email: zhongwei@tongji.edu.cn----------------------'
    write(*,*) '---------------------------Date: Dec. 20 2023---------------------------'
    write(*,*) '------------------------------------------------------------------------'
    
endif

!-----------------------all nodes setting------------------------

folder_path = './CONTROL.WPPT'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists) THEN
        
    open(1,file='CONTROL.WPPT')
    do i=1,3
    read(1,*)
    enddo
    read(1,*) char0,num_cell
    read(1,*) char0,n_k
    do i=1,2
    read(1,*)
    enddo
    read(1,*) char0,calcul_mode
    close(1)

ELSE

    if(myid==root_proc)then
        write(*,*)"ERROR: The file of CONTROL does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop

END IF

!----------------------parallel message setting and sending------------------

n_f=3*num_cell
n_calculation=n_k*n_f

allocate(sendcount(numprocs))
allocate(displs(numprocs))

nbins=n_calculation/numprocs
nlest=mod(n_calculation,numprocs)

do i=1,numprocs

    if(i>nlest) then
        sendcount(i)=nbins
    else
        sendcount(i)=nbins+1
    endif

enddo

displs(1)=0
do i=2,numprocs
     displs(i)=displs(i-1)+sendcount(i-1)
     displs(i)=displs(i)
enddo

do i=1,numprocs
    
    if(myid==(i-1)) then
        allocate(vector(sendcount(i)))
        vector(:)=0.0
        recvcount=sendcount(i)
    endif
enddo

allocate(q_mode(n_k*n_f,2))
q_mode=0.0
k=0
do i=1,n_k
    do j=1,n_f
        k=k+1
        q_mode(k,1)=float(i)
        q_mode(k,2)=float(j)
    enddo
enddo

if (myid==root_proc) then
   
    allocate(q_position(n_k*n_f))
    q_position=0.0

    k=0
    do i=1,n_k
        do j=1,n_f
            k=k+1
            q_position(k)=float(k)
        enddo
    enddo
 
endif

!------------------------make choice and goto calculation blocks----------------

if(calcul_mode=='T') then

    !goto temporal_coh
    if(myid==root_proc) then
    write(*,*) '********Temporal normal mode projection is set to be performed**********'
    write(*,*) '------------------------------------------------------------------------'
    endif
    goto 100

elseif(calcul_mode=='S') then

    !goto spatial_coh
    if(myid==root_proc) then
    write(*,*) '*********Spatial normal mode projection is set to be performed**********'
    write(*,*) '------------------------------------------------------------------------'
    endif
    goto 200

else

    if(myid==root_proc) then
    write(*,*)"ERROR: Made a bad choice..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop

endif

!--------------------------------- calculation block-I ----------------------------

100 continue

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,Ln
read(1,*) char0,num
do i=1,5
read(1,*)
enddo
read(1,*) char0,dt

read(1,*) char0,NMAw_cw
read(1,*) char0,SED_cw
read(1,*) char0,DOS_cw
read(1,*) char0,n_sample
close(1)

!----------------------------------files checking-------------------------------

folder_path = './periodicity.dat'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists) THEN
    open(1,file='periodicity.dat')
    read(1,*) num1

    if(num1.ne.num) then
    if(myid==root_proc) then
    write(*,*)"ERROR: The number of atoms is inconsistent between CONTROL and periodicity.dat..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
    endif
    close(1)
ELSE
    if(myid==root_proc) then
    write(*,*)"ERROR: The file of periodicity.dat does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

open(1,file='periodicity.dat')
read(1,*)
num_cell1=1
do i=1,100
read(1,*) para,para,para,para
if(para>num_cell1) then
num_cell1=int(para)
endif
enddo
close(1)
if(num_cell1.ne.num_cell) then
    if(myid==root_proc) then
    write(*,*)"ERROR: The number of atoms in primitive cell atoms inconsistent..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
endif

folder_path = './Eigenvector.dat'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists==.false.) THEN
    if(myid==root_proc) then
    write(*,*)"ERROR: The file of Eigenvector.dat does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

folder_path = './Kpoints.dat'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists) THEN
    open(1,file='Kpoints.dat')
    read(1,*) n_k1

    if(n_k1.ne.n_k) then
    if(myid==root_proc) then
    write(*,*)"ERROR: The number of Kpoints is inconsistent between CONTROL and Kpoints.dat..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
    endif
    close(1)
ELSE
    if(myid==root_proc) then
    write(*,*)"ERROR: The file of Kpoints.dat does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

folder_path = './vx.bin'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists==.false.) THEN
    if(myid==root_proc) then
    write(*,*)"ERROR: The files of v*.bin do not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

if(myid==root_proc) then

    folder_path = './NMAt'
    INQUIRE(DIRECTORY=folder_path, EXIST=folder_exists)
     
    IF (folder_exists) THEN
        call system("rm -fr NMAt/*")
    ELSE
        call system("mkdir NMAt")
    END IF
 
    if(NMAw_cw=='T') then

        folder_path = './NMAw'
        INQUIRE(DIRECTORY=folder_path, EXIST=folder_exists)
         
        IF (folder_exists) THEN
            call system("rm -fr NMAw/*")
        ELSE
            call system("mkdir NMAw")
        END IF

    endif

    if(SED_cw=='T') then

        folder_path = './SED'
        INQUIRE(DIRECTORY=folder_path, EXIST=folder_exists)
         
        IF (folder_exists) THEN
            call system("rm -fr SED/*")
        ELSE
            call system("mkdir SED")
        END IF

    endif

    if(DOS_cw=='T') then

        folder_path = './DOS'
        INQUIRE(DIRECTORY=folder_path, EXIST=folder_exists)
         
        IF (folder_exists) THEN
            call system("rm -fr DOS/*")
        ELSE
            call system("mkdir DOS")
        END IF

    endif
    write(*,*) ''
    write(*,*) "######## I. Files checked!"
    write(*,*) " "
    write(*,*) '------------------------------------------------------------------------'

endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)

npieces=int(Ln/int(huge(i)/num)*10+1)
if(npieces<10) then
    npieces=10  !default
endif
Ln2=int(Ln/npieces)

if(myid==root_proc)then
    write(*,*) " "
    write(*,*) 'Number of cpus: ',numprocs
    write(*,*) 'Number of block:',npieces
endif

!---------------------------------files reading-------------------------------

open(1,file='periodicity.dat')
read(1,*)
allocate(cell_position(num,5))           ! See the format in manual
do i=1,num
    read(1,*) cell_position(i,:)
end do
close(1)

n_f=num_cell*3
allocate(eigenvector(n_k,n_f,num_cell,3))
open(1,file='Eigenvector.dat')
do i=1,n_k
    do j=1,n_f
        do k=1,num_cell
            do x=1,3
                read(1,*) eigen(:)
                eigenvector(i,j,k,x)=cmplx(eigen(1),-eigen(2))
            enddo
        enddo
    enddo
enddo
close(1)
 
open(1,file='Kpoints.dat')
read(1,*)
allocate(wavevector(3,n_k))
wavevector=0.0
do i=1,n_k
    read(1,*) wavevector(:,i)
enddo
close(1)
 
if (myid==root_proc) then
  
    write(*,*) ' '
    write(*,*) 'The atoms in MD inputs:',num,'; in primitive cell:',num_cell,'.'
    write(*,*) 'The nstep,timestep:',Ln,',',dt,'ps.'
    write(*,'(A20,3I5,A3)') 'The supercell size:',int(maxval(cell_position(:,1))),int(maxval(cell_position(:,2))),int(maxval(cell_position(:,3))),'.'
    write(*,*) 'Wavevector number:',n_k,'; branch number:',n_f,'.'
    write(*,*) ' '
    write(*,*) "######## II. Pre-prepared files have been successfully read"
    write(*,*) "########       and settings are configured!"
    write(*,*) " "
    write(*,*) '------------------------------------------------------------------------'
    write(*,*) " "
endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)
!------------------------------------------------------------------

if (myid==root_proc) then
    write(*,*) "######## III. The calculation is beginning..."
endif

allocate(ve_x(num,Ln2))
allocate(ve_y(num,Ln2))
allocate(ve_z(num,Ln2))
 
Ln=npieces*Ln2
allocate(q(sendcount(myid+1),Ln))
q=cmplx(0.0d0,0.0d0)
 
call MPI_SCATTERV(q_position(:),sendcount,displs,MPI_DOUBLE_PRECISION,vector(:),recvcount,&
                        MPI_DOUBLE_PRECISION,root_proc,MPI_COMM_WORLD,ierr)

do ipiece=1,npieces

        if(myid==root_proc) then
           write(*,'(A31,2I7)') '          In calculating block',ipiece,npieces
        endif

        call MPI_Barrier(MPI_COMM_WORLD,IERR)
        offset=num*Ln2*(ipiece-1)*8
        call MPI_File_open(MPI_COMM_WORLD,'vx.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_File_read_at_all(fh,offset,ve_x,num*Ln2, MPI_DOUBLE_PRECISION, status,IERR)
        call MPI_File_close(fh, IERR)

        call MPI_File_open(MPI_COMM_WORLD,'vy.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_File_read_at_all(fh,offset,ve_y,num*Ln2, MPI_DOUBLE_PRECISION, status,IERR)
        call MPI_File_close(fh, IERR)

        call MPI_File_open(MPI_COMM_WORLD,'vz.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_File_read_at_all(fh,offset,ve_z,num*Ln2, MPI_DOUBLE_PRECISION, status,IERR)
        call MPI_File_close(fh, IERR)

        do i=1,sendcount(myid+1)
 
            m=int(q_mode(int(vector(i)),1))
            n=int(q_mode(int(vector(i)),2))
             
            do i_step=1, Ln2				  ! time step

                do j_num=1, num			  ! all atoms position
          
                    temp1=cmplx(0.0d0,wavevector(1,m)*cell_position(j_num,1)+wavevector(2,m)*cell_position(j_num,2)+wavevector(3,m)*cell_position(j_num,3))
                            
                    temp2=exp(2.0*pi*temp1)
                            
                    k=int(cell_position(j_num,4))
                      
                    q(i,(ipiece-1)*Ln2+i_step)=q(i,(ipiece-1)*Ln2+i_step)+&
                                ve_x(j_num,i_step)*temp2*(eigenvector(m,n,k,1))*sqrt(cell_position(j_num,5))&
                                +ve_y(j_num,i_step)*temp2*(eigenvector(m,n,k,2))*sqrt(cell_position(j_num,5))&
                                +ve_z(j_num,i_step)*temp2*(eigenvector(m,n,k,3))*sqrt(cell_position(j_num,5))
           
                enddo					  ! all position
           
            end do                        !time step

        enddo

        call MPI_Barrier(MPI_COMM_WORLD,ierr)

enddo

q=q/sqrt(float(num/num_cell))

deallocate(eigenvector)
deallocate(ve_x)
deallocate(ve_y)
deallocate(ve_z)

do i=1,sendcount(myid+1)
 
    m=int(q_mode(int(vector(i)),1))
    n=int(q_mode(int(vector(i)),2))
        
    write(file_id, '(i0,A1,i0)') m,'_',n
    file_name = 'NMAt/NMAt_' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(file_name))

    write(1,*) 't-Modal velocity.'
    write(1,*) 'wavevector:',m,'Branch:',n
    write(1,*) 't_evo     real    imag'
    do i_step=1, Ln
        write(1,*) float(i_step)*dt,real(q(i,i_step)),imag(q(i,i_step))
    enddo
  	close(1)

enddo

n_interval=int(Ln/5.0)
n_omega=int(n_interval/2.0)

allocate(q0(n_interval))
allocate(q1(n_interval))
allocate(sed_omega(n_omega))

if(NMAw_cw=='T') then

    do i=1,sendcount(myid+1)

        sed_omega=0.0
        do i_sample=1,n_sample

            n1=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+1
            n2=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+n_interval

            q0(:)=q(i,n1:n2);q0(:)=q0(:)-sum(q0(:))/float(n_interval)
            q1=cmplx(0.0d0,0.0d0)

            call dfftw_plan_dft_1d(plan,n_interval,q0,q1,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute_dft(plan,q0,q1)
            call dfftw_destroy_plan(plan)
            
            sed_omega=sed_omega+real(q1(1:n_omega))**2/float(n_interval)+imag(q1(1:n_omega))**2/float(n_interval)
             
        enddo

        m=int(q_mode(int(vector(i)),1))
        n=int(q_mode(int(vector(i)),2))
            
        write(file_id, '(i0,A1,i0)') m,'_',n

        file_name = 'NMAw/NMAw_' // trim(adjustl(file_id)) // '.dat'
        open(1,file = trim(file_name))

        write(1,*) 'SED for each mode can used for lifetime calculation'
        write(1,*) 'wavevector:',m,'Branch:',n
        write(1,*) 'omega     SED'
        do i_step=1, n_omega
            write(1,*) 1.0/dt/float(n_interval)*float(i_step),sed_omega(i_step)
        enddo
        close(1)

    enddo

endif
 
call MPI_Barrier(MPI_COMM_WORLD,ierr)
 
deallocate(q)
deallocate(q0)
deallocate(q1)
 
if(myid==root_proc) then

    if(SED_cw=='T') then

        open(21,file='SED/SED_k_w.dat')

        allocate(omega(n_omega))
        allocate(sed(n_omega))

        write(21,*) 'kx    ky    kz    omega    SED'

        do i=1,n_k

            sed=0.0
            do j=1,n_f

                write(file_id, '(i0,A1,i0)') i,'_',j
                file_name = 'NMAw/NMAw_' // trim(adjustl(file_id)) // '.dat'
                open(1,file = trim(file_name))
                read(1,*)
                read(1,*)
                read(1,*)
                do i_step=1, n_omega
                    read(1,*) omega(i_step),sed_omega(i_step)
                enddo
                close(1)

                sed=sed+sed_omega

            enddo

            do i_step=1, n_omega
                write(21,'(5E20.10)') wavevector(:,i),omega(i_step),sed(i_step)
            enddo
 
        enddo

        close(21)
        deallocate(omega)
        deallocate(sed)

    endif

    if(DOS_cw=='T') then
 
        allocate(DOS(n_omega))
        allocate(omega(n_omega))
        DOS=0.0
        do i=1,n_k
            do j=1,n_f

                write(file_id, '(i0,A1,i0)') i,'_',j
                file_name = 'NMAw/NMAw_' // trim(adjustl(file_id)) // '.dat'
                open(1,file = trim(file_name))
                read(1,*)
                read(1,*)
                read(1,*)
                do i_step=1, n_omega
                    read(1,*) omega(i_step),sed_omega(i_step)
                enddo
                close(1)
                DOS=DOS+sed_omega

            enddo
        enddo

        open(22,file='DOS/DOS_w.dat')
        write(22,*) 'omega    DOS'
        do i_step=1, n_omega
        write(22,'(2E20.10)') omega(i_step),DOS(i_step)
        enddo
        close(22)
        deallocate(DOS)
        deallocate(omega)

    endif

!------------------------------------------------------------------------
write(*,*) "Temporal NMA calculation is successfully done!"
write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
write(*,*) "######## IV. The output files:"
write(*,'(A17)',advance='no') '1. NMAt_*_*.dat'
n=1
if(NMAw_cw=='T') then
n=n+1
write(*,'(I2,A14)',advance='no') n,'. NMAw_*_*.dat'
endif
if(SED_cw=='T') then
n=n+1
write(*,'(I2,A13)',advance='no') n,'. SED_k_w.dat'
endif
if(DOS_cw=='T') then
n=n+1
write(*,'(I2,A11)',advance='no') n,'. DOS_w.dat'
endif
write(*,*) '.'
endif

goto 300

!--------------------------------- calculation block-II ----------------------------

200 continue

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,Ln
read(1,*) char0,num
read(1,*) char0,num_cell
read(1,*) char0,n_k
close(1)

!----------------------------------files checking-------------------------------

folder_path = './periodicity.dat'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists) THEN
    open(1,file='periodicity.dat')
    read(1,*) num1

    if(num1.ne.num) then
    if(myid==root_proc) then
    write(*,*)"ERROR: The number of atoms is inconsistent between CONTROL and periodicity.dat..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
    endif
    close(1)
ELSE
    if(myid==root_proc) then
    write(*,*)"ERROR: The file of periodicity.dat does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

open(1,file='periodicity.dat')
read(1,*)
num_cell1=1
do i=1,100
read(1,*) para,para,para,para
if(para>num_cell1) then
num_cell1=int(para)
endif
enddo
close(1)

if(num_cell1.ne.num_cell) then
    if(myid==root_proc) then
    write(*,*)"ERROR: The number of atoms in primitive cell atoms inconsistent..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
endif

folder_path = './Eigenvector.dat'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists==.false.) THEN
    if(myid==root_proc) then
    write(*,*)"ERROR: The file of Eigenvector.dat does not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF

folder_path = './vx.bin'
INQUIRE(FILE=folder_path, EXIST=folder_exists)
 
IF(folder_exists==.false.) THEN
    if(myid==root_proc) then
    write(*,*)"ERROR: The files of v*.bin do not exist..."
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
END IF
 
if (myid==root_proc) then

    folder_path = './Spatial'
    INQUIRE(DIRECTORY=folder_path, EXIST=folder_exists)
     
    IF (folder_exists) THEN
        call system("rm -fr Spatial/*")
    ELSE
        call system("mkdir Spatial")
    END IF
    write(*,*) ''
    write(*,*) "######## I. Files checked!"
    write(*,*) " "
    write(*,*) '------------------------------------------------------------------------'

endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)

!---------------------------------files reading-------------------------------

open(1,file='periodicity.dat')
read(1,*)
allocate(cell_position(num,5))           ! See the format in manual
do i=1,num
    read(1,*) cell_position(i,:)
end do
close(1)

n_f=num_cell*3
allocate(eigenvector(n_k,n_f,num_cell,3))
open(1,file='Eigenvector.dat')
do i=1,n_k
    do j=1,n_f
        do k=1,num_cell
            do x=1,3
                read(1,*) eigen(:)
                eigenvector(i,j,k,x)=cmplx(eigen(1),-eigen(2))
            enddo
        enddo
    enddo
enddo
close(1)
 
if (myid==root_proc) then
  
    write(*,*) ' '
    write(*,*) 'The atoms in MD inputs:',num,'; in primitive cell:',num_cell,'.'
    write(*,*) 'The nstep,timestep:',Ln,',',dt,'ps.'
    write(*,'(A20,3I5,A3)') 'The supercell size:',int(maxval(cell_position(:,1))),int(maxval(cell_position(:,2))),int(maxval(cell_position(:,3))),'.'
    write(*,*) 'Wavevector number:',n_k,'; branch number:',n_f,'.'
    write(*,*) ' '
    write(*,*) "######## II. Pre-prepared files have been successfully read"
    write(*,*) "########       and settings are configured!"
    write(*,*) ' '
    write(*,*) '------------------------------------------------------------------------'

endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_SCATTERV(q_position(:),sendcount,displs,MPI_DOUBLE_PRECISION,vector(:),recvcount,&
                        MPI_DOUBLE_PRECISION,root_proc,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------------------------------------

if (myid==root_proc) then
    write(*,*) ' '
    write(*,*) "######## III. The calculation is beginning..."
endif
 
!---------------------------------------------------------------------------------------------

allocate(ve_x(num,Ln))
allocate(ve_y(num,Ln))
allocate(ve_z(num,Ln))

allocate(q(1,1))
q=cmplx(0.0d0,0.0d0)
  
call MPI_Barrier(MPI_COMM_WORLD,IERR)
call MPI_File_open(MPI_COMM_WORLD,'vx.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_x,num*Ln, MPI_DOUBLE_PRECISION, status,ierr)
call MPI_File_close(fh, ierr)

call MPI_File_open(MPI_COMM_WORLD,'vy.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_y,num*Ln, MPI_DOUBLE_PRECISION, status,ierr)
call MPI_File_close(fh, ierr)

call MPI_File_open(MPI_COMM_WORLD,'vz.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_z,num*Ln, MPI_DOUBLE_PRECISION, status,ierr)
call MPI_File_close(fh, ierr)
  
do i=1,sendcount(myid+1)

    if(myid==root_proc) then
        write(*,'(A31,2I7)') '          In calculating block',i,sendcount(myid+1)
    endif

    m=int(q_mode(int(vector(i)),1))
    n=int(q_mode(int(vector(i)),2))
    
    write(file_id, '(i0,A1,i0)') m,'_',n
    file_name = 'Spatial/Spatial_' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(file_name))

    write(1,*) 'Spatial-time Modal velocity.'
    write(1,'(A11,I4,A7,I5,I8)') 'wavevector:',m,'Branch:',n,num
    write(1,*) 'ai real imag'

    do i_step=1, Ln                  ! time step

        do j_num=1, num              ! all atoms position
      
            k=int(cell_position(j_num,4))
                
            q(1,1)= ve_x(j_num,i_step)*(eigenvector(m,n,k,1))*sqrt(cell_position(j_num,5))&
                          +ve_y(j_num,i_step)*(eigenvector(m,n,k,2))*sqrt(cell_position(j_num,5))&
                          +ve_z(j_num,i_step)*(eigenvector(m,n,k,3))*sqrt(cell_position(j_num,5))

            write(1,'(I8,2E15.7)') j_num,real(q(1,1)),imag(q(1,1))

        enddo                      ! all position

        write(1,*)

    end do                         !time step
    close(1)

enddo

deallocate(ve_x)
deallocate(ve_y)
deallocate(ve_z)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

!------------------------------------------------------------------------
if(myid==root_proc) then


write(*,*) "NMA_spatial calculation is successfully done!"
write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
write(*,*) "######## IV. The output files:"
write(*,'(A19)',advance='no') '1. Spatial_*_*.dat'
write(*,*) '.'

endif

300 continue

if(myid==root_proc) then

write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
write(*,*) "Thank you for your using!"

call CPU_TIME(time1)
write(*,*) 'The used times:',time1-time0,'seconds.'

endif

call MPI_FINALIZE(ierr)
 
end program
