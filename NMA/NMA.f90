!-----to calculate and write the modal velocity from a parallel running code

program main
use mpi
include '/home/zhongwei/software/fftw/include/fftw3.f'

integer i,j,k,m,n,dimen
integer Ln,num,num_cell,i_v,n_omega,n_calculation                      !settings
integer n_sample,n_interval,nbins,nlest,recvcount
double precision,allocatable:: vector(:),q_position(:),q_mode(:,:),wavevector(:,:)
!double precision,allocatable:: modal_real(:,:),modal_imag(:,:),modal_ve_real(:,:),modal_ve_imag(:,:)
double precision,allocatable:: cell_position(:,:),sed_omega(:),sed(:),DOS(:),ve_x(:,:),ve_y(:,:),ve_z(:,:)
double precision:: eigen(2)
double complex::  temp1,temp2
double complex,allocatable:: eigenvector(:,:,:,:),q(:),q0(:)
real, parameter :: pi=3.1415926,mass=1.993e-26,Atom=1.0e-10,pstos=1.0e-12
character tmp 
integer	status(MPI_STATUS_SIZE)
integer(kind=MPI_OFFSET_KIND)	offset
integer  fh
integer,allocatable:: sendcount(:),displs(:)
character*16 infile1,infile2,infile3
character(len=10) :: file_id
character(len=50) :: file_name
character char0
character*4 NMAw_cw,SED_cw,DOS_cw         ! The control characters to calculate (c) and write (w) properties

integer ierr,my_id,num_procs,root_proc,filetype
integer :: plan

!----------------------------------------------------------------

integer :: tx
character( Len = 85 ) :: cStr

call CPU_TIME(time0)
    
root_proc=0
    
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
 
if(my_id==root_proc)then

    write(*,*) '------------------------------------------------------------'
    write(*,*) '---------Time-dependent modal velocity for phonons----------'
    write(*,*) '---------------The NMA module in WPPT package---------------'
    write(*,*) '---------------------By Zhongwei Zhang----------------------'
    write(*,*) '---------------Email: zhongwei@tongji.edu.cn----------------'
    write(*,*) '--------------------Date: Oct. 11 2020---------------------'
    write(*,*) '------------------------------------------------------------'

endif
!-----------------------all nodes setting------------------------

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,Ln
read(1,*) char0,dt
read(1,*)
read(1,*)
read(1,*) char0,num
read(1,*) char0,num_cell
read(1,*) char0,dimen
read(1,*) char0,n_k
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*) char0,NMAw_cw
read(1,*) char0,SED_cw
read(1,*) char0,DOS_cw
close(1)

!--------------------------read periodicity--------------------------

open(1,file='periodicity.dat')
read(1,*) num1

if(my_id==root_proc)then
if(num1.ne.num) then
write(*,*)"The number of atoms inconsistent in CONTROL and periodicity.dat"
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)
stop
endif
endif

allocate(cell_position(num,4))           ! See the format in manual
do i=1,num   
	read(1,*) cell_position(i,:)		 
end do
close(1)
 	  
!---------------------------read eigenvector----------------------

num_cell1=int(maxval(cell_position(:,4)))

if(my_id==root_proc)then
if(num_cell1.ne.num_cell) then
write(*,*)"The number of primitive cell atoms inconsistent in CONTROL and periodicity.dat"
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)
stop
endif
endif

n_f=num_cell*3
allocate(eigenvector(n_k,n_f,num_cell,3))

open(1,file='Eigenvector')
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
read(1,*) n_k1

if(my_id==root_proc)then
if(n_k1.ne.n_k) then
write(*,*)"The number of Kpoints inconsistent in CONTROL and periodicity.dat"
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)
stop
endif
endif

allocate(wavevector(3,n_k))
wavevector=0.0
do i=1,n_k
    read(1,*) wavevector(:,i)
enddo
close(1)

if (my_id==root_proc) then
  
    write(*,*) ' '
    write(*,*) "-----------Successfully read periodicity and eigenvector files!"
    write(*,*) "------------------------------------------------------------"
    write(*,*) ' '
    write(*,*) "--------------------------------------------------"
    write(*,*) 'The atom number in MD simulation:',num,'primitive cell:',num_cell
    write(*,*) 'The number of step and timestep:',Ln,dt,'ps'
    write(*,*) "*******************"
    write(*,*) 'The supercell number in MD:',int(maxval(cell_position(:,1))),int(maxval(cell_position(:,2))),int(maxval(cell_position(:,3)))
    write(*,*) 'The wavevector number:',n_k
    write(*,*) "*******************"
    write(*,*) "--------------------------------------------------"
    
    write(*,*) "------------------The calculation is begin:"

endif

!----------------------parallel message sending------------------

n_calculation=n_k*n_f

allocate(sendcount(num_procs))
allocate(displs(num_procs)) 

nbins=n_calculation/num_procs
nlest=mod(n_calculation,num_procs)

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

!-------------------------read modal velicity ------------------------

if(dimen==1) then

    CALL getarg(1,infile1)

    allocate(ve_x(num,Ln))

    if (my_id==root_proc) then
    write(*,*) 'mpirun reading..'
    endif

    call MPI_Barrier(MPI_COMM_WORLD,IERROR)

    call MPI_File_open(MPI_COMM_WORLD,infile1,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

    call MPI_File_read_all(fh,ve_x,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    !call MPI_File_read_ordered(fh,ve_id,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    call MPI_File_close(fh, ierror)

else

    CALL getarg(1,infile1)
    CALL getarg(2,infile2)
    CALL getarg(3,infile3)

    allocate(ve_x(num,Ln))
    allocate(ve_y(num,Ln))
    allocate(ve_z(num,Ln))

    if (my_id==root_proc) then
    write(*,*) 'mpirun reading..'
    endif

    call MPI_Barrier(MPI_COMM_WORLD,IERROR)

    call MPI_File_open(MPI_COMM_WORLD,infile1,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

    call MPI_File_read_all(fh,ve_x,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    !call MPI_File_read_ordered(fh,ve_id,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    call MPI_File_close(fh, ierror)

    call MPI_File_open(MPI_COMM_WORLD,infile2,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

    call MPI_File_read_all(fh,ve_y,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    !call MPI_File_read_ordered(fh,ve_id,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    call MPI_File_close(fh, ierror)

    call MPI_File_open(MPI_COMM_WORLD,infile3,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

    call MPI_File_read_all(fh,ve_z,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    !call MPI_File_read_ordered(fh,ve_id,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
    call MPI_File_close(fh, ierror)

endif

if (my_id==root_proc) then
	write(*,*) '------------------------------------------------------------'
	write(*,*) 'Finish the atomic velocity reading..'

    allocate(q_position(n_k*n_f))
    q_position=0.0

    k=0
    do i=1,n_k
        do j=1,n_f
            k=k+1
            q_position(k)=float(k)
        enddo
    enddo
 
    write(*,*) "--------------------------------------------------"
    write(*,*) "------------------The calculation is begin:"

endif

!------------------------------------------------------------------

call MPI_Barrier(MPI_COMM_WORLD,IERROR)
call CPU_TIME(time1)

if (my_id==root_proc) then
	write(*,*) '------------------------------------------------------------'
	write(*,*) 'mpirun calculating..'
	write(*,*) '------------------------------------------------------------'
!write(*,*) "The total time cost is: ", time1-time0, " s."
endif

! Distribute the matrix A into different processes and store in the local matrix A1 

call MPI_SCATTERV(q_position(:),sendcount,displs,MPI_DOUBLE_PRECISION,vector(:),recvcount,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)


!--------------------------calculation module---------------------------------------

allocate(q(Ln))
allocate(q0(Ln))
allocate(sed_omega(Ln))

do i=1,sendcount(my_id+1)

!--------------------open files--------------------------

	m=int(q_mode(int(vector(i)),1))
	n=int(q_mode(int(vector(i)),2))
    
    q=cmplx(0.0d0,0.0d0)
 
	do i_step=1, Ln				  ! time step

		do j_num=1, num			  ! all atoms position  
  
			temp1=cmplx(0.0d0,wavevector(1,m)*cell_position(j_num,1)+wavevector(2,m)*cell_position(j_num,2)+wavevector(3,m)*cell_position(j_num,3))
					
			temp2=exp(2.0*pi*temp1)
					
			k=int(cell_position(j_num,4))
            
                if(dimen==1)then
                    q(i_step)=q(i_step)+ve_x(j_num,i_step)*temp2*(eigenvector(m,n,k,1))
                else
                    q(i_step)=q(i_step)+ve_x(j_num,i_step)*temp2*(eigenvector(m,n,k,1))&
                        +ve_y(j_num,i_step)*temp2*(eigenvector(m,n,k,2))&
                        +ve_z(j_num,i_step)*temp2*(eigenvector(m,n,k,3))
                endif

		enddo					 ! all position
   
	end do                               !time step

    write(file_id, '(i0,A1,i0)') m,'_',n
    file_name = 'NMAt_' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(file_name))

    write(1,*) 't-Modal velocity.'
    write(1,*) 'wavevector:',m,'Branch:',n
    do i_step=1, Ln
        write(1,*) float(i_step)*dt,real(q(i_step))/sqrt(float(num)),imag(q(i_step))/sqrt(float(num))
    enddo
  	close(1)

    if(NMAw_cw=='Ture') then
        
        call dfftw_plan_dft_1d(plan,Ln,q,q0,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,q,q0)
        call dfftw_destroy_plan(plan)
        
        sed_omega=real(q0(:))**2+imag(q0(:))**2
        
        file_name = 'NMAw_' // trim(adjustl(file_id)) // '.dat'
        open(1,file = trim(file_name))

        write(1,*) 'SED for each mode can used for lifetime calculation'
        write(1,*) 'wavevector:',m,'Branch:',n
        do i_step=1, int(Ln/2)
            write(1,*) 1.0/dt/float(Ln)*float(i_step),sed_omega(i_step)
        enddo
        close(1)

    endif

enddo

call MPI_Barrier(MPI_COMM_WORLD,IERROR)

if(my_id==root_proc) then

n_omega=int(Ln/2.0)

if(SED_cw=='Ture') then

    open(2,file='SED_k_w.dat')

    allocate(sed(n_omega))
    do i=1,n_k
        sed=0.0
        do j=1,n_f

            write(file_id, '(i0,A1,i0)') i,'_',j
            file_name = 'NMAt_' // trim(adjustl(file_id)) // '.dat'
            open(1,file = trim(file_name))
            read(1,*)
            read(1,*)
            do i_step=1, n_omega
                write(1,*) char0,sed_omega(i_step)
            enddo
            close(1)
            sed=sed+sed_omega

        enddo
        
        write(2,'(5E20.10)') wavevector(:,i),1.0/dt/float(Ln)*float(i_step),sed(i_step)
    enddo
    close(2)
endif

if(DOS_cw=='Ture') then

    open(2,file='dos.dat')

    allocate(DOS(n_omega))
    DOS=0.0
    do i=1,n_k
        do j=1,n_f

            write(file_id, '(i0,A1,i0)') i,'_',j
            file_name = 'NMAt_' // trim(adjustl(file_id)) // '.dat'
            open(1,file = trim(file_name))
            read(1,*)
            read(1,*)
            do i_step=1, n_omega
                write(1,*) char0,sed_omega(i_step)
            enddo
            close(1)
            DOS=DOS+sed_omega

        enddo
        write(2,'(2E20.10)') 1.0/dt/float(Ln)*float(i_step),DOS(i_step)
    enddo
    close(2)
endif

!------------------------------------------------------------------------
write(*,*)
write(*,*) "NMA calculation is successfully done!"
write(*,*) "Thank you for your using!"

endif

call MPI_FINALIZE(ierr)

end program
