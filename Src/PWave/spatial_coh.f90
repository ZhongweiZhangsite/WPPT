!GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the phonon coherence from a parallel running code
 
module spatial

   contains
   subroutine spatial_coh

use wl_trans_spatial
use cudafor

implicit none
integer:: i,j,k,Ln,Ln2,num_time,width_number,ndimen
integer:: Min_time,Max_time,Step_time,displace,Step_displace
integer:: nk,n_k,num_cell,n_f,na,nb,nc,n,ik,ib,nstep
integer:: n1,n2
integer,allocatable:: do_id(:)

double precision,allocatable:: wavevector(:,:),cell_position(:,:),cell_period(:,:),mass_id(:)
double precision:: width,dt,pos(2)
double precision:: time0,time1,x,rate
double precision:: vector_target(3),lc_up,lc_down,sc_factor1,Size_WP
double precision:: Min_width,Max_width,Step_width
double precision,allocatable:: wavelet_results_node(:,:),coh_density(:),width_coh(:)
double complex,allocatable::modal_velocity(:)

double precision,device:: vector_targetd(3)
double precision,device,allocatable:: modal_velocityd(:),cell_positiond(:,:),cell_periodd(:,:),width_cohd(:)
double precision,device,allocatable:: wavelet_results_noded(:,:)
integer,device,allocatable:: do_idd(:)

character char0,Spatial_cw,SCD_WC
character(len=10) :: file_id
character(len=50) :: file_name

integer::tPB=32
type(dim3) :: grid, tBlock
logical isExist

!-----------------------parameters setting------------------------

open(1,file='CONTROL.WPPT')
read(1,*)
read(1,*) char0,nstep
read(1,*) char0,Ln
read(1,*) char0,num_cell
read(1,*) char0,n_k
n_f=num_cell*3
do i=1,11
read(1,*)
enddo
read(1,*) char0,rate
read(1,*) char0,lc_up
read(1,*) char0,lc_down
read(1,*) char0,width_number
read(1,*) char0,Size_WP
read(1,*)
read(1,*)
read(1,*) char0,SCD_WC
read(1,*)
read(1,*)
read(1,*) char0,ndimen
close(1) 

!-----------------------------------Initialization-----------------------------------------
 
inquire(file='Kpoints.dat',exist=isExist)

if(isExist) then

    open(1,file='Kpoints.dat')
    read(1,*) n_k

    allocate(wavevector(3,n_k))
    wavevector=0.0
    do i=1,n_k
        read(1,*) wavevector(:,i)
    enddo
    close(1)

else
    write(*,*) "ERROR: The file of Kpoints.dat does not exist..."
    stop
endif

inquire(file='model.xyz',exist=isExist)

if(isExist) then

    Open(1,file='model.xyz')
    allocate(cell_position(Ln,3))           ! See the format in manual
    allocate(cell_positiond(Ln,3))
    read(1,*)
    read(1,*)
    do i=1,Ln
        read(1,*) char0,cell_position(i,1:3)
    enddo
    close(1)

else
 
    inquire(file='position.dat',exist=isExist)
    if(isExist) then

    Open(1,file='position.dat')
    allocate(cell_position(Ln,3))           ! See the format in manual
    allocate(cell_positiond(Ln,3))
    do i=1,Ln
        read(1,*) j,j,cell_position(i,1:3)
    enddo
    close(1)

    else
    
    write(*,*) "ERROR: The file of position.dat does not exist..."
    stop
    endif
endif

cell_position(:,1)=cell_position(:,1)-minval(cell_position(:,1))
cell_position(:,2)=cell_position(:,2)-minval(cell_position(:,2))
cell_position(:,3)=cell_position(:,3)-minval(cell_position(:,3))
cell_positiond=cell_position

inquire(file='periodicity.dat',exist=isExist)

if(isExist) then

    Open(1,file='periodicity.dat')
    
    allocate(mass_id(Ln))
    allocate(cell_period(Ln,3))           ! See the format in manual
    allocate(cell_periodd(Ln,3))
    
    read(1,*)
    do i=1,Ln
        read(1,*) cell_period(i,1:3),x,mass_id(i)
    enddo
    close(1)

else
 
    write(*,*) "ERROR: The file of periodicity.dat does not exist..."
    stop
    
endif

cell_periodd=cell_period

allocate(do_id(Ln))           ! See the format in manual
allocate(do_idd(Ln))
do_id=0

call init_random_seed()

Ln2=0
do i=1,Ln
     
    call random_number(x)
    if(x<(1.0/rate) .or. x==(1.0/rate)) then
        Ln2=Ln2+1
        do_id(Ln2)=i
    endif

enddo

do_idd=do_id

!----------------------parallel message sending------------------

write(*,*) "######## I. Successfully read pre-prepared files..."
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
write(*,'(A18,I9,A1,A17,I9,A1)') ' The input natoms:',Ln,',',' The used natoms:',Ln2,'.'
write(*,'(A23,I4,A11,I4,A1)') ' Number of wavevectors:',n_k,'; branches:',n_f,'.'
write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'

write(*,*) '######## II. Setup for the phonon wavelet transform...'
write(*,*) '------------------------------------------------------------------------'
write(*,*) ' '
!-------------------------------------frequency/vector----------------------------------

Min_width=lc_down !/(2.0*sqrt(2*log(2.0)))      !1.0/maxval(vector_target(:2))*sc_factor1
Max_width=lc_up !/(2.0*sqrt(2*log(2.0)))        !1.0/maxval(vector_target(:2))*lc_up  ! considering the maximum coherence time
Step_width=(Max_width-Min_width)/float(width_number)

allocate(width_coh(width_number))
allocate(width_cohd(width_number))

do i=1,width_number
    width_coh(i)=Min_width+float(i-1)*Step_width 
enddo

width_cohd=width_coh

if(Max_width*2.0*Size_WP>maxval(cell_position(:,:))) then
write(*,*) 'ERROR: The coherence length is too long...'
stop
endif

write(*,'(A30,f13.3,A4)') ' The minimun coherence length:',lc_down, ' Angstrom.'
write(*,'(A30,f13.3,A4)') ' The maximun coherence length:',lc_up, ' Angstrom.'
write(*,'(A21,I7,A2)') ' The coherence steps:',width_number,'.'

write(*,*) ' '
write(*,*) '------------------------------------------------------------------------'
write(*,*) "######## III. GPU parallel calculation is beginning..."
write(*,*) ' '

!---------------------------------files checking--------------------------------

if(SCD_WC=='T') then

    inquire(file='./Sp_C_D',exist=isExist)
    if(isExist) then
        call system("rm -fr Sp_C_D/*")
    else
        call system("mkdir Sp_C_D")
    endif
else

    inquire(file='./Sp_C_D',exist=isExist)
    if(isExist) then
        call system("rm -fr Sp_C_D")
    endif
endif

inquire(file='./LenC_Dens',exist=isExist)
if(isExist) then
    call system("rm -fr LenC_Dens/*")
else
    call system("mkdir LenC_Dens")
endif
 
!----------------------------allocation---------------------------------

tBlock = dim3(tPB,1,1)
grid = dim3(int(real(Ln2)/tBlock%x),1,1)

Ln2=tBlock%x*grid%x

allocate(wavelet_results_node(width_number,Ln2))
allocate(wavelet_results_noded(width_number,Ln2))
wavelet_results_node=0.0d0
wavelet_results_noded=0.0

allocate(modal_velocity(Ln))
allocate(modal_velocityd(Ln))
allocate(coh_density(width_number))

!---------------------------------------------------------

if(ndimen==1) then
  
    do ik=1,n_k
        
        write(*,'(A30,2I7)') '    In calculating Kpoints', ik, n_k

        vector_target=wavevector(:,ik)
        vector_targetd=vector_target

        do ib=1,n_f

            !----------------------------------------read modal velocity from file------------------------------

            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'Spatial/Spatial_' // trim(adjustl(file_id)) // '.dat'

            inquire(file=trim(file_name),exist=isExist)
            if(isExist) then
                open(1,file=trim(file_name))
                do i=1,3
                read(1,*)
                enddo
                do i=1,(nstep-1)*(Ln+1)
                read(1,*)
                enddo

                do i=1,Ln
                    read(1,*) j,pos(:)
                    modal_velocity(i)=cmplx(pos(1),pos(2))
                enddo
                close(1)
                modal_velocityd=modal_velocity
            else
                write(*,*) "ERROR: The files of Spatial/Spatial_* do not exist..."
                stop
            endif

            !--------------------------------------------1D transform----------------------------------------

            call Cal_spatial_1D<<<grid, tBlock>>>(modal_velocityd,cell_positiond,cell_periodd,width_cohd,vector_targetd,do_idd,wavelet_results_noded)

            wavelet_results_node=wavelet_results_noded

            if(SCD_WC=='T') then

                write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
                file_name = 'Sp_C_D/spa_len_' // trim(adjustl(file_id)) // '.dat'

                open(1,file=trim(file_name))
                write(1,*) 'x  y  z  l_c  Ph_number'

                do i=1,Ln2
                    do j=1,width_number
                        write(1,'(5E15.7)') cell_position(do_id(i),1),cell_position(do_id(i),2),cell_position(do_id(i),3),width_coh(j),wavelet_results_node(j,i)
                    enddo
                enddo

                close(1)
        
            endif

            coh_density=0.0d0
            do i=1,Ln2
                coh_density=coh_density+wavelet_results_node(:,i)
            enddo
            coh_density=coh_density/sum(coh_density)

            write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
            file_name = 'LenC_Dens/coh_' // trim(adjustl(file_id)) // '.dat'
            
            open(1,file=trim(file_name))
            write(1,*) 'l_c  SAPND'
            do j=1,width_number
                write(1,*)  width_coh(j),coh_density(j)
            enddo
            close(1)

        enddo
    enddo
 
elseif(ndimen==2) then

    do ik=1,n_k
        
        write(*,'(A27,2I7)') ' In calculating Kpoints', ik, n_k

        vector_target=wavevector(:,ik)
        vector_targetd=vector_target

        do ib=1,n_f

!----------------------------------------read modal velocity from file------------------------------

            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'Spatial/Spatial_' // trim(adjustl(file_id)) // '.dat'
            open(1,file=trim(file_name))
            do i=1,3
            read(1,*)
            enddo

            do i=1,(nstep-1)*(Ln+1)
            read(1,*)
            enddo

            do i=1,Ln
                read(1,*) j,pos(:)
                modal_velocity(i)=cmplx(pos(1),pos(2))
            enddo
            close(1)
            modal_velocityd=modal_velocity

!--------------------------------------------2D transform----------------------------------------

            call Cal_spatial_2D<<<grid, tBlock>>>(modal_velocityd,cell_positiond,cell_periodd,width_cohd,vector_targetd,do_idd,wavelet_results_noded)

            wavelet_results_node=wavelet_results_noded
            
            if(SCD_WC=='T') then

                write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
                file_name = 'Sp_C_D/spa_len' // trim(adjustl(file_id)) // '.dat'
     
                open(1,file=trim(file_name))
                write(1,*) 'x  y  z  l_c  Ph_number'

                    do i=1,Ln2
                        do j=1,width_number
                            write(1,'(5E15.7)') cell_position(do_id(i),1),cell_position(do_id(i),2),cell_position(do_id(i),3),width_coh(j),wavelet_results_node(j,i)
                        enddo
                    enddo
                 
                close(1)

            endif

            coh_density=0.0d0
            do i=1,Ln2
                coh_density=coh_density+wavelet_results_node(:,i)
            enddo
            coh_density=coh_density/sum(coh_density)
            
            write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
            file_name = 'LenC_Dens/coh_' // trim(adjustl(file_id)) // '.dat'
            
            open(1,file=trim(file_name))
            write(1,*) 'l_c   SAPND'
            do j=1,width_number
                write(1,*)  width_coh(j),coh_density(j)
            enddo
            close(1)

        enddo
    enddo

elseif(ndimen==3) then

    do ik=1,n_k
        
        write(*,'(A27,2I7)') ' In calculating Kpoints', ik, n_k

        vector_target=wavevector(:,ik)
        vector_targetd=vector_target

        do ib=1,n_f

!----------------------------------------read modal velocity from file------------------------------

            write(file_id, '(i0,A1,i0)') ik,'_',ib
            file_name = 'Spatial/Spatial_' // trim(adjustl(file_id)) // '.dat'
            open(1,file=trim(file_name))
            do i=1,3
            read(1,*)
            enddo

            do i=1,(nstep-1)*(Ln+1)
            read(1,*)
            enddo

            do i=1,Ln
                read(1,*) j,pos(:)
                modal_velocity(i)=cmplx(pos(1),pos(2))
            enddo
            close(1)
            modal_velocityd=modal_velocity

!--------------------------------------------3D transform----------------------------------------

            call Cal_spatial_3D<<<grid, tBlock>>>(modal_velocityd,cell_positiond,cell_periodd,width_cohd,vector_targetd,do_idd,wavelet_results_noded)

            wavelet_results_node=wavelet_results_noded

            if(SCD_WC=='T') then

                write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
                file_name = 'Sp_C_D/spa_len' // trim(adjustl(file_id)) // '.dat'

                open(1,file=trim(file_name))
                write(1,*) 'x  y  z  l_c  Ph_number'

                    do i=1,Ln2
                        do j=1,width_number
                            write(1,'(5E15.7)') cell_position(do_id(i),1),cell_position(do_id(i),2),cell_position(do_id(i),3),width_coh(j),wavelet_results_node(j,i)
                        enddo
                    enddo
                close(1)

            endif

            coh_density=0.0d0
            do i=1,Ln2
                coh_density=coh_density+wavelet_results_node(:,i)
            enddo
            coh_density=coh_density/sum(coh_density)
            
            write(file_id, '(i0,A1,i0,A1,i0)') ik,'_',ib,'_',nstep
            file_name = 'LenC_Dens/coh_' // trim(adjustl(file_id)) // '.dat'
            
            open(1,file=trim(file_name))
            write(1,*) 'l_c   SAPND'
            do j=1,width_number
                write(1,*)  width_coh(j),coh_density(j)
            enddo
            close(1)

        enddo
    enddo

endif

!------------------------------------------------------------------------
 write(*,*) ' '
 write(*,*) "Coherence calculation is successfully done!"
 write(*,*) ' '
 write(*,*) '------------------------------------------------------------------------'
 write(*,*) "######## IV. The output files:"
 write(*,'(A23)',advance='no') ' 1. LenC_Dens/coh_*.dat'
 n=1

if(SCD_WC=='T') then
     n=n+1
     write(*,'(I2,A22)',advance='no') n,'. Sp_C_D/spa_len_*.dat'
endif

write(*,*) '.'

end subroutine spatial_coh

subroutine init_random_seed()
   integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
   
   call random_seed(size = n)
   allocate(seed(n))
   
   call system_clock(count=clock)
   
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(put = seed)
   
   deallocate(seed)
end subroutine init_random_seed

end module spatial
