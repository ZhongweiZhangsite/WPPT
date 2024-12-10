!GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the phonon coherence from a parallel running code


module wl_trans_spatial

   contains  

attributes(global)subroutine Cal_spatial_1D(modal_velocity,cell_position,cell_period,width_coh,&
     vector_target,do_id,wavelet_results_node)
     
implicit none

double precision,intent(in):: modal_velocity(:),cell_position(:,:),cell_period(:,:)
double precision,intent(in):: width_coh(:),vector_target(:)
integer,intent(in):: do_id(:)

integer:: displace0,displace,width_number
integer:: j,ia,Ln
double precision:: width,temp1,temp3,temp4,x1(3),x2(3),a1,a2
double precision,parameter:: Size_WP=2.0
double complex:: tmp1,tmp2,tmp3
double precision,intent(out):: wavelet_results_node(:,:)

displace0=blockDim%x*(blockIdx%x-1)+threadIdx%x   ! position, time
displace=do_id(displace0)

width_number=size(width_coh)
width=maxval(width_coh)
  
x1(:)=cell_position(displace,:)

a1=maxval((/cell_position(displace,1)-1.1*Size_WP*width,0.0d0/))
a2=minval((/cell_position(displace,1)+1.1*Size_WP*width,maxval(cell_position(:,1))/))
 
Ln=size(cell_position(:,1))

if((a2-a1)>2.0*Size_WP*width) then

  do j=1,width_number

      width=width_coh(j)

      temp1=0.0;tmp3=cmplx(0.0,0.0)
      
      do ia=1,Ln
           
          temp3=vector_target(1)*(cell_period(ia,1)-cell_period(displace,1))+&
                vector_target(2)*(cell_period(ia,2)-cell_period(displace,2))+&
                vector_target(3)*(cell_period(ia,3)-cell_period(displace,3))
                          
          x2(:)=cell_position(ia,:)

          temp4=(x2(1)-x1(1))**2
                          
          if(sqrt(temp4)<Size_WP*width) then

                  tmp1=exp(cmplx(0.0,2.0*3.1415926*temp3))*exp(-temp4/(width)**2*4.0*log(2.0))
                  tmp2=tmp1*modal_velocity(ia)
                      
                  temp1=temp1+real(tmp1)**2+imag(tmp1)**2
                  tmp3=tmp3+tmp2
                      
          endif
      
       enddo
      
      temp1=sqrt(temp1)
      
      wavelet_results_node(j,displace0)=real(tmp3/temp1)**2+imag(tmp3/temp1)**2
      wavelet_results_node(j,displace0)=wavelet_results_node(j,displace0)/sqrt(width)
             
  enddo
  
else

     wavelet_results_node(:,displace0)=0.0d0

endif


end subroutine Cal_spatial_1D


attributes(global)subroutine Cal_spatial_2D(modal_velocity,cell_position,cell_period,width_coh,&
     vector_target,do_id,wavelet_results_node)
     
implicit none

double precision,intent(in):: modal_velocity(:),cell_position(:,:),cell_period(:,:)
double precision,intent(in):: width_coh(:),vector_target(:)
integer,intent(in):: do_id(:)

integer:: displace0,displace,width_number
integer:: j,ia,Ln
double precision:: width,temp1,temp3,temp4,x1(3),x2(3),a1,a2,b1,b2
double precision,parameter:: Size_WP=2.0
double complex:: tmp1,tmp2,tmp3
double precision,intent(out):: wavelet_results_node(:,:)

displace0=blockDim%x*(blockIdx%x-1)+threadIdx%x   ! position, time
displace=do_id(displace0)

width_number=size(width_coh)
width=maxval(width_coh)
  
x1(:)=cell_position(displace,:)

a1=maxval((/cell_position(displace,1)-1.1*Size_WP*width,0.0d0/))
a2=minval((/cell_position(displace,1)+1.1*Size_WP*width,maxval(cell_position(:,1))/))

b1=maxval((/cell_position(displace,2)-1.1*Size_WP*width,0.0d0/))
b2=minval((/cell_position(displace,2)+1.1*Size_WP*width,maxval(cell_position(:,2))/))

Ln=size(cell_position(:,1))

if((a2-a1)>2.0*Size_WP*width .and. (b2-b1)>2.0*Size_WP*width) then

  do j=1,width_number

      width=width_coh(j)

      temp1=0.0;tmp3=cmplx(0.0,0.0)
      
      do ia=1,Ln
           
          temp3=vector_target(1)*(cell_period(ia,1)-cell_period(displace,1))+&
                vector_target(2)*(cell_period(ia,2)-cell_period(displace,2))+&
                vector_target(3)*(cell_period(ia,3)-cell_period(displace,3))
                          
          x2(:)=cell_position(ia,:)

          temp4=(x2(1)-x1(1))**2+(x2(2)-x1(2))**2
                          
          if(sqrt(temp4)<Size_WP*width) then

                  tmp1=exp(cmplx(0.0,2.0*3.1415926*temp3))*exp(-temp4/(width)**2*4.0*log(2.0))
                  tmp2=tmp1*modal_velocity(ia)
                      
                  temp1=temp1+real(tmp1)**2+imag(tmp1)**2
                  tmp3=tmp3+tmp2
                      
          endif
      
       enddo
      
      temp1=sqrt(temp1)
      
      wavelet_results_node(j,displace0)=real(tmp3/temp1)**2+imag(tmp3/temp1)**2
      wavelet_results_node(j,displace0)=wavelet_results_node(j,displace0)/sqrt(width)
             
  enddo
  
else

     wavelet_results_node(:,displace0)=0.0d0

endif

end subroutine Cal_spatial_2D


attributes(global)subroutine Cal_spatial_3D(modal_velocity,cell_position,cell_period,width_coh,&
     vector_target,do_id,wavelet_results_node)
     
implicit none

double precision,intent(in):: modal_velocity(:),cell_position(:,:),cell_period(:,:)
double precision,intent(in):: width_coh(:),vector_target(:)
integer,intent(in):: do_id(:)

integer:: displace0,displace,width_number
integer:: j,ia,Ln
double precision:: width,temp1,temp3,temp4,x1(3),x2(3),a1,a2,b1,b2,c1,c2
double precision,parameter:: Size_WP=2.0
double complex:: tmp1,tmp2,tmp3
double precision,intent(out):: wavelet_results_node(:,:)

displace0=blockDim%x*(blockIdx%x-1)+threadIdx%x   ! position, time
displace=do_id(displace0)

width_number=size(width_coh)
width=maxval(width_coh)
  
x1(:)=cell_position(displace,:)

a1=maxval((/cell_position(displace,1)-1.1*Size_WP*width,0.0d0/))
a2=minval((/cell_position(displace,1)+1.1*Size_WP*width,maxval(cell_position(:,1))/))

b1=maxval((/cell_position(displace,2)-1.1*Size_WP*width,0.0d0/))
b2=minval((/cell_position(displace,2)+1.1*Size_WP*width,maxval(cell_position(:,2))/))

c1=maxval((/cell_position(displace,2)-1.1*Size_WP*width,0.0d0/))
c2=minval((/cell_position(displace,2)+1.1*Size_WP*width,maxval(cell_position(:,2))/))

Ln=size(cell_position(:,1))

if((a2-a1)>2.0*Size_WP*width .and. (b2-b1)>2.0*Size_WP*width .and. (c2-c1)>2.0*Size_WP*width) then

  do j=1,width_number

      width=width_coh(j)

      temp1=0.0;tmp3=cmplx(0.0,0.0)
      
      do ia=1,Ln
           
          temp3=vector_target(1)*(cell_period(ia,1)-cell_period(displace,1))+&
                vector_target(2)*(cell_period(ia,2)-cell_period(displace,2))+&
                vector_target(3)*(cell_period(ia,3)-cell_period(displace,3))
                          
          x2(:)=cell_position(ia,:)

          temp4=(x2(1)-x1(1))**2+(x2(2)-x1(2))**2+(x2(3)-x1(3))**2
                          
          if(sqrt(temp4)<Size_WP*width) then

                  tmp1=exp(cmplx(0.0,2.0*3.1415926*temp3))*exp(-temp4/(width)**2*4.0*log(2.0))
                  tmp2=tmp1*modal_velocity(ia)
                      
                  temp1=temp1+real(tmp1)**2+imag(tmp1)**2
                  tmp3=tmp3+tmp2
                      
          endif
      
       enddo
      
      temp1=sqrt(temp1)
      
      wavelet_results_node(j,displace0)=real(tmp3/temp1)**2+imag(tmp3/temp1)**2
      wavelet_results_node(j,displace0)=wavelet_results_node(j,displace0)/sqrt(width)
             
  enddo
  
else

     wavelet_results_node(:,displace0)=0.0d0

endif

end subroutine Cal_spatial_3D

end module wl_trans_spatial
