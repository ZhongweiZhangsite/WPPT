 !GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the phonon coherence from a parallel running code

module wl_trans_temp

   contains
   attributes(global)subroutine Cal_temp_1D(modal_velocity,width_coh,&
    omega_target,Size_WP,time,dt,Step_displace,wavelet_results_node)

implicit none
double complex,intent(in):: modal_velocity(:)
double precision,intent(in):: width_coh(:),omega_target
double precision,intent(in):: Size_WP,time(:),dt
integer,intent(in):: Step_displace

integer:: displace,a1,a2,width_number
integer:: j,ia
double precision:: width,temp1,temp3,temp4
double complex:: tmp,tmp1,tmp2,tmp3
double precision,intent(out):: wavelet_results_node(:,:)

displace=blockDim%x*(blockIdx%x-1)+threadIdx%x

width_number=size(width_coh)
width=maxval(width_coh)

a1=int(maxval((/(time(displace*Step_displace)-1.1*Size_WP*width)/dt,1.0d0/)))
a2=int(minval((/(time(displace*Step_displace)+1.1*Size_WP*width)/dt,maxval(time)/dt/)))

if(float(a2-a1)*dt>2.0*Size_WP*width) then

    do j=1,width_number
         
        width=width_coh(j)
         
        !a1=int(maxval((/(time(displace*Step_displace)-1.15*Size_WP*width)/dt,1.0d0/)))
        !a2=int(minval((/(time(displace*Step_displace)+1.15*Size_WP*width)/dt,maxval(time)/dt/)))

        temp1=0.0;tmp3=cmplx(0.0,0.0)
        
        do ia=a1,a2
              
            temp3=omega_target*(time(ia)-time(displace*Step_displace))
            temp4=(time(ia)-time(displace*Step_displace))**2
                            
            if(sqrt(temp4)<Size_WP*width) then

                tmp1=exp(cmplx(0.0,2.0*3.1415926*temp3))*exp(-temp4/(width)**2*4.0*log(2.0))
                tmp2=tmp1*modal_velocity(ia)
                    
                temp1=temp1+real(tmp1)**2+imag(tmp1)**2
                tmp3=tmp3+tmp2
                        
            endif

        enddo
        
        temp1=sqrt(temp1)
        wavelet_results_node(j,displace)=real(tmp3/temp1)**2+imag(tmp3/temp1)**2
        wavelet_results_node(j,displace)=wavelet_results_node(j,displace)/sqrt(width)
               
    enddo
         
else

    wavelet_results_node(:,displace)=0.0d0

endif

end subroutine Cal_temp_1D
end module wl_trans_temp
