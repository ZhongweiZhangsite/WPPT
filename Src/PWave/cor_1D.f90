 !GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

!-----to calculate and write the phonon coherence from a parallel running code
  
module cor_transform

   contains
   attributes(global)subroutine Cal_cor1D(q1,cor_node)

implicit none
double precision,intent(in):: q1(:)

integer:: i,j,n,displace

double precision:: temp1
double precision,intent(out):: cor_node(:)
 
displace=blockDim%x*(blockIdx%x-1)+threadIdx%x

n=size(q1)
j=0;temp1=0.0

do i=1,n-displace
    j=j+1
    temp1=temp1+q1(i)*q1(i+displace)
enddo

if(j>0) then
    cor_node(displace)=temp1/float(j)
else
    cor_node(displace)=0.0d0
endif
 
end subroutine Cal_cor1D
end module cor_transform
