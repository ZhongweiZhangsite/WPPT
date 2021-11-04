program main
  
integer ii
double precision recording(3)

!----------------------------------------------------------------
  
open(1,file='velocity.out')
Open(2, file='vx.bin',access='direct',form='unformatted',recl=2)
Open(3, file='vy.bin',access='direct',form='unformatted',recl=2)
Open(4, file='vz.bin',access='direct',form='unformatted',recl=2)

ii=0
do while(.true.)		
 
	read(1,*,iostat=ierr) recording(:)
	
	if(ierr/=0)exit
	ii=ii+1	
	write(2,rec=ii) recording(1)
	write(3,rec=ii) recording(2)
	write(4,rec=ii) recording(3)

enddo
close(1)
close(2)
close(3)
close(4)

end program

