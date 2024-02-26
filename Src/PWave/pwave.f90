 !GNU GENERAL PUBLIC LICENSE
!                      Version 3, 29 June 2007

!Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
!Everyone is permitted to copy and distribute verbatim copies
!of this license document, but changing it is not allowed.

program PWave

use temporal
use spatial

implicit none
integer:: i
double precision:: time0,time1
character:: char0,calcul_mode

CHARACTER(256) :: folder_path
LOGICAL :: folder_exists

!----------------------------------------------------------------
call CPU_TIME(time0)

write(*,*) '------------------------------------------------------------------------'
write(*,*) '----------***----******---***---*******---*******--********-------------'
write(*,*) '-----------***---******---***--***---**--***---**-----**----------------'
write(*,*) '-----------***--***-***--***---********--********-----**----------------'
write(*,*) '------------***-***--***-***---***-------***----------**----------------'
write(*,*) '------------*******--******----***-------***----------**----------------'
write(*,*) '------------*******--******----**--------**-----------**---V.3.1--------'
write(*,*) '------------------------------------------------------------------------'
write(*,*) '-------The PWave module for coherence calculations in WPPT package------'
write(*,*) '------------------------------GPU version-------------------------------'
write(*,*) '---------------------------By Zhongwei Zhang----------------------------'
write(*,*) '---------------------Email: zhongwei@tongji.edu.cn----------------------'
write(*,*) '---------------------------Date: Dec. 20 2023---------------------------'
write(*,*) '------------------------------------------------------------------------'
 
!-----------------------parameters setting------------------------

folder_path = './CONTROL.WPPT'
INQUIRE(FILE=folder_path, EXIST=folder_exists)

IF(folder_exists) THEN
    open(1,file='CONTROL.WPPT')
    do i=1,15
    read(1,*)
    enddo
    read(1,*) char0,calcul_mode
    close(1)
ELSE
    write(*,*) "ERROR: The file of CONTROL.WPPT does not exist..."
    stop
ENDIF

write(*,*) "########################################################################"

if(calcul_mode=='T') then
    write(*,*) "Temporal coherence in calculating..."
else if(calcul_mode=='S') then
    write(*,*) "Spatial coherence in calculating..."
else
    write(*,*) "ERROR: The calculation mode is a bad choice..."
    stop
endif

write(*,*) "########################################################################"

if(calcul_mode=='T') then
    call temporal_coh
else
    call spatial_coh
endif

write(*,*) ' '

write(*,*) '------------------------------------------------------------'
write(*,*) "             Thank you for using WPPT package!"
write(*,*) ' '

call CPU_TIME(time1)

write(*,*) 'The used times:',time1-time0,'seconds.'

end program PWave
