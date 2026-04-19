   integer function opengl_proc()
   opengl_proc=2;end
   integer opengl_proc,window;external opengl_proc
   include <opengl.ins>,nolist
parameter(id=400,jd=300)
integer mg(0:id,0:jd+8),idm,jdm,iter,itr,Kshow
real u(0:id,0:jd),f(0:id,0:jd)
real h,Dt,Vm,cal,sigma,scale
 !ivan.platonychev@math.msu.ru / iv2017plat@gmail.com
DATA cal/.5/,sigma/0.5/,iter/1000/,Kshow/1/,scale/12./

   ier=winio@('%sp%ww[no_border]%pv%^og[double]%lw',0,0,id,jd+8,opengl_proc,window)

jdm=jd-1; idm=id-1 
h = 1./jd
mg=0;u=0.;f=0.01
call image;call image

1 CONTINUE;WRITE(*,*) ' 0-EXIT/1-EXE/2-u=0/3-cal/4-sigma/7-ITER/8-SCALE'
READ (*,*) key;SELECT CASE(key)
CASE(1)!MAIN PROGRAMM!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
   mg=0;do itr=1,iter
call ptm
call image
!read(*,*)
        enddo
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CASE(2); u=0.;call image
CASE(3); WRITE(*,*) cal,' - calcul parameter'; READ (*,*) cal
CASE(4); WRITE(*,*) sigma,' - calcul parameter'; READ (*,*) sigma
CASE(7); WRITE(*,*) iter; READ (*,*) iter
CASE(77);WRITE(*,*) Kshow; READ (*,*) Kshow
CASE(8); WRITE(*,*) scale,' - SCALE'; READ (*,*) scale;call image
CASE(0); GO TO 100       
END SELECT
GO TO 1

100   CONTINUE!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CONTAINs

subroutine PTM;real uu !Delta(u)+g=0, sigma=2*h*h/dt
uu(i,j)=(u(i,j)*sm2+u(i+1,j)+u(i,j+1)+u(i-1,j)+u(i,j-1)+f(i,j))*osp2
sm2=sigma-2.;osp2=1./(sigma+2.)
Do j = 1, jdm; do i = 1, idm; u(i,j)=uu(i,j); enddo; EndDo
Do j=jdm,1,-1; do i=idm,1,-1; u(i,j)=uu(i,j); enddo; EndDo 
Do j=jdm,1,-1; do i = 1, idm; u(i,j)=uu(i,j); enddo; EndDo
Do j = 1, jdm; do i=idm,1,-1; u(i,j)=uu(i,j); enddo; EndDo 
endsubroutine

subroutine image;kwh=127+ishft(127,8)+ishft(127,16)
Do j=0,jd; do i=0,id
kol=mod(nint(u(i,j)*scale),128)
!if(kol>0)then;mg(i,j)=kol;else;mg(i,j)=ishft(-kol,16);endif !Black Phone
if(kol>0)then;mg(i,j)=kwh-ishft(kol,8)-ishft(kol,16);
         else;mg(i,j)=kwh+kol+ishft(kol,8);endif !White Phone
enddo; EndDo
i=itr*id/iter;mg(i,jd+2:jd+7)=ishft(127,8)
call glDrawPixels(id+1,jd+9,GL_rgba,GL_byte,mg);call swap_opengl_buffers()
endsubroutine
END!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX