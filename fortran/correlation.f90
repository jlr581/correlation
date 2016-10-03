program guass

use correlation
use bootstrap

implicit none

integer :: i,j,n,m
real, allocatable :: x(:),y(:),tx(:),ty(:)
real :: r(3),tmp1,tmp2
character*255 :: f1,f2

print *,'Citation'
print *,'Roberts et al. (2016), Correlation confidence limits for unevenly sampled data, Computers & Geosciences'
print *
print *,'Enter filename of first data series'
read *,f1
print *,'Enter filename of second data series'
read *,f2

open(11,file=trim(f1),form="formatted")
n=0
do 
  read (11,*,end=10,err=10)tmp1,tmp2
  n=n+1
enddo
10 close(11)

open(11,file=trim(f2),form="formatted")
m=0
do 
  read (11,*,end=20,err=20)tmp1,tmp2
  m=m+1
enddo
20 close(11)

allocate (x(n),tx(n),y(m),ty(m))

open(11,file=trim(f1),form="formatted")
do i=1,n
  read(11,*)tx(i),x(i)
enddo
close(11)

open(11,file=trim(f2),form="formatted")
do i=1,m
  read(11,*)ty(i),y(i)
enddo
close(11)

r(1)=correlate_gaussian(x(1:n),y(1:m),tx(1:n),ty(1:m))
r(2:3)=bootstrap_ci(x(1:n),y(1:m),tx(1:n),ty(1:m))
print *,r(1),'[',r(2),":",r(3),"]"

deallocate (x,tx,y,ty)



end program
