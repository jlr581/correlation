module correlation

! JLR 30/6/2014
!
! uneven sampled correlation using Gaussian kernel using method of
! Rehfeld et al, Nonlinear Processes in Geophysics, 2011

implicit none

contains

function correlate_gaussian(x,y,tx,ty)
real :: correlate_gaussian

real, intent(in) :: x(:),y(:),tx(:),ty(:)
integer :: n,m,i,j
real :: b,h,d,xmean,ymean,txmean,tymean,delta_t,pi
real :: sdx(size(x)),sdy(size(x)),num(size(x))
real :: sdxg,sdyg,numg,deng

n=size(x,1)
m=size(y,1)

pi=abs(atan2(0.0,-1.0))


! calc mean for x series
xmean=0.0
do i=1,n
  xmean=xmean+x(i)
enddo
xmean=xmean/n
txmean=(tx(n)-tx(1))/(n-1)

! calc mean for y series
ymean=0.0
do j=1,m
  ymean=ymean+y(j)
enddo
ymean=ymean/m
tymean=(ty(m)-ty(1))/(m-1)

sdx=0.0
sdy=0.0
num=0.0
delta_t=max(txmean,tymean)
h=delta_t/4

!$acc kernels loop gang, vector(64)
do i=1,n
  do j=1,m
    d=ty(j)-tx(i)
    b=exp(-d**2/(2*h**2))/sqrt(2*pi*h)
    num(i)=num(i)+(x(i)-xmean)*(y(j)-ymean)*b
    sdx(i)=sdx(i)+b*(x(i)-xmean)**2
    sdy(i)=sdy(i)+b*(y(j)-ymean)**2
  enddo
enddo
!$acc end kernels

sdxg=0.0
sdyg=0.0
numg=0.0

do i=1,n
  sdxg=sdxg+sdx(i)
  sdyg=sdyg+sdy(i)
  numg=numg+num(i)
enddo

correlate_gaussian=numg/sqrt(sdxg*sdyg)

end function

end module
