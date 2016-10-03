module bootstrap

! JLR 25/8/2015
!
! Bootstrap confidence interval for Gaussian Kernel Correlation

use correlation

implicit none

contains

real function cpdf(x)
implicit none
real, intent(in) :: x

! cumulative probability density function
!
! Abramowitz and Stegun 26.2.17

real*8 :: z

if (x.ge.0) then
  z=1d0-0.5d0/(1d0+0.196854*x+0.115194*x**2+0.000344*x**3+0.019527*x**4)**4
else
  z=0.5d0/(1d0-0.196854*x+0.115194*x**2-0.000344*x**3+0.019527*x**4)**4
endif 

cpdf=z

end function

real function icpdf(x)
implicit none
real, intent(in) :: x

! inverse cumulative probability density function
!
! Algorithm AS 241 Wichura 1988


real*8 :: q,r,zp
real*8, parameter :: a0=3.3871327179d0, a1=50.434271938d0, a2=159.29113202d0
real*8, parameter :: a3=59.109374720d0, b1=17.895169469d0, b2=78.757757664d0
real*8, parameter :: b3=67.187563600d0, c0=1.4234372777d0, c1=2.7568153900d0
real*8, parameter :: c2=1.3067284816d0, c3=0.17023821103d0, d1=0.73700164250d0
real*8, parameter :: d2=0.12021132975d0, e0=6.6579051150d0, e1=3.0812263860d0
real*8, parameter :: e2=0.42868294337d0, e3=0.017337203997d0
real*8, parameter :: f1=0.24197894225d0, f2=0.012258202635d0

q=x-0.5d0

if (abs(q).le.0.425d0) then
  r=0.425d0**2-q**2
  zp=q*(((a3*r+a2)*r+a1)*r+a0)/(((b3*r+b2)*r+b1)*r+1d0)
else
  if ((x.eq.0d0).or.(x.eq.1d0)) then
    r=27d0
  else
    r=sqrt(-log(min(x,1d0-x)))
    r=min(r,27d0)
  endif
  if (r.le.5.0d0) then
    r=r-1.6d0
    zp=(((c3*r+c2)*r+c1)*r+c0)/((d2*r+d1)*r+1d0)
  else
    r=r-5.0d0
    zp=(((e3*r+e2)*r+e1)*r+e0)/((f2*r+f1)*r+1d0)
  endif
  if (q.lt.0) zp=-zp
endif

icpdf=zp

end function

subroutine add_one_sort(x)

implicit none

! add an element to a sorted list)

real, intent(inout) :: x(:)
real :: temp
integer :: i,j,n

n=size(x,1)

do i=1,n-1
  if (x(n)<x(i)) then
    temp=x(n)
    do j=n-1,i,-1
      x(j+1)=x(j)
    enddo
    x(i)=temp
    exit
  endif
enddo

end subroutine

subroutine joint_sort(a,ia,n)

implicit none

! sort a (union of two sorted lists, joined at n) perform sam operation on ia

real, intent(inout) :: a(:)
integer, intent(inout) :: ia(:)
integer, intent(in) :: n

integer :: i,j,k,m
real :: tmp(size(a,1))
integer :: itmp(size(a,1))

tmp(:)=a(:)
itmp(:)=ia(:)

i=1
j=n+1
k=1
m=size(a,1)

a(:)=-999
ia(:)=-999

do
  if (tmp(i).le.tmp(j)) then
    a(k)=tmp(i)
    ia(k)=itmp(i)
    i=i+1
  else
    a(k)=tmp(j)
    ia(k)=itmp(j)
    j=j+1
  endif
  k=k+1
  if (i>n) then
    a(k:m)=tmp(j:m)
    ia(k:m)=itmp(j:m)
    exit
  endif
  if (j>m) then
    a(k:m)=tmp(i:n)
    ia(k:m)=itmp(i:n)
    exit
  endif
!  if ((i.eq.n).and.(j.eq.m)) exit
enddo

end subroutine



function bootstrap_ci(x,y,tx,ty)
real ::  bootstrap_ci(1:2)

real, intent(in) :: x(:),y(:),tx(:),ty(:)
integer :: n,m,nm,i,j,k,iter,n_low,ierr,ci1,ci2,outer
real :: taux,tauy,txmean,tymean,rnd,p,logp,tmean,dt
real :: orig,pd,z0,a,astar,asum,a_num,a_den,b0,b1
real, allocatable :: temp(:),tstar(:)
real, allocatable :: nx(:),ny(:),ntx(:),nty(:)
real, allocatable :: bse(:)
integer, allocatable :: tmp_indx(:),indx(:,:),seq(:)
integer, parameter :: max_iter=2000
real :: med_min(25),med_max(25)
logical :: init(2)

n=size(x,1)
m=size(y,1)

allocate (temp(m+n),tstar(n+m),tmp_indx(n+m),indx(n+m,2),seq(n+m+1))
allocate (nx(n+m),ny(n+m),ntx(n+m),nty(n+m),bse(max_iter))

orig=correlate_gaussian(x,y,tx,ty)

! estimate persistence times

txmean=(tx(n)-tx(1))/(n-1)
tymean=(ty(m)-ty(1))/(m-1)
tmean=max(txmean,tymean)

do i=1,n/16
  temp(1:n)=tx(1:n)+i*txmean
  taux=correlate_gaussian(x,x,tx,temp(1:n))
  if (abs(taux)<0.368) exit
enddo
taux=i*txmean

do i=1,m/16
  temp(1:m)=ty(1:m)+i*tymean
  tauy=correlate_gaussian(y,y,temp(1:m),ty)
  if (abs(tauy)<0.368) exit
enddo
tauy=i*tymean

p=1.0-tmean/(4*max(taux,tauy))
logp=log(p)

! build integrated time-list
temp(1:n)=tx(1:n)
temp(n+1:n+m)=ty(1:m)
! index elements for series x  are positive integers, and negative for y 
do i=1,n
  tmp_indx(i)=i
enddo
do i=1,m
  tmp_indx(n+i)=-i
enddo
call joint_sort(temp,tmp_indx,n)
j=1
indx(:,:)=0 ! set no data value
do i=1,n+m
  tstar(j)=temp(i)
  if (tmp_indx(i)>0) then
    indx(j,1)=tmp_indx(i)
  else
    indx(j,2)=tmp_indx(i)
  endif
  if (i.lt.n+m) then
    if (temp(i).ne.temp(i+1)) j=j+1
  endif
enddo
nm=j

! do the bootstraps
do outer=1,25 ! outer iteration
 n_low=0
 iter=1
 do ! inner iteration
  j=1
  init(:)=.false.
  ! random starting point
  call random_number(rnd)
  i=floor((nm-1)*rnd)+1
  if (indx(i,1).ne.0) then
    seq(j)=indx(i,1)
    j=j+1
    init(1)=.true.
  endif
  if (indx(i,2).ne.0) then
    seq(j)=indx(i,2)
    j=j+1
    init(2)=.true.
  endif
  ! loop for rest of points
  do
    call random_number(rnd)
    if (i.gt.1) then
      dt=tstar(i)-tstar(i-1)
    else
      dt=tstar(2)-tstar(1)
    endif
    if (rnd.gt.exp(dt/tmean*logp)) then ! start a new block
      if (.not.(init(1).and.init(2))) then ! if haven't already sampled from both series, restart iteration
        j=1
        init(:)=.false.
      endif
      call random_number(rnd)
      i=floor((nm-1)*rnd)+1
      if (indx(i,1).ne.0) then
        seq(j)=indx(i,1)
        j=j+1
      endif
      if (indx(i,2).ne.0) then
        seq(j)=indx(i,2)
        j=j+1
      endif
    else ! continue with existing block
      i=i+1
      if (i.gt.nm) i=1
      if (indx(i,1).ne.0) then
        seq(j)=indx(i,1)
        j=j+1
        init(1)=.true.
      endif
      if (indx(i,2).ne.0) then
        seq(j)=indx(i,2)
        j=j+1
        init(2)=.true.
      endif
    endif
    if (j.gt.n+m) exit ! have enough data
  enddo
  ! built up two series and correspondind time bases
  ! find start times for each series
  init(:)=.false.
  j=1
  k=1
  do i=1,n+m
    if (seq(i)>0) then
      if (init(1)) then ! already have time base
        nx(j)=x(seq(i))
        if (seq(i).gt.1) then
          dt=tx(seq(i))-tx(seq(i)-1)
        else
          dt=tx(2)-tx(1)
        endif
        ntx(j)=ntx(j-1)+dt
      else
        nx(j)=x(seq(i))
        init(1)=.true.
        if (init(2)) then
          if (ty(-seq(i-1)).gt.tx(seq(i))) then
            if (seq(i).gt.1) then
              dt=tx(seq(i))-tx(seq(i)-1)
            else
              dt=tx(2)-tx(1)
            endif
            ntx(j)=nty(k-1)+dt
          else
            ntx(j)=nty(k-1)+tx(seq(i))-ty(-seq(i-1))
          endif
        else
          ntx(j)=tx(seq(i))
        endif
      endif
      j=j+1
    else
      if (init(2)) then ! already have time base
        ny(k)=y(-seq(i))
        if (-seq(i).gt.1) then
          dt=ty(-seq(i))-ty(-seq(i)-1)
        else
          dt=ty(2)-ty(1)
        endif
        nty(k)=nty(k-1)+dt
      else
        ny(k)=y(-seq(i))
        init(2)=.true.
        if (init(1)) then
          if (tx(seq(i-1)).gt.ty(-seq(i))) then
            if (-seq(i).gt.1) then
              dt=ty(-seq(i))-ty(-seq(i)-1)
            else
              dt=ty(2)-ty(1)
            endif
            nty(k)=ntx(j-1)+dt
          else
            nty(k)=ntx(j-1)+ty(-seq(i))-tx(seq(i-1))
          endif
        else
          nty(k)=ty(-seq(i))
        endif
      endif
      k=k+1
    endif
  enddo
  j=j-1
  k=k-1
  bse(iter)=correlate_gaussian(nx(1:j),ny(1:k),ntx(1:j),nty(1:k))
  if (bse(iter)<orig) n_low=n_low+1
  if (iter>1) call add_one_sort(bse(1:iter))
  iter=iter+1
  if (iter.gt.max_iter) exit
 enddo

!calcalate bias parameter
 z0=icpdf(real(n_low)/max_iter)

!calculate accelation parameter
 asum=0.0
 do i=1,max_iter
   asum=asum+bse(i)
 enddo
 astar=0.0
 do i=1,max_iter
   astar=astar+(asum-bse(i))/(max_iter-1)
 enddo
 astar=astar/max_iter
 a_num=0.0
 a_den=0.0
 do i=1,max_iter
   a_num=a_num+(astar-(asum-bse(i))/(max_iter-1))**3
   a_den=a_den+(astar-(asum-bse(i))/(max_iter-1))**2
 enddo
 a=a_num/(6.0*exp(1.5*log(a_den)))

! calculate interval points
 b0=z0+(z0-1.960)/(1-a*(z0-1.960))
 b1=z0+(z0+1.960)/(1-a*(z0+1.960))
 
 ci1=min(max(floor(iter*cpdf(b0)),1),max_iter)
 ci2=max(min(ceiling(iter*cpdf(b1)),max_iter),1)

 med_min(outer)=bse(ci1)
 med_max(outer)=bse(ci2)

 if (outer>1) then
   call add_one_sort(med_min(1:outer))
   call add_one_sort(med_max(1:outer))
 endif

! print *,bse(ci1),bse(ci2)

enddo

bootstrap_ci(1)=med_min(13)
bootstrap_ci(2)=med_max(13)

deallocate (temp,tstar,indx,tmp_indx,seq)
deallocate (nx,ny,ntx,nty,bse)


end function


end module
