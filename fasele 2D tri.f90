Parameter (n=500,m=400)
double precision f(0:8,0:n,0:m),d,fstar(0:8,0:n,0:m),fd(0:8,0:n,0:m)
double precision feq(0:8,0:n,0:m),rho(-1:n,-1:m)
double precision ubfx(1:1360),ubfy(1:1360),gama(1:1360),x_p1(1:1360),y_p1(1:1360),x_p2(1:1360),y_p2(1:1360),x_p3(1:1360),y_p3(1:1360),x_p4(1:1360),y_p4(1:1360)
double precision w(0:8),cx(0:8),cy(0:8),x_w1(1:1360),y_w1(1:1360),x_w2(1:1360),y_w2(1:1360),x_w3(1:1360),y_w3(1:1360),x_w4(1:1360),y_w4(1:1360),x_w5(1:1360),y_w5(1:1360)
double precision u(0:n,0:m),v(0:n,0:m),uo,delta_1(1:1360),delta_2(1:1360),delta_3(1:1360),delta_4(1:1360),delta_5(1:1360),delta(1:1360)
double precision omega_total(0:n,0:m),nu0(0:n,0:m),pi(0:n,0:m),Q(0:n,0:m),S_bar(0:n,0:m),nu_total(0:n,0:m),tau_total(0:n,0:m),haj(0:n,0:m),bar(0:n,0:m),turbo(0:n,0:m)
double precision umean(0:n,0:m),ubar(0:n,0:m),vmean(0:n,0:m),vbar(0:n,0:m),uv(0:n,0:m),uvmean(0:n,0:m),uvbar(0:n,0:m),uprimvprimbar(0:n,0:m),uu(0:n,0:m),uumean(0:n,0:m),uubar(0:n,0:m),uprim2bar(0:n,0:m),vv(0:n,0:m),vvmean(0:n,0:m),vvbar(0:n,0:m),vprim2bar(0:n,0:m),kinetic(0:n,0:m)
integer i,j,k
real h
open(2,file="uvfield.plt")
open(3,file="ubar(i,10D).plt")
open(4,file="vbar(i,10D).plt")
open(5,file="ubar(10D,j).plt")
open(6,file="vbar(10D,j).plt")
open(7,file="ubar(9.5D,j).plt")
open(8,file="vbar(9.5D,j).plt")
open(9,file="uvbar(i,10D).plt")
open(10,file="uvbar(10D,j).plt")
open(11,file="uvbar(9.5D,j).plt")
open(12,file="kinetic(i,10D).plt")
open(13,file="kinetic(10D,j).plt")
open(14,file="kinetic(9.5D,j).plt")
open(15,file="uubar(i,10D).plt")
open(16,file="uubar(10D,j).plt")
open(17,file="uubar(9.5D,j).plt")
open(18,file="vvbar(i,10D).plt")
open(19,file="vvbar(10D,j).plt")
open(20,file="vvbar(9.5D,j).plt")
cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
!uo=0.02
sumvelo=0.0
rhoo=1.00
dx=1.0
dy=dx
dt=1.0
!alpha=0.018733333
alpha=0.004
!Re=uo*m/alpha
Re=10000
r=10
rr=2.0*r
uo=(alpha*Re)/(m)
!print *, "Re=", Re
!pause
omega1=1.0/(3.*alpha+0.5)
tau1=1/omega1
cs=0.15
delta=1

mstep=12000

do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
omega_total(i,j)=omega1
tau_total(i,j)=tau1
nu0(i,j)=alpha
haj(i,j)=tau1
end do
end do

do j=0,m
u(0,j)=uo
v(0,j)=0.0
end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Geometry of circle #1
s=0
x0=7.5*rr
y0=11.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(1)=x0+r
y_w1(1)=y0
x_w1(21)=x0
y_w1(21)=y0+r
x_w1(41)=x0-r
y_w1(41)=y0
x_w1(61)=x0
y_w1(61)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=1,80
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=1,80
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=2,20
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=22,40
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=42,60
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=62,80
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Geometry of circle #2
s=80
x0=10.5*rr
y0=11.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(81)=x0+r
y_w1(81)=y0
x_w1(101)=x0
y_w1(101)=y0+r
x_w1(121)=x0-r
y_w1(121)=y0
x_w1(141)=x0
y_w1(141)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=81,160
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=81,160
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=82,100
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=102,120
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=122,140
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=142,160
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!!Geometry of circle #3
s=160
x0=13.5*rr
y0=11.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(161)=x0+r
y_w1(161)=y0
x_w1(181)=x0
y_w1(181)=y0+r
x_w1(201)=x0-r
y_w1(201)=y0
x_w1(221)=x0
y_w1(221)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=161,240
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=161,240
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=162,180
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=182,200
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=202,220
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=222,240
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #4
s=240
x0=16.5*rr
y0=11.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(241)=x0+r
y_w1(241)=y0
x_w1(261)=x0
y_w1(261)=y0+r
x_w1(281)=x0-r
y_w1(281)=y0
x_w1(301)=x0
y_w1(301)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=241,320
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=241,320
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=242,260
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=262,280
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=282,300
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=302,320
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #5
s=320
x0=9.0*rr
y0=10.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(321)=x0+r
y_w1(321)=y0
x_w1(341)=x0
y_w1(341)=y0+r
x_w1(361)=x0-r
y_w1(361)=y0
x_w1(381)=x0
y_w1(381)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=321,400
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=321,400
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=322,340
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=342,360
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=362,380
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=382,400
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #6
s=400
x0=12.0*rr
y0=10.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(401)=x0+r
y_w1(401)=y0
x_w1(421)=x0
y_w1(421)=y0+r
x_w1(441)=x0-r
y_w1(441)=y0
x_w1(461)=x0
y_w1(461)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=401,480
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=401,480
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=402,420
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=422,440
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=442,460
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=462,480
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #7
s=480
x0=15.0*rr
y0=10.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(481)=x0+r
y_w1(481)=y0
x_w1(501)=x0
y_w1(501)=y0+r
x_w1(521)=x0-r
y_w1(521)=y0
x_w1(541)=x0
y_w1(541)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=481,560
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=481,560
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=482,500
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=502,520
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=522,540
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=542,560
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #8
s=560
x0=7.5*rr
y0=8.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(561)=x0+r
y_w1(561)=y0
x_w1(581)=x0
y_w1(581)=y0+r
x_w1(601)=x0-r
y_w1(601)=y0
x_w1(621)=x0
y_w1(621)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=561,640
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=561,640
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=562,580
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=582,600
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=602,620
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=622,640
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #9
s=640
x0=10.5*rr
y0=8.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(641)=x0+r
y_w1(641)=y0
x_w1(661)=x0
y_w1(661)=y0+r
x_w1(681)=x0-r
y_w1(681)=y0
x_w1(701)=x0
y_w1(701)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=641,720
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=641,720
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=642,660
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=662,680
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=682,700
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=702,720
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #10
s=720
x0=13.5*rr
y0=8.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(721)=x0+r
y_w1(721)=y0
x_w1(741)=x0
y_w1(741)=y0+r
x_w1(761)=x0-r
y_w1(761)=y0
x_w1(781)=x0
y_w1(781)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=721,800
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=721,800
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=722,740
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=742,760
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=762,780
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=782,800
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #11
s=800
x0=16.5*rr
y0=8.5*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(801)=x0+r
y_w1(801)=y0
x_w1(821)=x0
y_w1(821)=y0+r
x_w1(841)=x0-r
y_w1(841)=y0
x_w1(861)=x0
y_w1(861)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=801,880
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=801,880
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=802,820
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=822,840
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=842,860
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=862,880
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #12
s=880
x0=9.0*rr
y0=13.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(881)=x0+r
y_w1(881)=y0
x_w1(901)=x0
y_w1(901)=y0+r
x_w1(921)=x0-r
y_w1(921)=y0
x_w1(941)=x0
y_w1(941)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=881,960
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=881,960
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=882,900
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=902,920
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=922,941
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=942,960
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #13
s=960
x0=12.0*rr
y0=13.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(961)=x0+r
y_w1(961)=y0
x_w1(981)=x0
y_w1(981)=y0+r
x_w1(1001)=x0-r
y_w1(1001)=y0
x_w1(1021)=x0
y_w1(1021)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=961,1040
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=961,1040
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=962,980
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=982,1001
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=1002,1020
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=1022,1040
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #14
s=1040
x0=15.0*rr
y0=13.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(1041)=x0+r
y_w1(1041)=y0
x_w1(1061)=x0
y_w1(1061)=y0+r
x_w1(1081)=x0-r
y_w1(1081)=y0
x_w1(1101)=x0
y_w1(1101)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=1041,1120
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=1041,1120
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=1042,1060
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=1062,1080
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=1082,1100
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=1102,1120
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #15
s=1120
x0=9.0*rr
y0=7.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(1121)=x0+r
y_w1(1121)=y0
x_w1(1141)=x0
y_w1(1141)=y0+r
x_w1(1161)=x0-r
y_w1(1161)=y0
x_w1(1181)=x0
y_w1(1181)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=1121,1200
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=1121,1200
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=1122,1140
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=1142,1160
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=1162,1180
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=1182,1200
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #16
s=1200
x0=12.0*rr
y0=7.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(1201)=x0+r
y_w1(1201)=y0
x_w1(1221)=x0
y_w1(1221)=y0+r
x_w1(1241)=x0-r
y_w1(1241)=y0
x_w1(1261)=x0
y_w1(1261)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=1201,1280
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=1201,1280
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=1202,1220
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=1222,1240
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=1242,1260
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=1262,1280
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
!Geometry of circle #17
s=1280
x0=15.0*rr
y0=7.0*rr
do h=0,6.283,0.07953164556962025316455696202532
s=s+1
x_w1(s)=x0+r*cos(h)
y_w1(s)=y0+r*sin(h)
x_w1(1281)=x0+r
y_w1(1281)=y0
x_w1(1301)=x0
y_w1(1301)=y0+r
x_w1(1321)=x0-r
y_w1(1321)=y0
x_w1(1341)=x0
y_w1(1341)=y0-r
!print *,x_w1(s),y_w1(s),s
end do
!pause
!stop
!********************************************************************************
!„‘Œ’ ò—œ‰ ‰ﬁ«ÿ p1,p2,p3,p4 »—«? ò· œ«?—Â

do k=1281,1360
do i=x0-r,x0+r
d=i-x_w1(k)
if (d<=1) then
    x_p1(k)=i
    x_p2(k)=i
    x_p3(k)=i-1
    x_p4(k)=i-1
else
    continue
end if
end do
!print *,x_p3(k),x_p4(k),x_w1(k),x_p1(k),x_p2(k)
end do
!pause
!stop
    
do k=1281,1360
do j=y0-r,y0+r
d=j-y_w1(k)
if (d<=1) then
    y_p1(k)=j-1
    y_p2(k)=j
    y_p3(k)=j
    y_p4(k)=j-1
else
    continue
end if
end do
!print *,y_p1(k),y_p4(k),y_w1(k),y_p3(k),y_p2(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ ‰ﬁ«ÿ Ê œ· « »—«? —»⁄ «Ê·

do k=1282,1300
a=y_p4(k)-x_p4(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
if (y_p2(k)>=y_w1(k)>=y_p4(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p4(k)+x_w1(k)+y_p4(k)
end if
o=((x_p2(k)-x_w1(k))**2+(y_p2(k)-y_w1(k))**2)**0.5
c=((x_p2(k)-x_p4(k))**2+(y_p2(k)-y_p4(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)

x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*******************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ œÊ„

do k=1302,1320
a=y_p1(k)+x_p1(k)-y0
b=-1*(a+x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
if (y_p3(k)<=y_w1(k)<=y_p1(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p1(k)-x_w1(k)+y_p1(k)
end if
o=((x_p3(k)-x_w1(k))**2+(y_p3(k)-y_w1(k))**2)**0.5
c=((x_p3(k)-x_p1(k))**2+(y_p3(k)-y_p1(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w5(k)=x_w5(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p2(k)-x_w2(k))**2+(y_p2(k)-y_w2(k))**2)**0.5
c=((x_p2(k)-x_p1(k))**2+(y_p2(k)-y_p1(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
x_w4(k)=x_w4(k)
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p3(k)-x_w4(k))**2+(y_p3(k)-y_w4(k))**2)**0.5
c=((x_p3(k)-x_p4(k))**2+(y_p3(k)-y_p4(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?Ê Ê œ· «Â« »—«? —»⁄ ”Ê„

do k=1322,1340
a=y_p2(k)-x_p2(k)-y0
b=(a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
if (y_p4(k)<=y_w1(k)<=y_p2(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=-1*x_p2(k)+x_w1(k)+y_p2(k)
end if
o=((x_p4(k)-x_w1(k))**2+(y_p4(k)-y_w1(k))**2)**0.5
c=((x_p4(k)-x_p2(k))**2+(y_p4(k)-y_p2(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p4(k),y_p4(k),x_w1(k),y_w1(k),x_p2(k),y_p2(k),delta_1(k)

x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w1(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p4(k)-x_w5(k))**2+(y_p4(k)-y_w5(k))**2)**0.5
c=((x_p4(k)-x_p1(k))**2+(y_p4(k)-y_p1(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p3(k)-x_w3(k))**2+(y_p3(k)-y_w3(k))**2)**0.5
c=((x_p3(k)-x_p2(k))**2+(y_p3(k)-y_p2(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!pause
!stop
!*************************************************************
!»œ”  ¬Ê—œ‰ œ»·?ÊÂ« Ê œ· «Â« »—«? —»⁄ çÂ«—„

do k=1342,1360
a=y_p3(k)+x_p3(k)-y0
b=(-1*a-x0)
c=(a**2-r**2+x0**2)/2
del=(b**2)-4*c
x_w1(k)=((-1*b)-(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
if (y_p1(k)<=y_w1(k).and.y_w1(k)<=y_p3(k)) then
x_w1(k)=x_w1(k)
else
x_w1(k)=((-1*b)+(del**0.5))/2
y_w1(k)=x_p3(k)-x_w1(k)+y_p3(k)
end if
o=((x_p1(k)-x_w1(k))**2+(y_p1(k)-y_w1(k))**2)**0.5
c=((x_p1(k)-x_p3(k))**2+(y_p1(k)-y_p3(k))**2)**0.5
delta_1(k)=abs(o/c)
!print *,x_p1(k),y_p1(k),x_w1(k),y_w1(k),x_p3(k),y_p3(k),delta_1(k)


x_w5(k)=x0+((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
x_w5(k)=x_w5(k)
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
x_w5(k)=x0-((r**2)-(y_p1(k)-y0)**2)**0.5
y_w5(k)=y_p1(k)
if (x_p4(k)<x_w5(k).and.x_w5(k)<x_p1(k)) then
o=((x_p1(k)-x_w5(k))**2+(y_p1(k)-y_w5(k))**2)**0.5
c=((x_p1(k)-x_p4(k))**2+(y_p1(k)-y_p4(k))**2)**0.5
delta_5(k)=abs(o/c)
else
delta_5(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w5(k),y_w5(k),x_p1(k),y_p1(k),delta_5(k)


x_w3(k)=x0+((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
x_w3(k)=x_w3(k)
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
x_w3(k)=x0-((r**2)-(y_p2(k)-y0)**2)**0.5
y_w3(k)=y_p2(k)
if (x_p3(k)<x_w3(k).and.x_w3(k)<x_p2(k)) then
o=((x_p2(k)-x_w3(k))**2+(y_p2(k)-y_w3(k))**2)**0.5
c=((x_p2(k)-x_p3(k))**2+(y_p2(k)-y_p3(k))**2)**0.5
delta_3(k)=abs(o/c)
else
delta_3(k)=0
end if
end if
!print *,x_p3(k),y_p3(k),x_w3(k),y_w3(k),x_p2(k),y_p2(k),delta_3(k)


x_w2(k)=x_p1(k)
y_w2(k)=y0-((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
y_w2(k)=y_w2(k)
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
x_w2(k)=x_p1(k)
y_w2(k)=y0+((r**2)-(x_p1(k)-x0)**2)**0.5
if (y_p1(k)<y_w2(k).and.y_w2(k)<y_p2(k)) then
o=((x_p1(k)-x_w2(k))**2+(y_p1(k)-y_w2(k))**2)**0.5
c=((x_p1(k)-x_p2(k))**2+(y_p1(k)-y_p2(k))**2)**0.5
delta_2(k)=abs(o/c)
else
delta_2(k)=0
end if
end if
!print *,x_p1(k),y_p1(k),x_w2(k),y_w2(k),x_p2(k),y_p2(k),delta_2(k)

x_w4(k)=x_p4(k)
y_w4(k)=y0-((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
y_w4(k)=y_w4(k)
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
x_w4(k)=x_p4(k)
y_w4(k)=y0+((r**2)-(x_p4(k)-x0)**2)**0.5
if (y_p4(k)<y_w4(k).and.y_w4(k)<y_p3(k)) then
o=((x_p4(k)-x_w4(k))**2+(y_p4(k)-y_w4(k))**2)**0.5
c=((x_p4(k)-x_p3(k))**2+(y_p4(k)-y_p3(k))**2)**0.5
delta_4(k)=abs(o/c)
else
delta_4(k)=0
end if
end if
!print *,x_p4(k),y_p4(k),x_w4(k),y_w4(k),x_p3(k),y_p3(k),delta_4(k)
end do
!=================================================================
u_old=0
do kk=1,mstep
!do while (error<1e-7)
    u_old=v(n/2,m/2)
do i=0,n
do j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
pi(i,j)=0
!rhou(i,j)=0
do k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
pi(i,j)=pi(i,j)+cx(k)*cy(k)*(f(k,i,j)-feq(k,i,j))
!rhou(i,j)=rhou(i,j)+f(k,i,j)
end do
bar(i,j)=(2.0*pi(i,j)*pi(i,j))**0.5
!nonequilibrium mode
turbo(i,j)=0.5*((haj(i,j)**2+(18.0*bar(i,j)*(cs)**2))**0.5-haj(i,j))
!finite differance mode
!tau_turb(i,j)=((cs*delta)**2*Q_bar(i,j))/(2*rhou(i,j)*tau_total(i,j)*dt**2)
tau_total(i,j)=haj(i,j)+turbo(i,j)
omega_total(i,j)=1/tau_total(i,j)
do k=0,8
f(k,i,j)=omega_total(i,j)*feq(k,i,j)+(1.-omega_total(i,j))*f(k,i,j)
end do

end do
end do

do j=0,m
do i=n,1,-1 
f(1,i,j)=f(1,i-1,j)
end do
do i=0,n-1  
f(3,i,j)=f(3,i+1,j)
end do
end do
do j=m,1,-1
do i=0,n
f(2,i,j)=f(2,i,j-1)
end do
do i=n,1,-1
f(5,i,j)=f(5,i-1,j-1)
end do
do i=0,n-1
f(6,i,j)=f(6,i+1,j-1)
end do
end do
do j=0,m-1   
do i=0,n
f(4,i,j)=f(4,i,j+1)
end do
do i=0,n-1
f(7,i,j)=f(7,i+1,j+1)
end do
do i=n,1,-1
f(8,i,j)=f(8,i-1,j+1)
end do
end do


 do j=0,m
! bounce back on west boundary
rhow=(f(0,0,j)+f(2,0,j)+f(4,0,j)+2.*(f(3,0,j)+f(6,0,j)+f(7,0,j)))/(1.-uo)
f(1,0,j)=f(3,0,j)+2.*rhow*uo/3
f(5,0,j)=f(7,0,j)+rhow*uo/6-0.5*(f(4,0,j)-f(2,0,j))
f(8,0,j)=f(6,0,j)+rhow*uo/6+0.5*(f(4,0,j)-f(2,0,j))

 end do
! bounce back on south boundary
do i=0,n
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
end do
! bounce back on north boundary
do i=0,n
f(4,i,m)=f(2,i,m)
f(8,i,m)=f(6,i,m)
f(7,i,m)=f(5,i,m)
end do
do j=1,m
f(1,n,j)=2.*f(1,n-1,j)-f(1,n-2,j)
f(5,n,j)=2.*f(5,n-1,j)-f(5,n-2,j)
f(8,n,j)=2.*f(8,n-1,j)-f(8,n-2,j)
end do

!******************************************************************
!Calculating f f for circle #1
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=7.5*rr
y0=11.5*rr
do k=2,20
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=22,40
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=42,60
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=62,80
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(1)=1
x_w5(1)=x0+r
y_w5(1)=y0
x_p1(1)=x0+r+1
y_p1(1)=y0
x_p4(1)=x0+r
y_p4(1)=y0
k=1
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(21)=1
x_w2(21)=x0
y_w2(21)=y0+r
x_p1(21)=x0
y_p1(21)=y0+r
x_p2(21)=x0
y_p2(21)=y0+r+1
k=21
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(41)=1
x_w5(41)=x0-r
y_w5(41)=y0
x_p1(41)=x0-r
y_p1(41)=y0
x_p4(41)=x0-r-1
y_p4(41)=y0
k=41
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(61)=1
x_w2(61)=x0
y_w2(61)=y0-r
x_p1(61)=x0
y_p1(61)=y0-r-1
x_p2(61)=x0
y_p2(61)=y0-r
k=61
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f for circle #2
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=10.5*rr
y0=11.5*rr
do k=82,100
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=102,120
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=122,140
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=142,160
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(81)=1
x_w5(81)=x0+r
y_w5(81)=y0
x_p1(81)=x0+r+1
y_p1(81)=y0
x_p4(81)=x0+r
y_p4(81)=y0
k=81
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(101)=1
x_w2(101)=x0
y_w2(101)=y0+r
x_p1(101)=x0
y_p1(101)=y0+r
x_p2(101)=x0
y_p2(101)=y0+r+1
k=101
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_5(121)=1
x_w5(121)=x0-r
y_w5(121)=y0
x_p1(121)=x0-r
y_p1(121)=y0
x_p4(121)=x0-r-1
y_p4(121)=y0
k=121
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(141)=1
x_w2(141)=x0
y_w2(141)=y0-r
x_p1(141)=x0
y_p1(141)=y0-r-1
x_p2(141)=x0
y_p2(141)=y0-r
k=141
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
!******************************************************************
!Calculating f f for circle #3
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=13.5*rr
y0=11.5*rr
do k=162,180
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=182,200
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=202,220
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=222,240
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(161)=1
x_w5(161)=x0+r
y_w5(161)=y0
x_p1(161)=x0+r+1
y_p1(161)=y0
x_p4(161)=x0+r
y_p4(161)=y0
k=161
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(181)=1
x_w2(181)=x0
y_w2(181)=y0+r
x_p1(181)=x0
y_p1(181)=y0+r
x_p2(181)=x0
y_p2(181)=y0+r+1
k=181
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(201)=1
x_w5(201)=x0-r
y_w5(201)=y0
x_p1(201)=x0-r
y_p1(201)=y0
x_p4(201)=x0-r-1
y_p4(201)=y0
k=201
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(221)=1
x_w2(221)=x0
y_w2(221)=y0-r
x_p1(221)=x0
y_p1(221)=y0-r-1
x_p2(221)=x0
y_p2(221)=y0-r
k=221
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!******************************************************************
!Calculating f f for circle #4
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=16.5*rr
y0=11.5*rr
do k=242,260
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=262,280
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=282,300
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=302,320
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(241)=1
x_w5(241)=x0+r
y_w5(241)=y0
x_p1(241)=x0+r+1
y_p1(241)=y0
x_p4(241)=x0+r
y_p4(241)=y0
k=241
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(261)=1
x_w2(261)=x0
y_w2(261)=y0+r
x_p1(261)=x0
y_p1(261)=y0+r
x_p2(261)=x0
y_p2(261)=y0+r+1
k=261
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(281)=1
x_w5(281)=x0-r
y_w5(281)=y0
x_p1(281)=x0-r
y_p1(281)=y0
x_p4(281)=x0-r-1
y_p4(281)=y0
k=281
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(301)=1
x_w2(301)=x0
y_w2(301)=y0-r
x_p1(301)=x0
y_p1(301)=y0-r-1
x_p2(301)=x0
y_p2(301)=y0-r
k=301
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #5
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=9.0*rr
y0=10.0*rr
do k=322,340
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=342,360
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=362,380
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=382,400
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(321)=1
x_w5(321)=x0+r
y_w5(321)=y0
x_p1(321)=x0+r+1
y_p1(321)=y0
x_p4(321)=x0+r
y_p4(321)=y0
k=321
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(341)=1
x_w2(341)=x0
y_w2(341)=y0+r
x_p1(341)=x0
y_p1(341)=y0+r
x_p2(341)=x0
y_p2(341)=y0+r+1
k=341
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(361)=1
x_w5(361)=x0-r
y_w5(361)=y0
x_p1(361)=x0-r
y_p1(361)=y0
x_p4(361)=x0-r-1
y_p4(361)=y0
k=361
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(381)=1
x_w2(381)=x0
y_w2(381)=y0-r
x_p1(381)=x0
y_p1(381)=y0-r-1
x_p2(381)=x0
y_p2(381)=y0-r
k=381
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #6
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=12.0*rr
y0=10.0*rr
do k=402,420
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=422,440
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=442,460
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=462,480
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(401)=1
x_w5(401)=x0+r
y_w5(401)=y0
x_p1(401)=x0+r+1
y_p1(401)=y0
x_p4(401)=x0+r
y_p4(401)=y0
k=401
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(421)=1
x_w2(421)=x0
y_w2(421)=y0+r
x_p1(421)=x0
y_p1(421)=y0+r
x_p2(421)=x0
y_p2(421)=y0+r+1
k=421
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(441)=1
x_w5(441)=x0-r
y_w5(441)=y0
x_p1(441)=x0-r
y_p1(441)=y0
x_p4(441)=x0-r-1
y_p4(441)=y0
k=441
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(461)=1
x_w2(461)=x0
y_w2(461)=y0-r
x_p1(461)=x0
y_p1(461)=y0-r-1
x_p2(461)=x0
y_p2(461)=y0-r
k=461
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #7
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=15.0*rr
y0=10.0*rr
do k=482,500
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=502,520
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=522,540
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=542,560
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(481)=1
x_w5(481)=x0+r
y_w5(481)=y0
x_p1(481)=x0+r+1
y_p1(481)=y0
x_p4(481)=x0+r
y_p4(481)=y0
k=481
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(501)=1
x_w2(50)=x0
y_w2(501)=y0+r
x_p1(501)=x0
y_p1(501)=y0+r
x_p2(501)=x0
y_p2(501)=y0+r+1
k=501
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(521)=1
x_w5(521)=x0-r
y_w5(521)=y0
x_p1(521)=x0-r
y_p1(521)=y0
x_p4(521)=x0-r-1
y_p4(521)=y0
k=521
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(541)=1
x_w2(541)=x0
y_w2(541)=y0-r
x_p1(541)=x0
y_p1(541)=y0-r-1
x_p2(541)=x0
y_p2(541)=y0-r
k=541
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #8
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=7.5*rr
y0=8.5*rr
do k=562,580
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=582,600
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=602,620
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=622,640
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(561)=1
x_w5(561)=x0+r
y_w5(561)=y0
x_p1(561)=x0+r+1
y_p1(561)=y0
x_p4(561)=x0+r
y_p4(561)=y0
k=561
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(581)=1
x_w2(581)=x0
y_w2(581)=y0+r
x_p1(581)=x0
y_p1(581)=y0+r
x_p2(581)=x0
y_p2(581)=y0+r+1
k=581
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(601)=1
x_w5(601)=x0-r
y_w5(601)=y0
x_p1(601)=x0-r
y_p1(601)=y0
x_p4(601)=x0-r-1
y_p4(601)=y0
k=601
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(621)=1
x_w2(621)=x0
y_w2(621)=y0-r
x_p1(621)=x0
y_p1(621)=y0-r-1
x_p2(621)=x0
y_p2(621)=y0-r
k=621
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #9
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=10.5*rr
y0=8.5*rr
do k=642,660
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=662,680
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=682,700
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=702,720
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(641)=1
x_w5(641)=x0+r
y_w5(641)=y0
x_p1(641)=x0+r+1
y_p1(641)=y0
x_p4(641)=x0+r
y_p4(641)=y0
k=641
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(661)=1
x_w2(661)=x0
y_w2(661)=y0+r
x_p1(661)=x0
y_p1(661)=y0+r
x_p2(661)=x0
y_p2(661)=y0+r+1
k=661
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(681)=1
x_w5(681)=x0-r
y_w5(681)=y0
x_p1(681)=x0-r
y_p1(681)=y0
x_p4(681)=x0-r-1
y_p4(681)=y0
k=681
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(701)=1
x_w2(701)=x0
y_w2(701)=y0-r
x_p1(701)=x0
y_p1(701)=y0-r-1
x_p2(701)=x0
y_p2(701)=y0-r
k=701
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #10
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=13.5*rr
y0=8.5*rr
do k=722,740
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=742,760
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=762,780
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=782,800
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(721)=1
x_w5(721)=x0+r
y_w5(721)=y0
x_p1(721)=x0+r+1
y_p1(721)=y0
x_p4(721)=x0+r
y_p4(721)=y0
k=721
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(741)=1
x_w2(741)=x0
y_w2(741)=y0+r
x_p1(741)=x0
y_p1(741)=y0+r
x_p2(741)=x0
y_p2(741)=y0+r+1
k=741
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(761)=1
x_w5(761)=x0-r
y_w5(761)=y0
x_p1(761)=x0-r
y_p1(761)=y0
x_p4(761)=x0-r-1
y_p4(761)=y0
k=761
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(781)=1
x_w2(781)=x0
y_w2(781)=y0-r
x_p1(781)=x0
y_p1(781)=y0-r-1
x_p2(781)=x0
y_p2(781)=y0-r
k=781
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #11
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=16.5*rr
y0=8.5*rr
do k=802,820
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=822,840
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=842,860
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=862,880
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(801)=1
x_w5(801)=x0+r
y_w5(801)=y0
x_p1(801)=x0+r+1
y_p1(801)=y0
x_p4(801)=x0+r
y_p4(801)=y0
k=801
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(821)=1
x_w2(821)=x0
y_w2(821)=y0+r
x_p1(821)=x0
y_p1(821)=y0+r
x_p2(821)=x0
y_p2(821)=y0+r+1
k=821
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(841)=1
x_w5(841)=x0-r
y_w5(841)=y0
x_p1(841)=x0-r
y_p1(841)=y0
x_p4(841)=x0-r-1
y_p4(841)=y0
k=841
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(861)=1
x_w2(861)=x0
y_w2(861)=y0-r
x_p1(861)=x0
y_p1(861)=y0-r-1
x_p2(861)=x0
y_p2(861)=y0-r
k=861
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #12
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=9.0*rr
y0=13.0*rr
do k=882,900
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=902,920
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=922,940
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=942,960
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(881)=1
x_w5(881)=x0+r
y_w5(881)=y0
x_p1(881)=x0+r+1
y_p1(881)=y0
x_p4(881)=x0+r
y_p4(881)=y0
k=881
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(901)=1
x_w2(901)=x0
y_w2(901)=y0+r
x_p1(901)=x0
y_p1(901)=y0+r
x_p2(901)=x0
y_p2(901)=y0+r+1
k=901
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(921)=1
x_w5(921)=x0-r
y_w5(921)=y0
x_p1(921)=x0-r
y_p1(921)=y0
x_p4(921)=x0-r-1
y_p4(921)=y0
k=921
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(941)=1
x_w2(941)=x0
y_w2(941)=y0-r
x_p1(941)=x0
y_p1(941)=y0-r-1
x_p2(941)=x0
y_p2(941)=y0-r
k=941
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #13
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=12.0*rr
y0=13.0*rr
do k=962,980
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=982,1000
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=1002,1020
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=1022,1040
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(961)=1
x_w5(961)=x0+r
y_w5(961)=y0
x_p1(961)=x0+r+1
y_p1(961)=y0
x_p4(961)=x0+r
y_p4(961)=y0
k=961
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(981)=1
x_w2(981)=x0
y_w2(981)=y0+r
x_p1(981)=x0
y_p1(981)=y0+r
x_p2(981)=x0
y_p2(981)=y0+r+1
k=981
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(1001)=1
x_w5(1001)=x0-r
y_w5(1001)=y0
x_p1(1001)=x0-r
y_p1(1001)=y0
x_p4(1001)=x0-r-1
y_p4(1001)=y0
k=1001
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(1021)=1
x_w2(1021)=x0
y_w2(1021)=y0-r
x_p1(1021)=x0
y_p1(1021)=y0-r-1
x_p2(1021)=x0
y_p2(1021)=y0-r
k=1021
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #14
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=15.0*rr
y0=13.0*rr
do k=1042,1060
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=1062,1080
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=1082,1100
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=1102,1120
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(1041)=1
x_w5(1041)=x0+r
y_w5(1041)=y0
x_p1(1041)=x0+r+1
y_p1(1041)=y0
x_p4(1041)=x0+r
y_p4(1041)=y0
k=1041
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(1061)=1
x_w2(1061)=x0
y_w2(1061)=y0+r
x_p1(1061)=x0
y_p1(1061)=y0+r
x_p2(1061)=x0
y_p2(1061)=y0+r+1
k=1061
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(1081)=1
x_w5(1081)=x0-r
y_w5(1081)=y0
x_p1(1081)=x0-r
y_p1(1081)=y0
x_p4(1081)=x0-r-1
y_p4(1081)=y0
k=1081
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(1101)=1
x_w2(1101)=x0
y_w2(1101)=y0-r
x_p1(1101)=x0
y_p1(1101)=y0-r-1
x_p2(1101)=x0
y_p2(1101)=y0-r
k=1101
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #15
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=9.0*rr
y0=7.0*rr
do k=1122,1140
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=1142,1160
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=1162,1180
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=1182,1200
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(1121)=1
x_w5(1121)=x0+r
y_w5(1121)=y0
x_p1(1121)=x0+r+1
y_p1(1121)=y0
x_p4(1121)=x0+r
y_p4(1121)=y0
k=1121
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(1141)=1
x_w2(1141)=x0
y_w2(1141)=y0+r
x_p1(1141)=x0
y_p1(1141)=y0+r
x_p2(1141)=x0
y_p2(1141)=y0+r+1
k=1141
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(1161)=1
x_w5(1161)=x0-r
y_w5(1161)=y0
x_p1(1161)=x0-r
y_p1(1161)=y0
x_p4(1161)=x0-r-1
y_p4(1161)=y0
k=1161
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(1181)=1
x_w2(1181)=x0
y_w2(1181)=y0-r
x_p1(1181)=x0
y_p1(1181)=y0-r-1
x_p2(1181)=x0
y_p2(1181)=y0-r
k=1181
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #16
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=12.0*rr
y0=7.0*rr
do k=1202,1220
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=1222,1240
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=1242,1260
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=1262,1280
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(1201)=1
x_w5(1201)=x0+r
y_w5(1201)=y0
x_p1(1201)=x0+r+1
y_p1(1201)=y0
x_p4(1201)=x0+r
y_p4(1201)=y0
k=1201
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(1221)=1
x_w2(1221)=x0
y_w2(1221)=y0+r
x_p1(1221)=x0
y_p1(1221)=y0+r
x_p2(1221)=x0
y_p2(1221)=y0+r+1
k=1221
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(1241)=1
x_w5(1241)=x0-r
y_w5(1241)=y0
x_p1(1241)=x0-r
y_p1(1241)=y0
x_p4(1241)=x0-r-1
y_p4(1241)=y0
k=1241
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(1261)=1
x_w2(1261)=x0
y_w2(1261)=y0-r
x_p1(1261)=x0
y_p1(1261)=y0-r-1
x_p2(1261)=x0
y_p2(1261)=y0-r
k=1261
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************
!Calculating f f for circle #17
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ «Ê·
x0=15.0*rr
y0=7.0*rr
do k=1282,1300
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else
ubfx(k)=u(x_p2(k)+1,y_p2(k)+1)
ubfy(k)=v(x_p2(k)+1,y_p2(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p2(k)+1,y_p2(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if


if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
if (delta_5(k)>0) then
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if


if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

end do
!******************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ œÊ„
do k=1302,1320
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else
ubfx(k)=u(x_p3(k)-1,y_p3(k)+1)
ubfy(k)=v(x_p3(k)-1,y_p3(k)+1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p3(k)-1,y_p3(k)+1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p2(k),y_p2(k)+1)
ubfy(k)=v(x_p2(k),y_p2(k)+1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k)+1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p3(k),y_p3(k)+1)
ubfy(k)=v(x_p3(k),y_p3(k)+1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p3(k),y_p3(k)+1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ ”Ê„
do k=1322,1340
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else
ubfx(k)=u(x_p4(k)-1,y_p4(k)-1)
ubfy(k)=v(x_p4(k)-1,y_p4(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p4(k)-1,y_p4(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p4(k)-1,y_p4(k))
ubfy(k)=v(x_p4(k)-1,y_p4(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k)-1,y_p4(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p3(k),y_p3(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p3(k),y_p3(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k),y_p3(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p3(k)-1,y_p3(k))
ubfy(k)=v(x_p3(k)-1,y_p3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p3(k)-1,y_p3(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p3(k),y_p3(k))**2+v(x_p3(k),y_p3(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p3(k),y_p3(k))*cx(s)+v(x_p3(k),y_p3(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p3(k),y_p3(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p3(k),y_p3(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!»œ”  ¬Ê—œ‰ «› »—«? —»⁄ çÂ«—„
do k=1342,1360
if (delta_1(k)>=0.5) then
ubfx(k)=(2*delta_1(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_1(k))
ubfy(k)=(2*delta_1(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_1(k))
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else
ubfx(k)=u(x_p1(k)+1,y_p1(k)-1)
ubfy(k)=v(x_p1(k)+1,y_p1(k)-1)
gama(k)=(2*delta_1(k)-1)/(tau_total(x_p1(k)+1,y_p1(k)-1)-2)
end if
if (delta_1(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_2(k)>=0.5) then
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_2(k).and.delta_2(k)<0.5) then
ubfx(k)=u(x_p1(k),y_p1(k)-1)
ubfy(k)=v(x_p1(k),y_p1(k)-1)
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k)-1)-2)
else 
continue
end if
if (delta_2(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do
else 
    continue
end if

if (delta_4(k)>=0.5) then
ubfx(k)=(2*delta_4(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_4(k))
ubfy(k)=(2*delta_4(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_4(k))
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
else if (0<delta_4(k).and.delta_4(k)<0.5) then
ubfx(k)=u(x_p4(k),y_p4(k)-1)
ubfy(k)=v(x_p4(k),y_p4(k)-1)
gama(k)=(2*delta_4(k)-1)/(tau_total(x_p4(k),y_p4(k)-1)-2)
else 
continue
end if
if (delta_4(k)>0) then
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if

if (delta_5(k)>=0.5) then
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
else if (0<delta_5(k).and.delta_5(k)<0.5) then
ubfx(k)=u(x_p1(k)+1,y_p1(k))
ubfy(k)=v(x_p1(k)+1,y_p1(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k)+1,y_p1(k))-2)
else 
continue
end if
if (delta_5(k)>0) then
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do
else 
    continue
end if

if (delta_3(k)>=0.5) then
ubfx(k)=(2*delta_3(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_3(k))
ubfy(k)=(2*delta_3(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_3(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
else if (0<delta_3(k).and.delta_3(k)<0.5) then
ubfx(k)=u(x_p2(k)+1,y_p2(k))
ubfy(k)=v(x_p2(k)+1,y_p2(k))
gama(k)=(2*delta_3(k)-1)/(tau_total(x_p2(k)+1,y_p2(k))-2)
else 
continue
end if
if (delta_3(k)>0) then
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p3(k),y_p3(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p3(k),y_p3(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p3(k),y_p3(k))
end do
else 
    continue
end if
end do
!**************************************************************************************
!Ê«—œ ò—œ‰ œ” ? ‰ﬁ«ÿ „—“? Ê „Õ«”»Â «› »—«? ¬‰Â«
delta_5(1281)=1
x_w5(1281)=x0+r
y_w5(1281)=y0
x_p1(1281)=x0+r+1
y_p1(1281)=y0
x_p4(1281)=x0+r
y_p4(1281)=y0
k=1281
ubfx(k)=(2*delta_5(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p4(k),y_p4(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p4(k),y_p4(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p4(k),y_p4(k))
end do

delta_2(1301)=1
x_w2(1301)=x0
y_w2(1301)=y0+r
x_p1(1301)=x0
y_p1(1301)=y0+r
x_p2(1301)=x0
y_p2(1301)=y0+r+1
k=1301
ubfx(k)=(2*delta_2(k)-3)*u(x_p2(k),y_p2(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p2(k),y_p2(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p2(k),y_p2(k))+0.5)
m1=u(x_p2(k),y_p2(k))**2+v(x_p2(k),y_p2(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p2(k),y_p2(k))*cx(s)+v(x_p2(k),y_p2(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p2(k),y_p2(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p2(k),y_p2(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do



delta_5(1321)=1
x_w5(1321)=x0-r
y_w5(1321)=y0
x_p1(1321)=x0-r
y_p1(1321)=y0
x_p4(1321)=x0-r-1
y_p4(1321)=y0
k=1321
ubfx(k)=(2*delta_5(k)-3)*u(x_p4(k),y_p4(k))/(2*delta_5(k))
ubfy(k)=(2*delta_5(k)-3)*v(x_p4(k),y_p4(k))/(2*delta_5(k))
gama(k)=(2*delta_5(k)-1)/(tau_total(x_p4(k),y_p4(k))+0.5)
m1=u(x_p4(k),y_p4(k))**2+v(x_p4(k),y_p4(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p4(k),y_p4(k))*cx(s)+v(x_p4(k),y_p4(k))*cy(s)
fstar(s,x_p1(k),y_p1(k))=w(s)*rho(x_p4(k),y_p4(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p1(k),y_p1(k))=(1-gama(k))*f(s,x_p4(k),y_p4(k))+gama(k)*fstar(s,x_p1(k),y_p1(k))
end do

delta_2(1341)=1
x_w2(1341)=x0
y_w2(1341)=y0-r
x_p1(1341)=x0
y_p1(1341)=y0-r-1
x_p2(1341)=x0
y_p2(1341)=y0-r
k=1341
ubfx(k)=(2*delta_2(k)-3)*u(x_p1(k),y_p1(k))/(2*delta_2(k))
ubfy(k)=(2*delta_2(k)-3)*v(x_p1(k),y_p1(k))/(2*delta_2(k))
gama(k)=(2*delta_2(k)-1)/(tau_total(x_p1(k),y_p1(k))+0.5)
m1=u(x_p1(k),y_p1(k))**2+v(x_p1(k),y_p1(k))**2
do s=0,8
m2=ubfx(k)*cx(s)+ubfy(k)*cy(s)
m3=u(x_p1(k),y_p1(k))*cx(s)+v(x_p1(k),y_p1(k))*cy(s)
fstar(s,x_p2(k),y_p2(k))=w(s)*rho(x_p1(k),y_p1(k))*(1+3*m2+4.5*m3*m3-1.5*m1)
f(s,x_p2(k),y_p2(k))=(1-gama(k))*f(s,x_p1(k),y_p1(k))+gama(k)*fstar(s,x_p2(k),y_p2(k))
end do

!******************************************************************


do j=0,m
do i=0,n
ssum=0.0
do k=0,8
ssum=ssum+f(k,i,j)
!print *,f(k,i,j)
end do
rho(i,j)=ssum
!print *,rho(i,j)
end do
end do
do i=0,n
rho(i,m)=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
end do
do i=1,n
do j=1,m-1
usum=0.0
vsum=0.0
do k=0,8
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
end do
u(i,j)=usum/rho(i,j)
v(i,j)=vsum/rho(i,j)
end do
end do
!if (kk>10000) then
!error=abs(v(n/2,m/2)-u_old)
!end if
!if (error<0.00000001.and.error>0.0) then
!exit
!end if

! Circle #1 
x0=7.5*rr
y0=11.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #2
x0=10.5*rr
y0=11.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #3 
x0=13.5*rr
y0=11.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #4 
x0=16.5*rr
y0=11.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #5
x0=9.0*rr
y0=10.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #6
x0=12.0*rr
y0=10.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #7 
x0=15.0*rr
y0=10.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #8
x0=7.5*rr
y0=8.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #9
x0=10.5*rr
y0=8.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #10
x0=13.5*rr
y0=8.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #11
x0=16.5*rr
y0=8.5*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #12
x0=9.0*rr
y0=13.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #13
x0=12.0*rr
y0=13.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #14
x0=15.0*rr
y0=13.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #15
x0=9.0*rr
y0=7.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #16 
x0=12.0*rr
y0=7.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

! Circle #17
x0=15.0*rr
y0=7.0*rr
do i=x0-r,x0+r
do j=y0-r,y0+r
d=((i-x0)*(i-x0))+((j-y0)*(j-y0))
if (d<=(r*r)) then
u(i,j)=0.0
v(i,j)=0.0
else
    continue
end if
end do
end do

do j=0,m
v(n,j)=0.0
end do

do i=0,n
do j=0,m
umean(i,j)=umean(i,j)+u(i,j)
ubar(i,j)=umean(i,j)/kk

vmean(i,j)=vmean(i,j)+v(i,j)
vbar(i,j)=vmean(i,j)/kk

uv(i,j)=u(i,j)*v(i,j)
uvmean(i,j)=uvmean(i,j)+uv(i,j)
uvbar(i,j)=uvmean(i,j)/kk

uprimvprimbar(i,j)=uvbar(i,j)-ubar(i,j)*vbar(i,j)

uu(i,j)=u(i,j)*u(i,j)
uumean(i,j)=uumean(i,j)+uu(i,j)
uubar(i,j)=uumean(i,j)/kk

uprim2bar(i,j)=uubar(i,j)-ubar(i,j)*ubar(i,j)

vv(i,j)=v(i,j)*v(i,j)
vvmean(i,j)=vvmean(i,j)+vv(i,j)
vvbar(i,j)=vvmean(i,j)/kk

vprim2bar(i,j)=vvbar(i,j)-vbar(i,j)*vbar(i,j)

kinetic(i,j)=0.5*rho(i,j)*(uprim2bar(i,j)+vprim2bar(i,j))
end do 
end do

print *,mstep-kk,ubar(13.5*rr,10.0*rr),vbar(13.5*rr,10.0*rr)!,error!,tau_total(n/2,m/2)
!pause
end do


write(2,*)"VARIABLES=I,J,U,V,Ubar,Vbar,UUbar,VVbar,UVbar,UprimVprimbar,Uprim2bar,Vprim2bar,Kinetic"
write(2,*)"ZONE ","I=",n+1,"J=",m+1,",","F=Block"
do j=0,m
write(2,*) (i,i=0,n)
end do
do j=0,m
write(2,*) (j,i=0,n)
end do
do j=0,m
write(2,*) (u(i,j),i=0,n)
end do
do j=0,m
write(2,*) (v(i,j),i=0,n)
end do
do j=0,m
write(2,*)(ubar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(vbar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(uubar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(vvbar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(uvbar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(uprimvprimbar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(uprim2bar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(vprim2bar(i,j),i=0,n)
end do
do j=0,m
write(2,*)(kinetic(i,j),i=0,n)
end do

do i=9.5*rr,11.5*rr
    write(3,*) (ubar(i,10.0*rr)/uo)
end do
do i=9.5*rr,11.5*rr
    write(4,*) (vbar(i,10.0*rr)/uo)
end do

do j=8.5*rr,11.5*rr
    write(5,*) (ubar(10*rr,j)/uo)
end do
do j=8.5*rr,11.5*rr
    write(6,*) (vbar(10*rr,j)/uo)
end do

do j=8.5*rr,11.5*rr
    write(7,*) (ubar(9.5*rr,j)/uo)
end do
do j=8.5*rr,11.5*rr
    write(8,*) (vbar(9.5*rr,j)/uo)
end do

do i=9.5*rr,11.5*rr
    write(9,*) (uvbar(i,10.0*rr)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(10,*) (uvbar(10*rr,j)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(11,*) (uvbar(9.5*rr,j)/(uo*uo))
end do

do i=9.5*rr,11.5*rr
    write(12,*) (kinetic(i,10.0*rr))
end do
do j=8.5*rr,11.5*rr
    write(13,*) (kinetic(10*rr,j))
end do
do j=8.5*rr,11.5*rr
    write(14,*) (kinetic(9.5*rr,j))
end do

do i=9.5*rr,11.5*rr
    write(15,*) (uubar(i,10.0*rr)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(16,*) (uubar(10*rr,j)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(17,*) (uubar(9.5*rr,j)/(uo*uo))
end do

do i=9.5*rr,11.5*rr
    write(18,*) (vvbar(i,10.0*rr)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(19,*) (vvbar(10*rr,j)/(uo*uo))
end do
do j=8.5*rr,11.5*rr
    write(20,*) (vvbar(9.5*rr,j)/(uo*uo))
end do
end program