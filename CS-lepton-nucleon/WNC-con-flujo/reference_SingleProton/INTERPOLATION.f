
c   R. Gonzalez-Jimenez --  26 Junio 2012

C Este programa hace una interpolacion a partir de puntos experimentales (un archivo externo, dos columnas: 'x' e 'y').
c El numero de puntos que añade entre cada dos puntos verdaderos viene dado por el parametro Finura de la siguiente forma:
c Finura representa el 'Delta x' que queremos establecer, por ejemplo, si entre cada dos puntos experimentales tenemos  
c un 'Delta x = 0.050', podriamos fijar 'Finura=0.001', con lo que estariamos introduciendo 49 nuevos puntos entre los dos 
c puntos experimentales. La interpolacion se hace con la ecuacion de una recta, asi que meter mas de 4 o 5 puntos es 
c estupido. 
	program Interpolation
	implicit real*8(a-h,p-z)
	real*8 x(300), y(300)
	real*8 xm(300), xn(300)
	character*40 ipivot

	open(16,file='flujo_interpolado.in')
	open(234,file='flujo.in',status='old')
	read(234,*)ipivot
	nmax = 0
	do 9879 ifi = 2,300
	read(234,*,end=9889)x(ifi),y(ifi)
	nmax = ifi
9879	continue
9889	continue
	print*,'nmax=',nmax

	x(1) = 0.d0
	y(1) = 0.d0

	do 10 i = 2,nmax
	xm(i-1) = (y(i) - y(i-1)) / (x(i) - x(i-1))
	xn(i-1) = y(i-1) - xm(i-1) * x(i-1)
10	continue

	Finura = 0.005
	do 300 ii=2,nmax
	Dx = 0.d0

	    do 20 jj = 1,1000
            xx = x(ii-1) + Dx
	    if(xx.le.x(ii).and.xx.gt.x(ii-1))then
	    yy = xm(ii-1)*xx + xn(ii-1)
	    write(16,*)xx,yy
	    endif
	    if(xx.gt.x(ii))goto 300
	    Dx = Dx + Finura
20	    continue

300	continue

	stop
	end 

