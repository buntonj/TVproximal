program inexactproxgradient
use TVproximal

implicit none
double precision, dimension(:,:), allocatable :: x, y, K, grd, z, x_k_plus_one
!double precision, dimension(:), allocatable :: x, y, K, grd, z, x_k_plus_one
integer :: i, j, m, n, iters = 1000, numiterations
logical :: spew = .true.
double precision :: tol, L, interval, start, finish, t

!n = 10000
! PENGUIN DIMS
m = 330
n = 220

! CHESS DIMS
!m = 600
!n = 600
spew = .false.
call random_seed()

! Error tolerance
tol = 10D-5
L = 1.0

!allocate(x(m,n), x_k_plus_one(m,n), y(m,n), K(m,m), z(m,n), grd(m,n))

allocate(x(m,n),y(m,n))

!open(unit=500, file = 'tests.txt')
!do i = 3,20
!	n = nint(10**(float(i)/3))
!	allocate(x(n),y(n))
!	call random_number(L)
!	L = nint(L*50)
!	call random_number(y)
!	y = 2.0*L*(y-0.5)
!	! place y in interval [-2l, 2l]
!	call cpu_time(start)
!	call onedimtv1prox(x,y,tol,L,.false.,numiterations)
!	deallocate(x,y)
!	call cpu_time(finish)
!	print*,n
!	write(500,*) n, finish-start, numiterations
!end do
!close(500)

open(unit=50, file = 'graypenguin.txt')
do i = 1,m
	read(50,*) (Y(i,j), j =1,n)
end do
close(50)

!K(:,:) = 0.0d0

!do i = 1,m-3
!	K(i:(i+3),i:(i+3)) = 1.0d0/8.0
!end do
!Y = matmul(K,Y)
! initialize a signal vector
do i = 1,m
	do j = 1,n
		Y(i,j) = max(0.0,y(i,j) + normal_sample(0.0d0,0.3d0))
	end do
end do

!do i = 1,n
!	y(i) = sin(float(i)/10) + normal_sample(0.0d0, 0.2d0)
!end do

! initialize

!call onedimTV2prox(x,y,tol,L,spew)

! PROXIMAL GRADIENT METHOD
!open(unit=60,file='convergence.txt')
!t = 0.01d0
!do j = 1,50
!	tol = 1.0/(float(j)**2.0)
!	grd = 7.0*gradient_calc(X,Y,K)
!	X_k_plus_one = X - t*GRD
!	call twodimTV1prox(X, X_k_plus_one, tol, t, .false.)
!	!z = x + float(j-1)/float(j+2)*(x - x_k_plus_one)
!	write(60,*) sum(sum(X_k_plus_one-X,1))
!	print*,j
!end do
!close(60)
!do i = 1,60
!	print*,Y(220,i),','
!end do
L = 1.9
call twodimTV1prox(x,y,tol,L,spew)



open(unit=100,file='result.txt')
open(unit=200,file='truth.txt')
!do i = 1,n
!	write(100,*) x(i)
!	write(200,*) y(i)
do i = 1,m
	write(100,*) (x(i,j), j = 1,n)
	write(200,*) (y(i,j), j = 1,n)
end do
close(100)
close(200)

deallocate(x,y)
!deallocate(X,Y,K,grd,z,x_k_plus_one)

end program
