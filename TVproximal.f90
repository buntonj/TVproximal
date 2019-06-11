module TVproximal
! MODULE FOR COMPUTING PROX OPERATORS FOR TV NORM
implicit none
contains

! 1-D PROXIMAL FOR L1 TV NORM
! Given a vector x in R^d, computes:
! minimize: 1/2 ||y-x||_2^2 + L||Dy||_1
! within input 'tol' of error
subroutine onedimTV1prox(x, y, tol,L, spew,iters)
integer :: n
logical, intent(in) :: spew
double precision, intent(inout), dimension(:) :: x
double precision, intent(in), dimension(:) :: y
double precision :: L, gap, tol, t, t_k_plus_one, c1 = 0.5
double precision :: f_k, f_k_plus_one, f_k_prime
double precision, dimension(:), allocatable :: diag, offdiag, u, gradient, descent_direction, sol
logical, dimension(:), allocatable :: inactive
integer :: i, j, k, info, iters, start,gradientcalc,rHessian,lapack,linesearch, project, gaptime, reduced_dim
iters = 0
n = size(x)
allocate(inactive(n-1), diag(n-1), u(n-1), gradient(n-1), offdiag(n-2))

! compute Dy for given y
do i = 1,n-1
	u(i) = y(i+1) - y(i)
end do

diag = 2.0d0
offdiag = -1.0d0

! solve DDT u = Dy
call dpttrf(n-1, diag, offdiag, info)
call dpttrs(n-1, 1, diag, offdiag, u, n-1, info)


! if all constraints already satisfied
if (all(dabs(u) .le. L)) then
	! compute corresponding primal optimal
	x = dual_to_primal(u, y)
	gradient = one_dim_TV_gradient(x)
	gap = L1_duality_gap(u,gradient,L)
	if (spew) then
		print*,'Initially correct'
		print*,'L1 TV Proximal'
		print*,'final gap:', gap,'Calculated with no iterations'
	end if
	deallocate(inactive, diag, gradient, u, offdiag)
	return

else
! Project back onto constraints
! check which constraints are satisfied
	inactive = (dabs(u) < L) 	!inactive identifies inactive constraints
	where(.not. inactive)
		u = sign(L,u)
	end where
end if

x = dual_to_primal(u, y)
gradient = one_dim_TV_gradient(x)
f_k = one_dim_TV_objective(x,y)
! compute initial duality gap
gap = L1_duality_gap(u,gradient,L)

allocate(descent_direction(n-1), sol(n-1))
! Iterate until convergence
! Projected Newton loop
do while (gap > tol)
	call system_clock(COUNT=start)
	iters = iters + 1
	
	!construct reduced Hessian and gradient at current u
	k = 0
	do i = 1, n-2
		if ( (inactive(i)) .or. (u(i)*gradient(i) > 0) ) then
			k = k + 1
			sol(k) = gradient(i)
			if (inactive(i+1)) then
				offdiag(k) = -1.0d0
			else
				offdiag(k) = 0.0d0
			end if
		end if
	end do
	if ((inactive(n-1)) .or. (u(n-1)*gradient(n-1) > 0)) then
		k = k + 1
		sol(k) = gradient(n-1)
	end if
	diag(1:k) = 2.0d0

	call system_clock(COUNT = rHessian)
	
	! solve Hd = gradient with reduced hessian for descent_direction d
	! Use structure to only factor Hessian once
	call dpttrf(k, diag(1:k), offdiag(1:(k-1)), info)
	call dpttrs(k, 1, diag(1:k), offdiag(1:(k-1)), sol(1:k), k, info)
	
	! map back from reduced space
	k=0
	do i = 1,n-1
		if ((inactive(i)) .or. (u(i)*gradient(i) > 0)) then
			k = k+1
			descent_direction(i) = -sol(k)
		else
			descent_direction(i) = 0.0
		end if
	end do

	call system_clock(COUNT=lapack)

	! begin backtracking search with quadratic interp about t = 1
	! attempt t = 1 step
	t = 1.0d0
	f_k_prime = dot_product(descent_direction, gradient) ! gradient wrt step size
	
	! Compute potential next step
	sol = u + t*descent_direction
	! project onto constraints
	inactive = (dabs(sol) < L)	! indicates inactive constraints
	where (.not. inactive)
		sol = sign(L,sol)
	end where
	x = dual_to_primal(sol,y)
	gradient = one_dim_TV_gradient(x)
	f_k_plus_one = one_dim_TV_objective(x,y)
	
	! if insufficient decrease
	if (f_k_plus_one > f_k + c1*t*f_k_prime) then
		t_k_plus_one =  0.5*(f_k_prime*t**2.0) / (f_k - f_k_plus_one + t*f_k_prime)	!set to quadratic interpolation minimizer
		sol = u + t_k_plus_one*descent_direction
		
		! identify inactive constraints and project
		inactive = (dabs(sol) < L)	! indicates inactive constraints
		where (.not. inactive)
			sol = sign(L,sol)
		end where
		x = dual_to_primal(sol,y)
		gradient = one_dim_TV_gradient(x)
		f_k_plus_one = one_dim_TV_objective(x,y)

		! if interpolated step is bigger than original, or insufficient decrease
		if ((t_k_plus_one > t)) then
			t_k_plus_one = t_k_plus_one * 0.5
			sol = u + t_k_plus_one*descent_direction
			
			! project on to constraints
			inactive = (dabs(sol) < L)	! indicates inactive constraints
			where (.not. inactive)
				sol = sign(L,sol)
			end where
			x = dual_to_primal(sol,y)
			gradient = one_dim_TV_gradient(x)
			f_k_plus_one = one_dim_TV_objective(x,y)
		end if
	else
		t_k_plus_one = t
	end if

	call system_clock(COUNT=linesearch)
	
	! take step of size t_k_plus_one in descent_direction
	u = sol
	f_k = f_k_plus_one
	! compute duality gap
	gap = L1_duality_gap(u,gradient,L)
	
	call system_clock(COUNT = gaptime)
	if (spew) then
		print*,'iter:',iters,'gap:',gap,'step:',t_k_plus_one
	end if
end do

if (spew) then
	print*,'L1 TV Proximal'
	print*,'Final gap:', gap,'Calculated in', iters, 'iterations'
end if
!if (iters > 0) then
!	print*,'TIMES:'
!	print*,'Computing gradient:', gradientcalc - start
!	print*,'Computing reduced Hessian:', rhessian -gradientcalc
!	print*,'Computing w/ LAPACK:', lapack - gradientcalc
!	print*,'Computing linesearch:', linesearch - lapack
!	print*,'Computing projection:', project - linesearch
!	print*,'Computing gap:', gaptime - project
!	print*,'TOTAL:', gaptime - start
!end if
deallocate(inactive, diag, offdiag, u, gradient, descent_direction, sol)
return
end subroutine


! 1-D PROXIMAL FOR L2 TV NORM
! Given a vector x in R^d, computes:
! minimize: 1/2 ||y-x||_2^2 + L||Dy||_2
! within input 'tol' of error
subroutine onedimTV2prox(x, y, tol, L, spew, iters)
double precision, dimension(:), intent(in) :: y
double precision, dimension(:), intent(out) :: x
double precision, intent(in) :: L, tol
double precision, allocatable, dimension(:) :: gradient, diag, offdiag, Dy, q, u
integer :: i, iters, n, info
double precision :: gap, alpha, normu, normq
logical :: spew
iters = 0
n = size(y,1)
allocate(diag(n-1), offdiag(n-2), gradient(n-1), Dy(n-1), u(n-1), q(n-1))

alpha = 0.0d0
gap = tol + 1.0

do i = 1,n-1
	Dy(i) = y(i+1) - y(i)
end do

! compute DD^T + alpha*I
diag(:) = 2.0d0 + alpha
offdiag(:) = -1.0d0

! compute cholesky factorization
call dpttrf(n-1, diag, offdiag, info)

! Solves (DD^T + alpha*I)u = Dy
u = Dy
call dpttrs(n-1,1,diag,offdiag,u,n-1,info)
open(unit=111,file='tv2convergence.txt')
gap = tol + 1.0d0
do while (gap > tol)
	iters = iters + 1
	! compute DD^T + alpha*I
	diag(:) = 2.0d0 + alpha
	offdiag(:) = -1.0d0

	! compute cholesky factorization
	call dpttrf(n-1, diag, offdiag, info)
	
	! Solves (DD^T + alpha*I)u = Dy
	u = Dy
	call dpttrs(n-1,1,diag,offdiag,u,n-1,info)

	! Compute cholesky matrix R^T
	do i = 1,n-2
		diag(i) = sqrt(diag(i))
		offdiag(i) = offdiag(i)*diag(i)
	end do
	diag(n-1) = sqrt(diag(n-1))

	! Solve R^Tq = u
	q(1) = u(1)/diag(1)
	do i = 2,n-1
		q(i) = (u(i) - offdiag(i)*q(i-1))/diag(i)
	end do
	
	normu = norm2(u)
	normq = norm2(q)

	alpha = max(0.0d0,alpha - (1.0d0 - normu/L)*(normu/normq)**2.0d0)

	! compute DD^T + alpha*I
	diag(:) = 2.0d0 + alpha
	offdiag(:) = -1.0d0

	! compute cholesky factorization
	call dpttrf(n-1, diag, offdiag, info)

	! Solves (DD^T + alpha*I)u = Dy
	u = Dy
	call dpttrs(n-1,1,diag,offdiag,u,n-1,info)
	
	! check stopping condition
	x = dual_to_primal(u,y)
	gradient = one_dim_TV_gradient(x)
	gap = L2_duality_gap(u,gradient,L)
	!print*,'gap:',gap
	write(111,*) gap
end do
close(111)
if (spew) then
	print*,'L2 TV Proximal'
	print*,'final gap:', gap,'computed in',iters,'iterations'
end if
deallocate(diag, offdiag, gradient, Dy, u, q)

end subroutine



! SUBROUTINE TO COMPUTE 2-D PROX of L1 TV NORM
! Given a (m x n) matrix Y
! solves:   minimize ||X-Y||_F^2 + L||X|_{tv1}
subroutine twodimTV1prox(X,Y,tol,L,spew)
double precision, dimension(:,:), intent(in) :: Y
double precision, dimension(:,:), intent(inout) :: X
double precision :: tol, L, mean_change, t
double precision, dimension(:,:), allocatable :: P, Q, Z, last
integer :: m, n, i, fill
logical :: spew

m = size(X,1)
n = size(X,2)
allocate(P(m,n), Q(m,n), Z(m,n),last(m,n))

X = Y
P = 0
Q = 0
t = 0
mean_change = 1
do while (mean_change > tol)
	last = X
	do i = 1,m
		call onedimTV1prox(Z(i,:), X(i,:) + P(i,:), tol, L, .false., fill)
	end do
	P = X + P - Z
	do i = 1,n
		call onedimTV1prox(X(:,i), Z(:,i) + Q(:,i), tol, L, .false., fill)
	end do
	Q = Z + Q - X
	last = dabs(last - X)
	mean_change = sum(sum(last,dim=1))/(m*n)
	if (spew) then
		print*,'mean change:',mean_change
	end if
	t = t + 1
end do
if (spew) then
	print*,'Converged in', t,'steps'
end if
deallocate(P,Q,Z,last)

end subroutine


! SUBROUTINE TO COMPUTE 2-D PROX of L2 TV NORM
! Given a (m x n) matrix Y
! solves:   minimize ||X-Y||_F^2 + L||X||_{tv2}
subroutine twodimTV2prox(X,Y,tol,L,spew)
double precision, dimension(:,:), intent(in) :: Y
double precision, dimension(:,:), intent(inout) :: X
double precision :: tol, L, mean_change, t
double precision, dimension(:,:), allocatable :: P, Q, Z, last
logical :: spew
integer :: m, n, i, fill

m = size(X,1)
n = size(X,2)
allocate(P(m,n), Q(m,n), Z(m,n),last(m,n))

X = Y
P = 0
Q = 0
t = 0
mean_change = 1
do while (mean_change > 10D-5)
	last = X
	do i = 1,m
		call onedimTV2prox(Z(i,:), X(i,:) + P(i,:), tol, L, .false.,fill)
	end do
	P = X + P - Z
	do i = 1,n
		call onedimTV2prox(X(:,i), Z(:,i) + Q(:,i), tol, L, .false.,fill)
	end do
	Q = Z + Q - X
	last = dabs(last - X)
	mean_change = sum(sum(last,dim=1))/(m*n)
	!if (spew) then
	!	print*,'mean change:',mean_change
	!end if
	t = t + 1
end do
if (spew) then
	print*,'Converged in', t,'steps'
end if

deallocate(P,Q,Z,last)

end subroutine





! FUNCTION TO COMPUTE L2 TV DUALITY GAP
! INPUT:
! dual variable u
! gradient g
! TV norm coefficient L
! returns duality gap between current dual and primal
double precision function L2_duality_gap(u,g,L)
double precision, dimension(:) :: g, u
double precision :: L

L2_duality_gap = dabs(L*norm2(g) + dot_product(u,g))

return
end function




! FUNCTION TO COMPUTE L1 TV DUALITY GAP
! INPUT:
! dual variable u
! current gradient g
! TV norm coefficient L
! Returns duality gap between primal and dual results
double precision function L1_duality_gap(u,g,L)
double precision, dimension(:) :: g, u
double precision :: L
L1_duality_gap = L*sum(dabs(g)) + dot_product(u,g)
return
end function



! FUNCTION TO COMPUTE OBJECTIVE VALUE OF 1D DUAL
! INPUT:
! current primal variable x
! TV norm input y
double precision function one_dim_TV_objective(x,y)
double precision, dimension(:) :: x,y

one_dim_TV_objective = 0.5*(dot_product(x,x) - dot_product(y, y))

return
end function



! FUNCTION TO COMPUTE GRADIENT OF OBJECTIVE FUNCTION OF 1D DUAL
! INPUT:
! current primal variable x
function one_dim_TV_gradient(x)
double precision, dimension(:) :: x
integer ::i,n
double precision, dimension(size(x,1)-1) :: one_dim_TV_gradient

n = size(x,1)

do i = 1,n-1
	one_dim_TV_gradient(i) = x(i) - x(i+1)
end do
return
end function



! FUNCTION TO COMPUTE CORRESPONDING PRIMAL x for 1D TV DUAL
! INPUT:
! current dual variable u
! TV norm input y
function dual_to_primal(u,y)
double precision, dimension(:) :: u, y
integer :: i, n
double precision, dimension(size(y,1)) :: dual_to_primal

n = size(y,1)

dual_to_primal(1) = y(1) + u(1)
do i = 2,n-1
	dual_to_primal(i) = y(i) - u(i-1) + u(i)
end do
dual_to_primal(n) = y(n) - u(n-1)
return
end function


function gradient_calc(X,Y,K)
	double precision, dimension(:,:) :: X, Y, K
	double precision, dimension(size(X,1),size(X,2)) :: gradient_calc
	gradient_calc = matmul(K,X) - Y
	return
end function

! given parameters mu and sigma, computes a sample from
! normal(mu,sigma)
double precision function normal_sample(mu, sigma)
	double precision :: rand1, rand2, mu, sigma

	call random_number(rand1)
	call random_number(rand2)
	
	! given (0,1)
	normal_sample = sigma*sqrt(-2*log(rand1))*cos(2*acos(-1.0)*rand2) + mu
return
end function


end module
