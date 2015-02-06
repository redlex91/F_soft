! modules
module prec_def
implicit none
save
integer, parameter :: long = selected_real_kind( 18,307 )
end module prec_def

module constants
use prec_def
implicit none
save
integer, parameter :: NOS = 8, MAX_ITER = 10
real( long ), parameter :: rho = 0.1, E = 1, sigma = 1
real( long ), parameter :: r_c = 2.5, F_c = -0.039, u_c = -0.0163
real(long), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_long, ETA_M = PI/4
end module constants




program soft_spheres

use prec_def
use constants

! VARIABLE DECLARATION
implicit none
real( long ) :: vol, v_mean, r;
real( long ), dimension( 1:NOS ) :: rx, ry, rz, vx, vy, vz, fx, fy, fz
real( long ) :: distance, comp ! delcaration of function distance
integer :: i, j, iter
integer, parameter :: SCALE_STEP = 10

interface

subroutine init_cond( x, y, z, vx, vy, vz )
use prec_def
use constants
implicit none
real( long ), dimension( 1: NOS ), intent( out ) :: x, y, z, vx, vy, vz
end subroutine init_cond

subroutine rescale_temp( d_temp, v_mean, vx, vy, vz )
	use prec_def
	use constants
	implicit none
	real( long ), intent( in ) :: d_temp, v_mean
	real( long ), dimension( 1:NOS ), intent( inout ) :: vx, vy, vz
endsubroutine rescale_temp

end interface 

vol = NOS / rho
call init_cond( rx, ry, rz, vx, vy, vz )
open( unit = 1, file = "posit.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of energies

do i= 1, NOS
	write( unit = 1, fmt = 1 ) rx( i ), ry( i ) , rz( i )
		1 format( f20.15, f20.15, f20.15 )
enddo

step: do iter = 1, MAX_ITER

	v_mean = ( dot_product( vx, vx ) + dot_product( vy, vy ) + dot_product( vz, vz ) ) / SCALE_STEP
	if( mod( iter, SCALE_STEP ) == 0 ) then
		call rescale_temp( 298.0_long, v_mean, vx, vy, vz )
		v_mean = 0
		print *, iter
	endif

	fx = 0; fy = 0; fz = 0
	do i = 1, NOS
		do j = 1, NOS
			r = distance( rx( i ), ry( i ), rx( i ), rx( j ), ry( j ), rz( j ) )
			if(  r <= r_c .and. i.ne.j ) then
				fx = fx + 24 * comp( rx( i ), rx( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
				fy = fy + 24 * comp( ry( i ), ry( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
				fz = fz + 24 * comp( rz( i ), rz( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
			end if
		enddo
	end do
enddo step


close( unit = 1 )
end program soft_spheres

subroutine init_cond( x, y, z, vx, vy, vz )
use prec_def
use constants
implicit none
real( long ), dimension( 1: NOS ), intent( out ) :: x, y, z, vx, vy, vz
integer :: n
real( long ) :: side, step
integer :: i

! var for random values
real( long ) :: randnum, K, vxc, vyc, vzc ! K =n kinetic energy per particle
integer, dimension( : ), allocatable:: seed
integer, dimension( 1:8 ) :: dt_seed
integer :: n_seed
n = nint( NOS**(1.0/3.0) )
side = ( NOS / rho )**(1.0/3.0)
step = side / n
! assigning positions
do i = 1, NOS
	x( i ) = step * ( 1/2.0 + mod( i - 1 , n ) )
	y( i ) = step * ( 1/2.0 + mod( ( i - 1 ) / n , n ) )
	z( i ) = step * ( 1/2.0 + ( i - 1 ) / (n*n) )
enddo

! setting seed for speen generator
call random_seed( size = n_seed )
allocate( seed( 1:n_seed ) )
call random_seed( get = seed )
call date_and_time( values = dt_seed )
seed( n_seed ) = dt_seed( 8 ); seed( 1 ) = dt_seed( 8 ) * dt_seed( 7 ) * dt_seed( 6 )
call random_seed( put = seed )
deallocate( seed )
! done setting seed

! assigning random velocities
do i= 1, NOS
	call random_number( randnum )
	vx( i ) = 2*randnum - 1
end do
do i= 1, NOS
	call random_number( randnum )
	vy( i ) = 2*randnum - 1
end do
do i= 1, NOS
	call random_number( randnum )
	vz( i ) = 2*randnum - 1
end do

! computing the cm velocitiy
vxc = 0; vyc = 0; vzc = 0
do i = 1, NOS
	vxc = vxc + vx( i ) / NOS
	vyc = vyc + vy( i ) / NOS
	vzc = vzc + vz( i ) / NOS
end do
! subtracting the cm velocity
do i = 1, NOS
	vx( i ) = vx( i ) -vxc
	vy( i ) = vy( i ) - vyc
	vz( i ) = vz( i ) - vzc
end do
! verifying cm velocity = 0
vxc = 0; vyc = 0; vzc = 0
do i = 1, NOS
	vxc = vxc + vx( i )
	vyc = vyc + vy( i )
	vzc = vzc + vz( i )
end do

K = ( dot_product( x, x ) + dot_product( y, y ) + dot_product( z, z ) ) / ( 2.0 * NOS )

end subroutine init_cond

function distance( x1, y1, z1, x2, y2, z2 )
use prec_def
implicit none
real( long ), intent( in ) :: x1, x2, y1, y2, z1, z2
real( long ) :: distance
! local var
real( long ) :: a1, a2, b1, b2, c1, c2
a1 = x1; a2 = x2; b1 = y1; b2 = y2; c1 = z1; c2 = z2

if( abs( a1-a2 ) > 0.5 ) a2 = a2 - sign( 1.0_long, a2 - a1 )
if( abs( b1-b2 ) > 0.5 ) b2 = b2 - sign( 1.0_long, b2 - b1 )
if( abs( c1-c2 ) > 0.5 ) c2 = c2 - sign( 1.0_long, c2 - c1 )

distance = sqrt( (a2 - a1 )**2 + ( b2 - b1 )**2 + ( c2 - c1 )**2 )

end function distance

function comp( q1, q2 )
use prec_def
implicit none
real( long ), intent( in ) :: q1, q2
real( long ) :: comp
! local var
real( long ) :: a1, a2, b1, b2, c1, c2
a1 = q1; a2 = q2

if( abs( a1-a2 ) > 0.5 ) a2 = a2 - sign( 1.0_long, a2 - a1 )
comp = a2 - a1

end function comp

subroutine rescale_temp( d_temp, v_mean, vx, vy, vz )
	use prec_def
	use constants
	implicit none
	real( long ), intent( in ) :: d_temp, v_mean
	real( long ), dimension( 1:NOS ), intent( inout ) :: vx, vy, vz
	! local var
	real( long ) :: a_temp
	integer :: k
	
	a_temp =  v_mean / ( 2 * NOS )
	
	vx = vx * sqrt( d_temp / a_temp )
	vy = vy * sqrt( d_temp / a_temp )
	vz = vz * sqrt( d_temp / a_temp )

end subroutine rescale_temp
