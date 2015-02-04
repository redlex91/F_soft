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
integer, parameter :: NOS = 15625
real( long ), parameter :: rho = 0.1, E = 1, sigma = 1
real(long), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_long, ETA_M = PI/4
end module constants




program soft_spheres

use prec_def
use constants

! VARIABLE DECLARATION
implicit none
real( long ) :: vol;
real( long ), dimension( 1:NOS ) :: rx, ry, rz, vx, vy, vz
integer :: i

interface

subroutine init_cond( x, y, z, vx, vy, vz )
use prec_def
use constants
implicit none
real( long ), dimension( 1: NOS ), intent( out ) :: x, y, z, vx, vy, vz
end subroutine init_cond

end interface 

vol = NOS / rho
call init_cond( rx, ry, rz, vx, vy, vz )
open( unit = 1, file = "posit.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of energies

do i= 1, NOS
	write( unit = 1, fmt = 1 ) rx( i ), ry( i ) , rz( i )
		1 format( f20.15, f20.15, f20.15 )
enddo

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
real( long ) :: randnum, K, vxc, vyc, vzc
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
end subroutine init_cond

