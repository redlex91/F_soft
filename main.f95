! modules
MODULE prec_def
  ! purpose: this module defines the precision on REAL type
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: long = selected_real_kind( 18,307 )
ENDMODULE prec_def

MODULE constants
  ! purpose: declare all the constants used in the program,
  !          so that they are available to both the main
  !          and all the subroutines
USE prec_def
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: NOS = 1000, MAX_ITER = 10 ! NOS = number of spheres, MAX_ITER = maximum number of iterations of the program
REAL( long ), PARAMETER :: rho = 0.7, deltat = 0.001, r_c = 2.5, r_L = 2.8, F_c = -0.039, u_c = -0.0163 ! density, time interval, core distance, effective core distance, force and potential energy at r_c
REAL(long), PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510_long, ETA_M = PI/4
endmodule constants


MODULE dynamic_list
 
  ! Purpose: it implements a dynamic list, whose dimension in a priori unknown,
  !          and all the subroutines needed to madeal with it (add node, print,
  !          destroy the entire list after use
  
  IMPLICIT NONE
  
  TYPE :: node
     INTEGER :: data 
     TYPE( node ), POINTER :: nxtPtr
  ENDTYPE node
  
CONTAINS

  SUBROUTINE add_node( headPtr, tailPtr, val )
    ! Purpose: this sbrt adds a new node at the end of the list (queue-type).
    IMPLICIT NONE

    TYPE( node ), POINTER :: headPtr, tailPtr
    TYPE( node ), POINTER :: newPtr => NULL()
    INTEGER, INTENT( IN ) :: val
    INTEGER :: flag ! flag == 0 iff allocation is ok

    ALLOCATE( newPtr, STAt = flag ) ! allocating memory of the new node
    IF( flag == 0 ) THEN ! allocation was ok
       newPtr%data = val
       NULLIFY( newPtr%nxtPtr )
       
       IF( .NOT. ASSOCIATED ( headPtr ) ) THEN ! if list is empty
          headPtr => newPtr 
       ELSE
          tailPtr%nxtPtr => newPtr
       ENDIF
       tailPtr => newPtr
    ELSE ! allocations was KO
       WRITE( *, * ) 'ERROR! ADDING NEW NODE HAS FAILED'
    ENDIF
  ENDSUBROUTINE add_node

  SUBROUTINE kill_list( headPtr, tailPtr )
    ! Purpose: this sbrt destroys the entire list after use in order to free
    !          memory

    IMPLICIT NONE
    TYPE( node ), POINTER :: headPtr, tailPtr
    TYPE( node ), POINTER :: currPtr => NULL()
    INTEGER :: flag = 0

    currPtr => headPtr
    DO
       IF( .NOT. ASSOCIATED( headPtr ) ) EXIT ! list is already empty
       currPtr => currPtr%nxtPtr
       DEALLOCATE( headPtr, STAT = flag )
       IF ( flag == 0 ) THEN
          headPtr => currPtr
       ELSE
          WRITE( *, * ) 'ERROR! DEALLOCATION HAS FAILED!'
       ENDIF
    ENDDO
    NULLIFY( tailPtr )
    WRITE( *, * ) 'List is NOW empty.'
    
  ENDSUBROUTINE kill_list

  SUBROUTINE print_list( headPtr, tailPtr )
    IMPLICIT NONE
    TYPE( node ), POINTER :: headPtr, tailPtr, currPtr => NULL()

    currPtr => headPtr

    IF( .NOT. ASSOCIATED( headPtr ) ) THEN ! list is empty
       WRITE( *, * ) 'List is empty.'
    ELSE
       WRITE( *, * ) 'Printing list:'
       DO
          WRITE( *, * ) currPtr%data
          IF ( ASSOCIATED( currPtr, tailPtr ))  EXIT ! reached end of list
          currPtr => currPtr%nxtPtr
       ENDDO
       WRITE( *, * ) 'End of list.'
    ENDIF
          
  ENDSUBROUTINE print_list
  
ENDMODULE dynamic_list

!*****************************************************************************************************************************************************************************************************************
!   THE MASSES, SIGMA AND EPSILON ARE ALL NORMALISED TO 1 SO THAT ALL DISTANCES ARE MEASURED IN UNITS OF SIGMA AND ALL ENERGIES ARE MEASURED IN UNITS OF EPSILON
!*****************************************************************************************************************************************************************************************************************
!*****************************************************************************************************************************************************************************************************************
!                                                         START OF MAIN PROGRAM
!*****************************************************************************************************************************************************************************************************************
!*****************************************************************************************************************************************************************************************************************

PROGRAM soft_spheres

use prec_def
use constants

! VARIABLE DECLARATION
IMPLICIT NONE
REAL( long ) :: v_mean, r
REAL( long ), dimension( 1:NOS ) :: rx, ry, rz, vx, vy, vz, fx, fy, fz 
REAL( long ) :: distance, comp ! delcaration of function distance
INTEGER :: i, j, iter
INTEGER, PARAMETER :: SCALE_STEP = 10

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

! vol = NOS / rho ! the volume is fixed by the densisty

! STEP 1: INITIALISATION

call init_cond( rx, ry, rz, vx, vy, vz ) ! first initialization of the system, no total energy / temperature assumed

open( unit = 1, file = "posit.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of positions

DO i = 1, NOS
   WRITE( unit = 1, fmt = 1 ) rx( i ), ry( i ) , rz( i )
   1 format( f20.15, f20.15, f20.15 )
ENDDO

!step: do iter = 1, MAX_ITER

!	v_mean = ( dot_product( vx, vx ) + dot_product( vy, vy ) + dot_product( vz, vz ) ) / SCALE_STEP
!	if( mod( iter, SCALE_STEP ) == 0 ) then
!		call rescale_temp( 298.0_long, v_mean, vx, vy, vz )
!		v_mean = 0
!		print *, iter
!	endif

!	fx = 0; fy = 0; fz = 0
!	do i = 1, NOS
!		do j = 1, NOS
!			r = distance( rx( i ), ry( i ), rx( i ), rx( j ), ry( j ), rz( j ) )
!			if(  r <= r_c .and. i.ne.j ) then
!				fx = fx + 24 * comp( rx( i ), rx( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
!				fy = fy + 24 * comp( ry( i ), ry( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
!				fz = fz + 24 * comp( rz( i ), rz( j ) ) * ( 2 * ( 1/r )**14 - ( 1/r )**8 + F_c / 24 )
!			end if
!		enddo
!	end do
!enddo step

close( unit = 1 )
end program soft_spheres

!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************
!                                                                                       MAIN PROGRAM ENDS HERE!
!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************


subroutine init_cond( x, y, z, vx, vy, vz )
  ! purpose: this subroutine initialises the system, assigning positions on a regular lattice and random velocities 
use prec_def
use constants
implicit none
real( long ), dimension( 1: NOS ), intent( out ) :: x, y, z, vx, vy, vz
integer :: n
real( long ) :: side, step
integer :: i

! variables for generating random values
real( long ) :: randnum, K, vxc, vyc, vzc ! K = kinetic energy per particle
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


!!! COORECT 0.5 WITH 0.5 * SIDE AND SO ON
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
  ! purpose: 
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
  ! purpose: this subroutine rescales the velocities, given a desired value of the temperature
  !          it takes as inputs the desired temperature and ???
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
