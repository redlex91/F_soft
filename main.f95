! Declaration of modules
MODULE prec_def
  ! Purpose: this module defines the precision on REAL type
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: long = selected_real_kind( 18,307 )
ENDMODULE prec_def

MODULE constants
  ! Purpose: declare all the constants used in the program,
  !          so that they are available to both the main
  !          and all the subroutines
USE prec_def
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: NOS = 9, MAX_ITER = 10 ! NOS = number of spheres, MAX_ITER = maximum number of iterations of the program
REAL( long ), PARAMETER :: rho = 0.3, deltat = 0.001, r_c = 2.5, r_L = 2.8, F_c = -0.039, u_c = -0.0163 ! density, time interval, core distance, effective core distance, force and potential energy at r_c
REAL(long), PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510_long, ETA_M = PI/4
END MODULE constants

MODULE global

  ! Purpose: this module contains all the global variables and the general-aim functions and procedures.

  USE prec_def
  USE constants

  IMPLICIT NONE
  SAVE
  REAL( long ) :: side, step ! side of the cube and initial separation of the particles
  INTEGER :: n ! cube root of NOS rounded to the nearest integer

CONTAINS

  FUNCTION distance( x1, y1, z1, x2, y2, z2 )

  ! Purpose: it takes as arguments the coordinates of two particles and returns the distance between them,
  !          accounting for the boundary conditions

  IMPLICIT NONE
  REAL( long ), INTENT( in ) :: x1, x2, y1, y2, z1, z2
  REAL( long ) :: distance
  ! local var
  REAL( long ) :: a1, a2, b1, b2, c1, c2
  a1 = x1; a2 = x2; b1 = y1; b2 = y2; c1 = z1; c2 = z2

  IF( ABS( a1-a2 ) > 0.5*side ) a2 = a2 - SIGN( 1._long, a2 - a1 )*side
  IF( ABS( b1-b2 ) > 0.5*side ) b2 = b2 - SIGN( 1._long, b2 - b1 )*side
  IF( ABS( c1-c2 ) > 0.5*side ) c2 = c2 - SIGN( 1._long, c2 - c1 )*side

  distance = sqrt( (a2 - a1 )**2 + ( b2 - b1 )**2 + ( c2 - c1 )**2 )
  RETURN
  
END FUNCTION distance

FUNCTION comp( q1, q2 )
  
  ! Purpose: this function computes the difference between the same components of a vector
  
  IMPLICIT NONE
  REAL( long ), INTENT( in ) :: q1, q2
  REAL( long ) :: comp
  ! local var
  REAL( long ) :: a1, a2, b1, b2, c1, c2
  a1 = q1; a2 = q2

  IF( ABS( a1-a2 ) > 0.5*side ) a2 = a2 - SIGN( 1.0_long, a2 - a1 )*side
  comp = a2 - a1
  RETURN

END FUNCTION comp

FUNCTION find_next( current, array )

  ! Purpose: it looks for the next non-NULL pointer in the npoint array from
  !          the current position and returns its position; if no such pointer
  !          can be found it returns NOS

  USE prec_def
  USE constants
  USE dynamic_list
  
  IMPLICIT NONE
  TYPE( ptr2node ), DIMENSION( 1: NOS ), INTENT( in ) :: array
  INTEGER, INTENT( in ) :: current
  INTEGER :: find_next
  INTEGER :: k

  find_next = NOS
  DO k = current + 1, NOS
     IF( ASSOCIATED( array( k )%Ptr ) ) THEN
        find_next = k
        EXIT
     ENDIF
  ENDDO
  RETURN

ENDFUNCTION find_next

FUNCTION force( x1, x2, r )

  ! Purpose: this function returns the value of the force along one component

  USE prec_def
  USE constants
  
  IMPLICIT NONE
  REAL( long ), INTENT( in ) :: x1, x2, r
  REAL( long ) :: force

  IF( r <= r_c ) THEN
     force = 24. * ( 2.*(1./r)**(13.) - (1./r)**(7.) )* comp( x1, x2 ) / r
  ELSE
     force = 0
  ENDIF
  RETURN
  
ENDFUNCTION force

ENDMODULE global

MODULE init

  ! Purpose: it contains all the subroutines dedicated to the initialisation of
  ! the system

  USE prec_def
  USE constants
  USE global
 
CONTAINS
  
SUBROUTINE init_cond( x, y, z, vx, vy, vz )
  ! Purpose: this subroutine initialises the system, assigning positions on a regular lattice and random velocities 

IMPLICIT NONE
real( long ), dimension( 1: NOS ), intent( out ) :: x, y, z, vx, vy, vz
integer :: i

! variables for generating random values
real( long ) :: randnum, K, vxc, vyc, vzc ! K = kinetic energy per particle
integer, dimension( : ), allocatable:: seed
integer, dimension( 1:8 ) :: dt_seed
integer :: n_seed


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

END SUBROUTINE init_cond

SUBROUTINE rescale_temp( d_temp, v_mean, vx, vy, vz )
  ! Purpose: this subroutine rescales the velocities, given a desired value of the temperature
  !          it takes as inputs the desired temperature and ???

  IMPLICIT NONE
  real( long ), intent( in ) :: d_temp, v_mean
  real( long ), dimension( 1:NOS ), intent( inout ) :: vx, vy, vz
  ! local var
  real( long ) :: a_temp
  integer :: k
	
	a_temp =  v_mean / ( 2 * NOS )
	
	vx = vx * sqrt( d_temp / a_temp )
	vy = vy * sqrt( d_temp / a_temp )
	vz = vz * sqrt( d_temp / a_temp )

END SUBROUTINE rescale_temp

END MODULE init

MODULE dynamic_list
 
  ! Purpose: it implements a dynamic list, whose dimension in a priori unknown,
  !          and all the subroutines needed to madeal with it (add node, print,
  !          destroy the entire list after use
  
  IMPLICIT NONE
  
  TYPE :: node
     INTEGER :: data 
     TYPE( node ), POINTER :: nxtPtr
  ENDTYPE node

  TYPE ptr2node
   ! I define a pointer to the node type in order to define an "array of pointers"
     TYPE( node ), POINTER :: Ptr
  END TYPE ptr2node

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

!****************************************************************************************************************************************************************************************************************
!   THE MASSES, SIGMA AND EPSILON ARE ALL NORMALISED TO 1 SO THAT ALL DISTANCES ARE MEASURED IN UNITS OF SIGMA AND ALL ENERGIES ARE MEASURED IN UNITS OF EPSILON
!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************
!                                                         START OF MAIN PROGRAM
!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************

PROGRAM soft_spheres

USE prec_def
USE constants
USE init
USE dynamic_list
USE global

! VARIABLE DECLARATION
IMPLICIT NONE
REAL( long ) :: v_mean, r, time
REAL( long ), dimension( 1:NOS ) :: rx, ry, rz, vx, vy, vz, fx, fy, fz, ax, ay, az 

INTEGER :: i, j, iter
INTEGER :: flag = 0

TYPE( ptr2node ) :: listH, listT  ! these variables contain the pointers to the head and the tail of the "list", see notes
TYPE( ptr2node ), DIMENSION( 1:NOS ) :: npoint ! same as in the notes but instead of containing the position inside the list it points to the first occurence
TYPE( ptr2node ) :: tmp ! this pointer is needed to skim the list

! vol = NOS / rho ! the volume is fixed by the densisty

! computing some quantities neeeded throughout the program
! these variables MUST NOT be touched elsewhere
n = nint( NOS**(1.0/3.0) )
side = ( NOS / rho )**(1.0/3.0) 
step = side / n

! initialisation of pointer variables: all pointers are nullified for safety reasons
NULLIFY( listH%Ptr, listT%Ptr, tmp%Ptr ) ! list is empty
FORALL( i = 1:NOS ) npoint( i )%Ptr => NULL() ! all pointers in the array are nullified

! time is set to zero before the start of the program
time = 0


! STEP 1: INITIALISATION

CALL init_cond( rx, ry, rz, vx, vy, vz ) ! first initialization of the system, no total energy / temperature assumed

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

! STEP 2: TERMALISATION
! the list should be destroyed after 10 time steps of simulation


build_list: DO i = 1, NOS
   flag = 0
   DO j = i+1,  NOS
      IF( distance( rx(i), ry(i), rz(i), rx(j), ry(j), rz(j) ) <= r_L ) THEN
         CALL add_node( listH%Ptr, listT%Ptr, j )
         IF( flag == 0 ) THEN ! it means that it is the first occurrence
            npoint( i )%Ptr => listT%Ptr
            flag = 1 
         ENDIF
      ENDIF
   ENDDO
   IF( flag == 0 ) NULLIFY( npoint( i )%Ptr ) ! in this case no particle j>i is at distance < r_L
ENDDO build_list

! compute the force at time: time

FORALL( i = 1 : NOS )
   fx( i ) = -F_c
   fy( i ) = -F_c
   fz( i ) = -F_c
   ax( i ) = 0
   ay( i ) = 0
   az( i ) = 0
ENDFORALL

tmp%Ptr => listH%Ptr  ! the temporary pointer is set at the beginning of the list


skim_array: DO i = 1, NOS
   
   IF( ASSOCIATED( npoint( i )%Ptr ) ) THEN ! in this case the i-th particle is interacting with another one
      ! skim the list untill next non-NULL npoint( j )%Ptr
      ! the last particle's interactions are accounted for in the previous elements of the list so that npoint( NOS )%Ptr => NULL() always
      skim_list: DO
         IF( (.NOT.ASSOCIATED(tmp%Ptr)) .OR. ASSOCIATED( tmp%Ptr, npoint( find_next( i, npoint ) )%Ptr ) )   EXIT
         ! it means that we have reached the next group of interacting particles, and the loop must be ended; if no next element can be found then the loop is eneded when tmp%Ptr reaches NULL()
         j = tmp%Ptr%data ! set the index of the second particle
         ! PRINT *, j
         r = distance( rx( i ), ry( i ), rx( i ), rx( j ), ry( j ), rz( j ) )
         fx( i ) = fx( i ) + force( rx( i ), rx( j ), r )
         fy( i ) = fy( i ) + force( ry( i ), ry( j ), r )
         fz( i ) = fz( i ) + force( rz( i ), rz( j ), r )
         ! we must also account for the tmp%Ptr%data particle which experiments an equal but opposite force
         fx( j ) = fx( j ) + force( rx( j ), rx( i ), r )
         fy( j ) = fy( j ) + force( ry( j ), ry( i ), r )
         fz( j ) = fz( j ) + force( rz( j ), rz( i ), r )
         ! go to the next element in the list
         tmp%Ptr => tmp%Ptr%nxtPtr
      ENDDO skim_list
   ENDIF
   ! otherwise jump to the next particle
ENDDO skim_array

! compute the acceleration at time: time
FORALL( i = 1 : NOS )
   ax( i ) = -fx( i )
   ay( i ) = -fy( i )
   az( i ) = -fz( i )
ENDFORALL

! VELOCITY-VERLET ALGORYTHM:
! 1. positions at time: time + delta_t
FORALL( i = 1 : NOS )
   rx( i ) = rx( i ) + vx( i )*deltat + (1./2.)*ax( i )*deltat**2
   ry( i ) = ry( i ) + vy( i )*deltat + (1./2.)*ay( i )*deltat**2
   rz( i ) = rz( i ) + vz( i )*deltat + (1./2.)*az( i )*deltat**2
ENDFORALL

! 2. velocities at atime: time + delta_t/2
FORALL( i = 1 : NOS )
   vx( i ) = vx( i ) + (1./2.)*ax( i )*deltat
   vy( i ) = vy( i ) + (1./2.)*ay( i )*deltat
   vz( i ) = vz( i ) + (1./2.)*az( i )*deltat
ENDFORALL

! 3. acceleration at time: time + delta_t
! 4. velocities at time: time + delta_t
! End of Velocity-Verlet algorythm



! CALL print_list( listH%Ptr, listT%Ptr )

! DO i= 1, NOS
!    IF( .NOT.ASSOCIATED( npoint(i)%Ptr) ) THEN
!       WRITE( *, * ) 0
!    ELSE
!       PRINT *, npoint( i )%Ptr%data
!    ENDIF
! ENDDO


CLOSE( unit = 1 )

ENDPROGRAM soft_spheres

!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************
!                                                                                       MAIN PROGRAM ENDS HERE!
!****************************************************************************************************************************************************************************************************************
!****************************************************************************************************************************************************************************************************************

