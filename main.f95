! Declaration of modules
MODULE prec_def
  ! Purpose: this module defines the precision of REAL type
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: sp = KIND( 1.0 ), &
       dp = SELECTED_REAL_KIND( 2 * PRECISION( 1.0_sp ) ), &
       qp = SELECTED_REAL_KIND( 2 * PRECISION( 1.0_dp ) ), &
       prec = dp!SELECTED_REAL_KIND( 6, 37 ) ! this is the chosen precision which
  ! will be used in the program
ENDMODULE prec_def

MODULE constants
  ! Purpose: declare all the constants used in the program,
  !          so that they are available to both the main
  !          and all the subroutines
  USE prec_def
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: NOS = 100, MAX_ITER = 500000, INTEG_STEP = 15000, EQ_STEP = 30000 ! NOS = number of spheres, number of iterations of the evolution step, number of steps needed for computing mean values, number of steps after which the system has reached equilibrium
  ! MAX_ITER = maximum number of iterations of the program
  REAL( prec ), PARAMETER :: e_D = -2.98_prec, T_D = 1.22_prec ! energy per particle to impose
  REAL( prec ), PARAMETER :: rho = 0.7_prec, deltat = 0.001_prec, r_c = 2.5_prec, r_L = 2.80_prec
  ! density, time interval, core distance, effective core distance
END MODULE constants

MODULE global

  ! Purpose: this module contains all the global variables and the general-aim
  !          functions and procedures.

  USE prec_def
  USE constants

  IMPLICIT NONE
  SAVE
  REAL( prec ) :: side, step, volume, & ! side of the cube, initial separation
       ! of the particles, and volume of the box
       F_c, u_c ! force and potential energy at r_c
  INTEGER :: n ! cube root of NOS rounded to the nearest integer

CONTAINS

  FUNCTION distance( x1, y1, z1, x2, y2, z2 )

    ! Purpose: it takes as arguments the coordinates of two particles and
    !          returns the distance between them, accounting for the boundary
    !          conditions

    IMPLICIT NONE
    REAL( prec ), INTENT( in ) :: x1, x2, y1, y2, z1, z2
    REAL( prec ) :: distance
    ! local var
    REAL( prec ) :: a1, a2, b1, b2, c1, c2
    a1 = x1; a2 = x2; b1 = y1; b2 = y2; c1 = z1; c2 = z2

    IF( ABS( a1-a2 ) > 0.5*side ) a2 = a2 - SIGN( 1._prec, a2 - a1 )*side
    IF( ABS( b1-b2 ) > 0.5*side ) b2 = b2 - SIGN( 1._prec, b2 - b1 )*side
    IF( ABS( c1-c2 ) > 0.5*side ) c2 = c2 - SIGN( 1._prec, c2 - c1 )*side

    distance = sqrt( (a2 - a1 )**2. + ( b2 - b1 )**2. + ( c2 - c1 )**2. )
    RETURN

  END FUNCTION distance

  FUNCTION comp( q1, q2 )

    ! Purpose: this function computes the difference between the same components of a vector

    IMPLICIT NONE
    REAL( prec ), INTENT( in ) :: q1, q2
    REAL( prec ) :: comp
    ! local var
    REAL( prec ) :: a1, a2, b1, b2, c1, c2
    a1 = q1; a2 = q2

    IF( ABS( a1-a2 ) > 0.5*side ) a2 = a2 - SIGN( 1._prec, a2 - a1 )*side
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
    TYPE( ptr2node ), DIMENSION( 1: NOS ), INTENT( IN ) :: array
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
    REAL( prec ), INTENT( in ) :: x1, x2, r
    REAL( prec ) :: force

    IF( r <= r_c ) THEN
       force = ( 24. * ( 2.*(1./r)**(13.) - (1./r)**(7.) ) + (-F_c) ) * comp( x2, x1 ) / r
    ELSE
       force = 0
    ENDIF
    RETURN

  ENDFUNCTION force

  FUNCTION pot_energy( r )

    USE prec_def
    USE constants

    IMPLICIT NONE
    REAL( prec ), INTENT( in ) :: r
    REAL( prec ) :: pot_energy

    IF( r <= r_c ) THEN
       pot_energy = 4. * ( (1./r)**12. - (1./r)**6. ) - u_c  - ( r - r_c ) * (-F_c)
    ELSE
       pot_energy = 0
    ENDIF
    RETURN
  ENDFUNCTION pot_energy

ENDMODULE global

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
     ! I define a pointer to the node type in order to define an
     ! "array of pointers"
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
    ! WRITE( *, * ) 'List is NOW empty.'

  ENDSUBROUTINE kill_list

  SUBROUTINE print_list( headPtr, tailPtr )
    ! Purpose: it prints the list
    !          it is not used in the program but can be useful for debug
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

!******************************************************************************

! THE MASSES, SIGMA AND EPSILON ARE ALL NORMALISED TO 1 SO THAT ALL DISTANCES
! ARE MEASURED IN UNITS OF SIGMA AND ALL ENERGIES ARE MEASURED IN UNITS OF
! EPSILON

!******************************************************************************
!****************************** START OF MAIN PROGRAM *************************
!******************************************************************************
!******************************************************************************

PROGRAM soft_spheres

  USE prec_def
  USE constants
  USE dynamic_list
  USE global

  ! VARIABLE DECLARATION
  IMPLICIT NONE
  REAL( prec ) :: time, &
       r ! r is used to store the distance between two particles
  REAL( prec ) :: kin, pot, mec ! total kinetic, potential
  ! and mechanical energy resp.
  REAL( prec ), dimension( 1:NOS ) :: rx, ry, rz,&
       rx0, ry0, rz0,&
       vx, vy, vz,&
       fx, fy, fz,&
       ax, ay, az,&
       u
  REAL( prec ) :: vxc, vyc, vzc, & ! CoM velocities
       fijx, fijy, fijz
  
  ! Variables needed for computing the temperaure
  REAL( prec ) :: kin_sum, & ! sum of the kinetic energy during the time sample
       integ_time, & ! time sample to compute the mean kinetic energy
       temp, temp_old, & ! temperature
       msqdisp = 0 ! mean square displacement
  ! Variables needed for computing the pressure
  REAL( prec ) :: w, & ! 1/dim * <sum ri*Fi >
       press ! pressure / (rho*temp) - 1 = w/(NOS * temp)

  INTEGER :: i, j, iter, integ_count = 0
  INTEGER :: flag = 0, ok = 0, msqd_first = 0 ! put ok = 1 if the temperature is enforced

  TYPE( ptr2node ) :: listH, listT  ! these variables contain the pointers to
  ! the head and the tail of the "list", see notes
  TYPE( ptr2node ), DIMENSION( 1 : NOS ) :: npoint ! same as in the notes but
  ! instead of containing the index position inside the list it points to the
  ! first occurence
  TYPE( ptr2node ) :: tmp ! this pointer is needed to skim the list

  ! variables for generating random values
  REAL( prec ) :: randnum 
  INTEGER, DIMENSION( : ), ALLOCATABLE :: seed
  INTEGER, DIMENSION( 1 : 8 ) :: dt_seed
  INTEGER :: n_seed

  ! these variables are temporary and sould be eliminated

  CHARACTER( LEN = 70 ) :: fn
  INTEGER, PARAMETER :: numfiles = 10
  INTEGER, PARAMETER :: outunit = 44
  INTEGER :: filenum
  
  !****************************************************************************
  !****************************** END OF PREAMBLE *****************************
  !****************************************************************************

  OPEN( UNIT = 1, FILE = "posit.dat", STATUS = "replace", &
       ACCESS = "sequential", position = "rewind" )
  OPEN( UNIT = 2, FILE = "energy.dat", STATUS = "replace", &
       ACCESS = "sequential", POSITION = "rewind" )
  OPEN( UNIT = 3, FILE = "obs.dat", STATUS = "replace", &
       ACCESS = "sequential", POSITION = "rewind" )

  OPEN( UNIT = 4, FILE = "speed.dat", STATUS = "replace", &
       ACCESS = "sequential", POSITION = "rewind" )

  OPEN( UNIT = 12, FILE = "msqd.dat", STATUS = "REPLACE", ACCESS = "SEQUENTIAL", POSITION = "REWIND" )
  
  ! computing some quantities neeeded throughout the program
  ! these variables MUST NOT be modified elsewhere
  n = CEILING( NOS**(1.0/3.0) ) 
  side = ( NOS / rho )**(1.0/3.0) 
  step = side / n
  volume = REAL( NOS ) / rho ! the volume is fixed by the density
  F_c = 24. * ( 2. * ( 1./(r_c)**13. ) - 1./(r_c)**7. )
  u_c = 4. * ( 1./(r_c)**12. -1./(r_c)**6. )

  !****************************************************************************
  
  ! initialisation of pointer variables:
  ! all pointers are nullified for safety reasons
  NULLIFY( listH%Ptr, listT%Ptr, tmp%Ptr ) ! list is empty
  FORALL( i = 1 : NOS ) npoint( i )%Ptr => NULL() ! all pointers in the array
  ! are nullified

  ! time is set to zero before the start of the program
  time = 0
  
  integ_time = 0
  kin_sum = 0
  temp = 0; temp_old = 0
  w = 0
  press = 0

  ! PRINT *, n, side

  !****************************************************************************

  ! STEP 1: INITIALISATION
  ! assigning positions
  DO i = 1, NOS
     rx( i ) = step * ( 1/2.0 + MOD( i - 1 , n ) )
     ry( i ) = step * ( 1/2.0 + MOD( ( i - 1 ) / n , n ) )
     rz( i ) = step * ( 1/2.0 + ( i - 1 ) / (n*n) )
  ENDDO

  ! setting seed for speen generator
  CALL RANDOM_SEED( size = n_seed )
  ALLOCATE( seed( 1 : n_seed ) )
  CALL RANDOM_SEED( get = seed )
  CALL DATE_AND_TIME( values = dt_seed )
  seed( n_seed ) = dt_seed( 8 )
  seed( 1 ) = dt_seed( 8 ) * dt_seed( 7 ) * dt_seed( 6 )
  CALL RANDOM_SEED( put = seed )
  DEALLOCATE( seed )
  ! done setting seed

  ! assigning random velocities
  DO i= 1, NOS
     CALL RANDOM_NUMBER( randnum )
     vx( i ) = 2. * randnum - 1.
  ENDDO

  DO i= 1, NOS
     CALL RANDOM_NUMBER( randnum )
     vy( i ) = 2. * randnum - 1.
  ENDDO

  DO i= 1, NOS
     CALL RANDOM_NUMBER( randnum )
     vz( i ) = 2. * randnum - 1.
  ENDDO

  ! computing the cm velocitiy
  vxc = 0; vyc = 0; vzc = 0
  DO i = 1, NOS 
     vxc = vxc + vx( i ) / NOS
     vyc = vyc + vy( i ) / NOS
     vzc = vzc + vz( i ) / NOS
  ENDDO

  ! subtracting the cm velocity
  FORALL( i = 1 : NOS )
     vx( i ) = vx( i ) - vxc
     vy( i ) = vy( i ) - vyc
     vz( i ) = vz( i ) - vzc
  ENDFORALL

  ! verifying cm velocity = 0
  vxc = 0; vyc = 0; vzc = 0
  DO i = 1, NOS 
     vxc = vxc + vx( i )
     vyc = vyc + vy( i )
     vzc = vzc + vz( i )
  ENDDO

  DO i = 1, NOS
     WRITE( unit = 1, fmt = 1 ) rx( i ), ry( i ) , rz( i )
1    FORMAT( f20.15, f20.15, f20.15 )
  ENDDO

  !****************************************************************************

  ! STEP 2: TERMALISATION
  ! the list should be destroyed after 10 time steps of simulation

 ! measure: DO filenum = 1, 000numfiles
     !WRITE( fn, FMT = '(i0,a)' ) filenum, '.dat'
  !OPEN( UNIT = outunit, file = fn, form = 'formatted' )

  ! First step
  !Build list for the first time
  
  build_list_zero: DO i = 1, NOS
     flag = 0
     DO j = i+1, NOS
        IF( distance( rx(i), ry(i), rz(i), rx(j), ry(j), rz(j) ) < r_L )&
             THEN
           CALL add_node( listH%Ptr, listT%Ptr, j )
           IF( flag == 0 ) THEN ! it means it is the first occurrence
              npoint( i )%Ptr => listT%Ptr
              flag = 1 
           ENDIF
        ENDIF
     ENDDO
     IF( flag == 0 ) NULLIFY( npoint( i )%Ptr ) ! in this case no
     ! particle j>i is at distance < r_L
  ENDDO build_list_zero

  ! computing potential energy and forces at time zero
  ! compute the force at time: time = 0
  FORALL( i = 1 : NOS )
     fx( i ) = 0.
     fy( i ) = 0.
     fz( i ) = 0.
     u( i ) = 0.
  ENDFORALL

  tmp%Ptr => listH%Ptr  ! the temporary pointer is set at the beginning
  ! of the list

  skim_array_zero: DO i = 1, NOS
     IF( ASSOCIATED( npoint( i )%Ptr ) ) THEN ! in this case the i-th
        ! particle is interacting with another one

        ! skim the list untill next non-NULL npoint( j )%Ptr
        ! the last particle's interactions are accounted for in the
        ! previous elements of the list so that npoint( NOS )%Ptr =>
        ! NULL() always
        skim_list_zero: DO
           IF( (.NOT.ASSOCIATED(tmp%Ptr)) .OR. ASSOCIATED( tmp%Ptr, &
                npoint( find_next( i, npoint ) )%Ptr ) )   EXIT
           ! it means that we have reached the next group of interacting
           ! particles, and the loop must be ended; if no next element
           ! can be found then the loop is eneded when tmp%Ptr reaches
           ! NULL()
           j = tmp%Ptr%data ! set the index of the second particle
           ! PRINT *, j
           r = distance( rx( i ), ry( i ), rz( i ), rx( j ), ry( j ), &
                rz( j ) )
           ! PRINT *, i, j, r

           fx( i ) = fx( i ) + force( rx( i ), rx( j ), r )
           fy( i ) = fy( i ) + force( ry( i ), ry( j ), r )
           fz( i ) = fz( i ) + force( rz( i ), rz( j ), r )
           u( i ) = u( i ) + pot_energy( r )
           ! we must also account for the tmp%Ptr%data particle which
           ! experiments an equal but opposite force
           fx( j ) = fx( j ) - force( rx( i ), rx( j ), r )
           fy( j ) = fy( j ) - force( ry( i ), ry( j ), r )
           fz( j ) = fz( j ) - force( rz( i ), rz( j ), r )
           ! u( j ) = u( j ) + pot_energy( r )
           ! go to the next element in the list

           ! PRINT *, 'Forces at time:', time
           ! DO k = 1, NOS
           !    PRINT *, fx(k), fy(k), fz(k)
           ! ENDDO

           tmp%Ptr => tmp%Ptr%nxtPtr
        ENDDO skim_list_zero
     ENDIF
     ! otherwise jump to the next particle
  ENDDO skim_array_zero

  ! compute the energy
  kin = (1./2.) * ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) &
       + DOT_PRODUCT( vz, vz ) )
  pot = SUM( u( 1: NOS ) )
  mec = kin + pot

  PRINT *, kin, pot/NOS, mec, e_D

  evolution: DO iter = 0, MAX_ITER 

     IF( MOD( iter, 100 ) == 0 ) PRINT *, "Completed:", &
          (100 * iter) / MAX_ITER, "%"

     refresh_list: IF( MOD( iter, 10 ) == 0 ) THEN ! refreshing list every 10
        ! time steps (as suggested in the notes)

        CALL kill_list( listH%Ptr, listT%Ptr )

        build_list: DO i = 1, NOS
           flag = 0
           DO j = i+1, NOS
              IF( distance( rx(i), ry(i), rz(i), rx(j), ry(j), rz(j) ) < r_L )&
                   THEN
                 CALL add_node( listH%Ptr, listT%Ptr, j )
                 IF( flag == 0 ) THEN ! it means it is the first occurrence
                    npoint( i )%Ptr => listT%Ptr
                    flag = 1 
                 ENDIF
              ENDIF
           ENDDO
           IF( flag == 0 ) NULLIFY( npoint( i )%Ptr ) ! in this case no
           ! particle j>i is at distance < r_L
        ENDDO build_list

     ENDIF refresh_list

     ! time_zero: IF( time == 0 ) THEN ! compute force, potential energy,
     !    ! acceleration at time = 0

     !    ! compute the force at time: time = 0
     !    FORALL( i = 1 : NOS )
     !       fx( i ) = 0.
     !       fy( i ) = 0.
     !       fz( i ) = 0.
     !       u( i ) = 0.
     !    ENDFORALL

     !    tmp%Ptr => listH%Ptr  ! the temporary pointer is set at the beginning
     !    ! of the list

     !    skim_array_zero: DO i = 1, NOS
     !       IF( ASSOCIATED( npoint( i )%Ptr ) ) THEN ! in this case the i-th
     !          ! particle is interacting with another one
              
     !          ! skim the list untill next non-NULL npoint( j )%Ptr
     !          ! the last particle's interactions are accounted for in the
     !          ! previous elements of the list so that npoint( NOS )%Ptr =>
     !          ! NULL() always
     !          skim_list_zero: DO
     !             IF( (.NOT.ASSOCIATED(tmp%Ptr)) .OR. ASSOCIATED( tmp%Ptr, &
     !                  npoint( find_next( i, npoint ) )%Ptr ) )   EXIT
     !             ! it means that we have reached the next group of interacting
     !             ! particles, and the loop must be ended; if no next element
     !             ! can be found then the loop is eneded when tmp%Ptr reaches
     !             ! NULL()
     !             j = tmp%Ptr%data ! set the index of the second particle
     !             ! PRINT *, j
     !             r = distance( rx( i ), ry( i ), rz( i ), rx( j ), ry( j ), &
     !                  rz( j ) )
     !             ! PRINT *, i, j, r

     !             fx( i ) = fx( i ) + force( rx( i ), rx( j ), r )
     !             fy( i ) = fy( i ) + force( ry( i ), ry( j ), r )
     !             fz( i ) = fz( i ) + force( rz( i ), rz( j ), r )
     !             u( i ) = u( i ) + pot_energy( r )
     !             ! we must also account for the tmp%Ptr%data particle which
     !             ! experiments an equal but opposite force
     !             fx( j ) = fx( j ) - force( rx( i ), rx( j ), r )
     !             fy( j ) = fy( j ) - force( ry( i ), ry( j ), r )
     !             fz( j ) = fz( j ) - force( rz( i ), rz( j ), r )
     !             ! u( j ) = u( j ) + pot_energy( r )
     !             ! go to the next element in the list

     !             ! PRINT *, 'Forces at time:', time
     !             ! DO k = 1, NOS
     !             !    PRINT *, fx(k), fy(k), fz(k)
     !             ! ENDDO

     !             tmp%Ptr => tmp%Ptr%nxtPtr
     !          ENDDO skim_list_zero
     !       ENDIF
     !       ! otherwise jump to the next particle
     !    ENDDO skim_array_zero

     !    ! compute the energy
     !    kin = (1./2.) * ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) &
     !         + DOT_PRODUCT( vz, vz ) )
     !    pot = SUM( u( 1: NOS ) )
     !    mec = kin + pot

     !    PRINT *, kin, pot/NOS, mec, e_D

     ! ENDIF time_zero

     ! SCALE VELOCITIES FOR A GIVEN VALUE OF ENERGY PER PARTICLE
     ! uncomment the ! following if the energy per particle is fixed
     set_energy: IF( e_D - pot/NOS >= 0 .AND. ok == 0 ) THEN
        ! the system must reach a sufficient potential energy per particle
        ! before setting the total energy to the desired value
        FORALL( i = 1 : NOS)
           vx( i ) = vx( i ) * sqrt( ( e_D - pot/NOS ) / ( kin/NOS ) )
           vy( i ) = vy( i ) * sqrt( ( e_D - pot/NOS ) / ( kin/NOS ) )
           vz( i ) = vz( i ) * sqrt( ( e_D - pot/NOS ) / ( kin/NOS ) )
           ! PRINT *, vx(i), vy(i), vz(i)
        ENDFORALL
        PRINT *, 'Energy has been succefully set to the desired value!', time
        ok = 1 ! switch flag for energy setting
        ! compute the energy after rescaling
        kin = (1./2.) * ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) &
             + DOT_PRODUCT( vz, vz ) )
        pot = SUM( u( 1 : NOS ) )
        mec = kin + pot
        PRINT *, kin, pot, mec, mec/NOS
     ENDIF set_energy

     ! block for the computation of observables
     IF( iter >= EQ_STEP ) THEN
        kin_sum = kin_sum + kin * deltat
        integ_time = integ_time + deltat
        ! temp_old = temp

        !DO i = 1, NOS
         !  w = w + 1./3. * ( comp( 0., rx( i ) ) * fx( i ) + comp( 0., ry( i ) ) * fy (i ) + &
         !       comp( 0., rz( i ) ) * fz( i ) ) * deltat
        !ENDDO
        ! w = w + (1./3.) * ( DOT_PRODUCT( comp( rx, 0. ), fx ) + DOT_PRODUCT( comp( ry, 0. ) , fy ) + DOT_PRODUCT( comp( rz, 0 ) , fz ) ) * deltat

        IF( msqd_first == 0 ) THEN ! store initial positions for the computation of the mean square displacement
           FORALL( i = 1 : NOS )
              rx0( i ) = rx( i )
              ry0( i ) = ry( i )
              rz0( i ) = rz( i )
           ENDFORALL
           msqd_first = 1
        ENDIF
        msqdisp = 0
        DO i = 1, NOS
           msqdisp = msqdisp + distance( rx( i ), ry( i ), rz( i ), rx0( i ), ry0( i ), rz0( i ) )**2 
        ENDDO
        msqdisp = msqdisp / NOS
        WRITE( UNIT = 12, FMT = '(2ES12.4E2)' ) time, msqdisp
        
        IF( integ_count == INTEG_STEP ) THEN
           temp = kin_sum / ( 2._prec * integ_time * NOS )
           press = w/( 3.0_prec * NOS * temp * integ_time )
           
           ! SCALE VELOCITIES FOR A GIVEN VALUE OF TEMPERATURE
           ! ! uncomment the following if the temperature is fixed
           ! IF( iter <= EQ_STEP ) THEN
           !    FORALL( i = 1 : NOS )
           !       vx( i ) = vx( i ) * SQRT( T_D / temp )
           !       vy( i ) = vy( i ) * SQRT( T_D / temp )
           !       vz( i ) = vz( i ) * SQRT( T_D / temp )
           !    ENDFORALL
           ! ENDIF                

           WRITE( UNIT = 3, FMT = '(2F9.3, 2ES12.4E2 )' ) time, integ_time, temp, press
           integ_time = 0.
           kin_sum = 0.
           w = 0.
           integ_count = 0
        ENDIF
        integ_count = integ_count + 1
        ! WRITE(*,*) time_sample, temp
     ENDIF

     ! print "time\t energy" on file 
     WRITE( UNIT = 2, FMT = '(F9.3, 4ES12.4E2 )' ) time, kin/NOS, pot/NOS, mec/NOS

     ! PRINT *, 't = ',time
     ! DO i = 1, NOS
     !    PRINT *, 'particle: ', i
     !    PRINT *, 'positions: ', rx( i ), ry( i ), rz( i )
     !    PRINT *, 'velocities: ', vx( i ), vy( i ), vz( i )
     !    PRINT *, 'accelerations: ', fx( i ), fy( i ), fz( i )
     ! ENDDO


     time = time + deltat
     kin = 0; pot = 0; mec = 0;

     ! VELOCITY-VERLET ALGORYTHM:

     ! 1. positions at time: time + deltat
     ! here the accelaration is evaluated at time: time (coming from the time_zero IF or from the previous execution of the evolution DO)

     FORALL( i = 1 : NOS )
        rx( i ) = rx( i ) + vx( i )*deltat + (1./2.)*fx( i )*deltat**2.
        ry( i ) = ry( i ) + vy( i )*deltat + (1./2.)*fy( i )*deltat**2.
        rz( i ) = rz( i ) + vz( i )*deltat + (1./2.)*fz( i )*deltat**2.
     ENDFORALL

     BC: DO i = 1, NOS 
        IF( rx( i ) > side .OR. rx( i ) < 0 ) rx( i ) = rx( i ) - FLOOR( rx( i ) / side ) * side
        IF( ry( i ) > side .OR. ry( i ) < 0 ) ry( i ) = ry( i ) - FLOOR( ry( i ) / side ) * side
        IF( rz( i ) > side .OR. rz( i ) < 0 ) rz( i ) = rz( i ) - FLOOR( rz( i ) / side ) * side
     ENDDO BC


     DO i = 1, NOS
        IF( rx( i ) > side .OR. rx( i ) < 0 ) CALL ABORT
        IF( ry( i ) > side .OR. ry( i ) < 0 ) CALL ABORT
        IF( rz( i ) > side .OR. rz( i ) < 0 ) CALL ABORT
     ENDDO

     ! 2. velocities at time: time + deltat/2

     FORALL( i = 1 : NOS )
        vx( i ) = vx( i ) + (1./2.) * fx( i ) * deltat
        vy( i ) = vy( i ) + (1./2.) * fy( i ) * deltat
        vz( i ) = vz( i ) + (1./2.) * fz( i ) * deltat
     ENDFORALL

     ! 3. acceleration at time: time + deltat
     ! compute the force at time: time + deltat

     FORALL( i = 1 : NOS )
        fx( i ) = 0.
        fy( i ) = 0.
        fz( i ) = 0.
        u( i ) = 0.
     ENDFORALL

     tmp%Ptr => listH%Ptr  ! the temporary pointer is set at the beginning
     ! of the list

     skim_array: DO i = 1, NOS
        IF( ASSOCIATED( npoint( i )%Ptr ) ) THEN ! in this case the i-th
           ! particle is interacting with another one

           ! skim the list untill next non-NULL npoint( j )%Ptr
           ! the last particle's interactions are accounted for in the previous
           ! elements of the list so that npoint( NOS )%Ptr => NULL() always
           skim_list: DO
              IF( (.NOT.ASSOCIATED(tmp%Ptr) .OR. ASSOCIATED( tmp%Ptr, &
                   npoint( find_next( i, npoint ) )%Ptr ) ) )   EXIT
              ! it means that we have reached the next group of interacting
              ! particles, and the loop must be ended; if no next element can
              ! be found then the loop is eneded when tmp%Ptr reaches NULL()

              j = tmp%Ptr%data ! set the index of the second particle

              r = distance( rx( i ), ry( i ), rz( i ), rx( j ), ry( j ), &
                   rz( j ) )
              ! PRINT *, i, j, r, side 
              fijx = force( rx( i ), rx( j ), r )
              fijy = force( ry( i ), ry( j ), r )
              fijz = force( rz( i ), rz( j ), r )

              fx( i ) = fx( i ) + fijx
              fy( i ) = fy( i ) + fijy
              fz( i ) = fz( i ) + fijz
              IF( integ_count == INTEG_STEP ) THEN
                 w = w &
                      + comp( rx( i ) , rx( j ) ) * fijx &
                      + comp( ry( i ) , ry( j ) ) * fijy &
                      + comp( rz( i ) , rz( j ) ) * fijz
              ENDIF
              
              ! fx( i ) = fx( i ) + force( rx( i ), rx( j ), r )
              ! fy( i ) = fy( i ) + force( ry( i ), ry( j ), r )
              ! fz( i ) = fz( i ) + force( rz( i ), rz( j ), r )
              u( i ) = u( i ) + pot_energy( r )
              ! we must also account for the tmp%Ptr%data particle which
              ! experiments an equal but opposite force
              fx( j ) = fx( j ) - fijx
              fy( j ) = fy( j ) - fijy
              fz( j ) = fz( j ) - fijz
              
              ! fx( j ) = fx( j ) - force( rx( i ), rx( j ), r )
              ! fy( j ) = fy( j ) - force( ry( i ), ry( j ), r )
              ! fz( j ) = fz( j ) - force( rz( i ), rz( j ), r )
              ! u( j ) = u( j ) + pot_energy( r )
              ! go to the next element in the list
              tmp%Ptr => tmp%Ptr%nxtPtr
           ENDDO skim_list
        ENDIF
        ! otherwise jump to the next particle
     ENDDO skim_array

     ! 4. velocities at time: time + deltat
     ! here the acceleration is evaluated at time:  time + deltat

     FORALL( i = 1 : NOS )
        vx( i ) = vx( i ) + (1./2.) * fx( i ) * deltat
        vy( i ) = vy( i ) + (1./2.) * fy( i ) * deltat
        vz( i ) = vz( i ) + (1./2.) * fz( i ) * deltat
     ENDFORALL

     ! End of VELOCITY-VERLET ALGORYTHM

     ! compute the energy
     kin = (1._prec/2._prec) * ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) + &
          DOT_PRODUCT( vz, vz ) )
     pot = SUM( u( 1: NOS ))
     mec = kin + pot
     ! IF( iter > 1000 .AND. &
     !      ABS( temp_old - temp ) < EPSILON( temp) ) THEN
     !    DO i = 1, NOS
     !       WRITE( UNIT = 4, FMT = '(F7.5)' ) SQRT( vx( i )**2. +  vy( i )**2.&
     !            + vz( i )**2. )
     !    ENDDO
     !    PRINT *, 'Velecities have been written!'
     !    flag_term = 1
     ! ENDIF
     
  ENDDO evolution

  CALL kill_list( listH%Ptr, listT%Ptr )



  !CLOSE( outunit )
!ENDDO measure

  CLOSE( UNIT = 1 )
  CLOSE( UNIT = 2 )
  CLOSE( UNIT = 3 )

  CLOSE( UNIT = 4 )
  CLOSE( UNIT = 12 )
ENDPROGRAM soft_spheres

!******************************************************************************
!******************************************************************************
!**************************** MAIN PROGRAM ENDS HERE! *************************
!******************************************************************************
!******************************************************************************
