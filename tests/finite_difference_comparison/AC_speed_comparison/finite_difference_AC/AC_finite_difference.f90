MODULE parameters
IMPLICIT NONE
SAVE

	CHARACTER(LEN=3), PARAMETER :: 	run_num = '001'
	CHARACTER(LEN=0), PARAMETER :: 	output_path = ''
	
	!! Grid size
	
	
	REAL(KIND=8), PARAMETER :: 	domain_size_x = 100.0d0 				!! nm
	REAL(KIND=8), PARAMETER :: 	domain_size_y = 100.0d0 				!! nm
	REAL(KIND=8), PARAMETER :: 	domain_size_z = 100.0d0 					!! nm
	
	INTEGER, PARAMETER :: 		xg = 128/2 !128 !224 !128 !184 !127 !181 !241 !181 !121
	INTEGER, PARAMETER :: 		yg = 128/2 !128 !224 !128 !256 !128 !184 !127 !181 !241 !181 !121
	INTEGER, PARAMETER :: 		zg = 128/2 !92 !92 !163 !91 !195 !195 !121 !91 !61 
	
	INTEGER, PARAMETER :: 		bb = 1
	
	REAL(KIND=8), PARAMETER :: 	dz = domain_size_z/DBLE(zg) 						!! Grid spacing (nm)
	REAL(KIND=8), PARAMETER :: 	dx = domain_size_x/DBLE(xg) 							!! Grid spacing (nm)
	REAL(KIND=8), PARAMETER :: 	dy = domain_size_y/DBLE(yg)					 							!! Grid spacing (nm)
	
	!! The number of subdomains in each direction for MPI
	INTEGER, PARAMETER ::		num_subdomain_x = 4 !8 !8 !8 !8 !8 !4 !8 !6
	INTEGER, PARAMETER ::		num_subdomain_y = 2 !4 !4 !4 !4 !4 !1 !1
	INTEGER, PARAMETER ::		num_subdomain_z = 2 !4 !2 !2 !3 !6 !4
	
	!! Particle initialization
	REAL(KIND=8), PARAMETER ::	particle_x = domain_size_x/2.0d0 !75.0d0									!! nm
	REAL(KIND=8), PARAMETER :: 	particle_y = domain_size_y/2.0d0 !75.0d0									!! nm
	REAL(KIND=8), PARAMETER ::	particle_z = domain_size_z/2.0d0						!! nm
	REAL(KIND=8), PARAMETER :: 	rad = 40.0d0									!! nm
	
	
	!! Phase field constants
	
	REAL(KIND=8), PARAMETER :: 	A_alpha = 24.7939d0
	REAL(KIND=8), PARAMETER :: 	B_alpha = -1.6752d0
	REAL(KIND=8), PARAMETER :: 	C_alpha = 1.9453d-6
	REAL(KIND=8), PARAMETER :: 	A_beta = 37.9316d0
	REAL(KIND=8), PARAMETER :: 	B_beta = -10.7373d0
	REAL(KIND=8), PARAMETER :: 	C_beta = 0.5401d0
	
	REAL(KIND=8), PARAMETER :: 	L = 1.0d0
	REAL(KIND=8), PARAMETER :: 	kappa = 4.0d0
	
	REAL(KIND=8), PARAMETER :: 	delta = 2.0d0 !0.215d0
	
	REAL(KIND=8), PARAMETER :: 	pi = ACOS(-1.0d0)
	
	INTEGER, PARAMETER :: 		max_iter = 30000
	REAL(KIND=8), PARAMETER :: 	dt = 1.0d-3											!! seconds
	
	!! Tags for MPI communication
	INTEGER, PARAMETER :: c1Tag=91,c2Tag=92,c3Tag=93,c4Tag=94, pTag=95
	INTEGER, PARAMETER :: r1Tag=101,l1Tag=102,f1Tag=103,b1Tag=104,t1Tag=105,d1Tag=106    
	INTEGER, PARAMETER :: r2Tag=111,l2Tag=112,f2Tag=113,b2Tag=114,t2Tag=115,d2Tag=116    
	INTEGER, PARAMETER :: Sr1Req=81,Sl1Req=82,Sf1Req=83,Sb1Req=84,St1Req=85,Sd1Req=86
	INTEGER, PARAMETER :: Rr1Req=71,Rl1Req=72,Rf1Req=73,Rb1Req=74,Rt1Req=75,Rd1Req=76
	INTEGER, PARAMETER :: Sr2Req=61,Sl2Req=62,Sf2Req=63,Sb2Req=64,St2Req=65,Sd2Req=66
	INTEGER, PARAMETER :: Rr2Req=51,Rl2Req=52,Rf2Req=53,Rb2Req=54,Rt2Req=55,Rd2Req=56
	
END MODULE parameters

!! CHAC_finite_difference: This program solves the coupled Cahn-Hilliard-Allen-Cahn system
!! of equations. It was developed for comparison to the PRISMS-PF finite element code to 
!! ensure that its performance was competitive with finite difference.
!!
!! Author: Stephen DeWitt (stvdwtt@umich.edu)
!! Date Started: 1/6/2016
!!
!! =============================================================================

PROGRAM CHAC_finite_difference
USE parameters
IMPLICIT NONE

INCLUDE 'mpif.h'
INTEGER :: len
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname

INTEGER errcode
INTEGER rank, nsize
INTEGER :: dims(1:3),coord(1:3),left,right,front,back,top,bot,ndims,cartcomm
LOGICAL :: reorder,periodic(1:3)
INTEGER L1, L2, L3, S1, S2, S3, Id1, Id2, Id3, Rm1, Rm2, Rm3

CALL outputDomainInfo()

!! Starting the MPI section of the code

CALL MPI_INIT (errcode)

CALL MPI_COMM_RANK (MPI_COMM_WORLD, rank, errcode)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nsize, errcode)
CALL MPI_GET_PROCESSOR_NAME(hostname, len, errcode)

IF (rank .eq. 0) PRINT *, 'number of partitions', num_subdomain_x*num_subdomain_y*num_subdomain_z

!! Start off by partitioning the domain

!! Number of dimensions (for the MPI Cartesian communicator)
ndims = 3

!! Set the number of partitions in each direction
dims(1) = num_subdomain_x	
dims(2) = num_subdomain_y
dims(3) = num_subdomain_z

periodic(1) = .TRUE.
periodic(2) = .TRUE.
periodic(3) = .FALSE.
reorder = .true.

!! Create the Cartesian communicator
call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periodic,reorder,cartcomm,errcode)
!PRINT *, 'rank', rank, 'cartcomm', cartcomm

!! Calculate the non-ghost-zone length of each partition (L#) and the remaining grid points (Rm#)
L1 = floor(xg/real(num_subdomain_x))                           
Rm1 = mod(xg,num_subdomain_x)                                  

L2 = floor(yg/real(num_subdomain_y))
Rm2 = mod(yg,num_subdomain_y)

L3 = floor(zg/real(num_subdomain_z))
Rm3 = mod(zg,num_subdomain_z)

!! Get coordinates of each partition using MPI_CART_COORDS
call MPI_CART_COORDS(cartcomm,rank,ndims,coord,errcode)

!! The index of the partition in all three directions
Id3 = coord(3)                       
Id2 = coord(2)
Id1 = coord(1)

!! Defines the 6 neighbor partitions
call MPI_CART_SHIFT(cartcomm,0,1,left,right,errcode) !+/- x shifts are left/right
call MPI_CART_SHIFT(cartcomm,1,1,front,back,errcode) !+/- y shifts are front/back
call MPI_CART_SHIFT(cartcomm,2,1,bot,top,errcode) !+/- z shifts are top/bottom

!! Find the starting point for each partition with respect to the global array location (S#), also find the length of each partition (L#)
if(Id1 == 0)then                                          ! starting point of domain is set to 1, end point is n1
   S1 = 1-bb
   L1=L1+2*bb
elseif(Id1 == num_subdomain_x-1)then   
   S1 = Id1*L1+1-bb
   L1 = L1+2*bb+Rm1
else
   S1 = Id1*L1+1-bb
   L1= L1+2*bb
endif
if(Id2 == 0)then                                          ! starting point of domain is set to 1, end point is n1
   S2 = 1-bb
   L2= L2+2*bb
elseif(Id2 == num_subdomain_y-1)then   
   S2 = Id2*L2+1-bb
   L2 = L2+2*bb+Rm2
else
   S2 = Id2*L2+1-bb
   L2 = L2+2*bb
endif
if(Id3 == 0)then                                          ! starting point of domain is set to 1, end point is n1
   S3 = 1-bb
   L3= L3+2*bb
elseif(Id3 == num_subdomain_z-1)then   
   S3 = Id3*L3+1-bb
   L3 = L3+2*bb+Rm3
else
   S3 = Id3*L3+1-bb
   L3= L3+2*bb
endif

CALL main_subroutine(Id1,Id2,Id3,S1,S2,S3,L1,L2,L3,left,right,front,back,top,bot,cartcomm)	
					
CALL MPI_FINALIZE(errcode)


END PROGRAM

!! ===================================================================================
!! Main subroutine where most of the calculations are done
!! ===================================================================================

SUBROUTINE main_subroutine(Id1,Id2,Id3,S1,S2,S3,L1,L2,L3,left,right,front,back,top,bot,cartcomm)
USE parameters
IMPLICIT NONE

INCLUDE 'mpif.h'
INTEGER :: L1, L2, L3, S1, S2, S3, Id1, Id2, Id3, Is1, Is2, Is3, Ie1, Ie2, Ie3, k, iter
INTEGER :: left,right,front,back,top,bot, cartcomm 	
REAL(KIND=8), DIMENSION(S1:S1+L1-1,S2:S2+L2-1,S3:S3+L3-1) :: eta, x, y, z
REAL(KIND=8) :: time_elapsed, dist_particle
INTEGER :: i,j, out_num, output_crit
INTEGER ::  errcode 
CHARACTER(len=3) :: counter

out_num = 0

time_elapsed = 0.0d0 

!! Real grid points inside the subdomain (without ghost cells) 
Is1 = S1 + bb
Ie1 = S1 + L1 - 1 - bb
Is2 = S2 + bb
Ie2 = S2 + L2 - 1 - bb
Is3 = S3 + bb
Ie3 = S3 + L3 - 1 - bb

!! Initialize the fields

DO k = Is3-1,Ie3+1
DO j = Is2-1,Ie2+1
DO i = Is1-1,Ie1+1

	x(i,j,k) = dx * DBLE(i) -dx/2.0d0 
	y(i,j,k) = dy * DBLE(j) -dy/2.0d0
	z(i,j,k) = dz * DBLE(k)	-dz/2.0d0 

	!! Distance functions for spherical particles
	!dist_particle = SQRT(DBLE((particle_1_x-x(i,j,k))**2)+DBLE((particle_1_y-y(i,j,k))**2)+DBLE((particle_1_z-z(i,j,k))**2))
	
	!! Distance functions for ellipsoidal particles (squished in y direction)
	dist_particle = SQRT(DBLE((particle_x-x(i,j,k))**2) + DBLE(((particle_y-y(i,j,k))*2.0d0)**2) &
								+ DBLE( (particle_z-z(i,j,k))**2) )
							
	!! Flat interface							
	!dist_particle = SQRT((particle_z-z(i,j,k))**2)

	!! 3D Sphere
	eta(i,j,k) = 0.5d0*(1.0d0 - ( TANH( (dist_particle - rad)/delta )))
		
END DO
END DO
END DO

!! Apply boundary conditions

call commuBC(Id1,Id2,Id3,L1,L2,L3,S1,S2,S3,eta,left,right,front,back,top,bot,cartcomm)  	  
CALL external_BCs_Zero_Derivative(eta,.TRUE.,1,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,1,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.TRUE.,2,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,2,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.TRUE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)

!! Initial Outputs
IF (Id1+Id2+Id3.eq. 0) PRINT *, 'Time step: ', dt
CALL writeLotsFiles(eta,'eta_r',out_num,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)
CALL writeLotsFiles(x,'xgrid',out_num,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)
CALL writeLotsFiles(y,'ygrid',out_num,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)
CALL writeLotsFiles(z,'zgrid',out_num,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)


out_num = out_num + 1


!! -------------------------------------------------------------------------
!! Iterate through time
!! -------------------------------------------------------------------------

DO iter = 1,max_iter

	CALL allenCahn_solver_explicit(eta,S1,S2,S3,L1,L2,L3,Is1,Ie1,Is2,Ie2,Is3,Ie3,Id1,Id2,Id3, &
										left,right,front,back,top,bot,cartcomm)
										
	IF (ANY(ISNAN(eta))) THEN
		PRINT *, 'Error: Element in eta is NAN'
		STOP
	END IF
	
	time_elapsed = time_elapsed + dt
	
	!! -------------------------------------------------------------------------
	!! Output files
	IF (max_iter <= 10) THEN
		output_crit = 0
	ELSE  
		output_crit = MOD(iter,max_iter)
	END IF 
	
	IF (output_crit == 0) THEN
		IF (Id1 .eq. 0 .and. Id2 .eq. 0 .and. Id3 .eq. 0) PRINT *, 'Output number ', out_num, ' complete'
		
		CALL writeLotsFiles(eta,'eta_r',out_num,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)
		
		out_num = out_num + 1
	
	END IF
	
END DO

IF (Id1+Id2+Id3.eq. 0) PRINT *, 'Total time elapsed: ', max_iter*dt

END SUBROUTINE

!! ===================================================================================
!! Explicit Allen-Cahn Solver
!! ===================================================================================

SUBROUTINE allenCahn_solver_explicit(eta,S1,S2,S3,L1,L2,L3,Is1,Ie1,Is2,Ie2,Is3,Ie3,Id1,Id2,Id3, &
										left,right,front,back,top,bot,cartcomm)
USE parameters
IMPLICIT NONE

INCLUDE 'mpif.h'
INTEGER :: L1, L2, L3, S1, S2, S3, Id1, Id2, Id3, Is1, Is2, Is3, Ie1, Ie2, Ie3
INTEGER :: left,right,front,back,top,bot, cartcomm 	
REAL(KIND=8), DIMENSION(S1:S1+L1-1,S2:S2+L2-1,S3:S3+L3-1) :: c, eta, mu, dfdeta, landau_term
REAL(KIND=8), DIMENSION(S1+1:S1+L1-2,S2+1:S2+L2-2,S3+1:S3+L3-2) :: kappa_Lap_eta

dfdeta = 4.0d0*eta*(eta-1.0d0)*(eta-0.5)

kappa_Lap_eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) = kappa*( (eta(Is1+1:Ie1+1,Is2:Ie2,Is3:Ie3) - 2.0d0*eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) &
			+ eta(Is1-1:Ie1-1,Is2:Ie2,Is3:Ie3))/dx**2 &
			+ (eta(Is1:Ie1,Is2+1:Ie2+1,Is3:Ie3) - 2.0d0*eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) &
			+ eta(Is1:Ie1,Is2-1:Ie2-1,Is3:Ie3))/dy**2 &
			+ (eta(Is1:Ie1,Is2:Ie2,Is3+1:Ie3+1) - 2.0d0*eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) &
			+ eta(Is1:Ie1,Is2:Ie2,Is3-1:Ie3-1))/dz**2 )
			
eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) = eta(Is1:Ie1,Is2:Ie2,Is3:Ie3) &
			- dt * L * ( dfdeta(Is1:Ie1,Is2:Ie2,Is3:Ie3) - kappa_Lap_eta)

call commuBC(Id1,Id2,Id3,L1,L2,L3,S1,S2,S3,eta,left,right,front,back,top,bot,cartcomm)  	  
CALL external_BCs_Zero_Derivative(eta,.TRUE.,1,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,1,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.TRUE.,2,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,2,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.TRUE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)
CALL external_BCs_Zero_Derivative(eta,.FALSE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)


END SUBROUTINE allenCahn_solver_explicit


!! ===================================================================================

SUBROUTINE applyBCs(Id1,Id2,Id3,L1,L2,L3,S1,S2,S3,phi1,U_BC_top,left,right,front,back,top,bot,cartcomm,direction)
use parameters
implicit none
INCLUDE 'mpif.h'

INTEGER S1, S2, S3, Id1, Id2, Id3,L1,L2,L3
INTEGER :: left,right,front,back,top,bot,cartcomm !ranks of processors to shift to and cartcomm
INTEGER :: i,j,k, direction
INTEGER :: errcode
REAL(KIND=8) :: phi1(S1:S1+L1-1,S2:S2+L2-1,S3:S3+L3-1), U_BC_top
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

REAL(KIND=8), DIMENSION(1:L2,1:L3) :: Buffer_x1, Buffer_x2, Buffer_x3, Buffer_x4
REAL(KIND=8), DIMENSION(1:L1,1:L3) :: Buffer_y1, Buffer_y2, Buffer_y3, Buffer_y4
REAL(KIND=8), DIMENSION(1:L1,1:L2) :: Buffer_z1, Buffer_z2, Buffer_z3, Buffer_z4
INTEGER,DIMENSION(12) :: req

!! Communication in n1 direction

IF ( direction .EQ. 0 .OR. direction .EQ. 1 ) THEN

	! x-direction
	Buffer_x3 = phi1(S1+L1-1-bb,S2:S2+L2-1,S3:S3+L3-1)
	Buffer_x2 = phi1(S1+bb,S2:S2+L2-1,S3:S3+L3-1)		
	
	CALL MPI_iSEND(Buffer_x3,L2*L3,MPI_DOUBLE_PRECISION,right,r1Tag,cartcomm,req(1),errcode)
	CALL MPI_iRECV(Buffer_x4,L2*L3,MPI_DOUBLE_PRECISION,right,l1Tag,cartcomm,req(2),errcode)
	CALL MPI_iSEND(Buffer_x2,L2*L3,MPI_DOUBLE_PRECISION,left ,l1Tag,cartcomm,req(4),errcode)
	CALL MPI_iRECV(Buffer_x1,L2*L3,MPI_DOUBLE_PRECISION,left ,r1Tag,cartcomm,req(3),errcode)
	call MPI_waitall(4,req(1:4),MPI_STATUS_IGNORE,errcode)
	
	phi1(S1,S2:S2+L2-1,S3:S3+L3-1) = Buffer_x1
	phi1(S1+L1-1,S2:S2+L2-1,S3:S3+L3-1) = Buffer_x4

ENDIF	
	
!! Communication in n2 direction

IF ( direction .EQ. 0 .OR. direction .EQ. 2  )	THEN
	! y-direction
	Buffer_y3 = phi1(S1:S1+L1-1,S2+L2-1-bb,S3:S3+L3-1)
	Buffer_y2 = phi1(S1:S1+L1-1,S2+bb,S3:S3+L3-1)	
	CALL MPI_iSEND(Buffer_y3,L1*L3,MPI_DOUBLE_PRECISION,back ,b1Tag,cartcomm,req(5),errcode)
	CALL MPI_iRECV(Buffer_y4,L1*L3,MPI_DOUBLE_PRECISION,back ,f1Tag,cartcomm,req(6),errcode)
	CALL MPI_iSEND(Buffer_y2,L1*L3,MPI_DOUBLE_PRECISION,front,f1Tag,cartcomm,req(8),errcode)
	CALL MPI_iRECV(Buffer_y1,L1*L3,MPI_DOUBLE_PRECISION,front,b1Tag,cartcomm,req(7),errcode)
	call MPI_waitall(4,req(5:8),MPI_STATUS_IGNORE,errcode)
	phi1(S1:S1+L1-1,S2,S3:S3+L3-1) = Buffer_y1
	phi1(S1:S1+L1-1,S2+L2-1,S3:S3+L3-1) = Buffer_y4

ENDIF

!! Communication in n3 direction

IF ( direction .EQ. 0 .OR. direction .EQ. 3 ) THEN
	! z-direction
	Buffer_z3=phi1(S1:S1+L1-1,S2:S2+L2-1,S3+L3-1-bb)
	Buffer_z2=phi1(S1:S1+L1-1,S2:S2+L2-1,S3+bb)
	IF ( Id3 .NE. num_subdomain_z-1 ) THEN
		CALL MPI_iSEND(Buffer_z3,L1*L2,MPI_DOUBLE_PRECISION,top,t1Tag,cartcomm,req(9),errcode)
		CALL MPI_iRECV(Buffer_z4,L1*L2,MPI_DOUBLE_PRECISION,top,b1Tag,cartcomm,req(10),errcode)
	ELSE
		!CALL external_BCs_Dirichlet(phi1,U_BC_top,.TRUE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)
		Buffer_z4 = U_BC_top
	ENDIF
	
	IF ( Id3 .NE. 0 ) THEN
		CALL MPI_iSEND(Buffer_z2,L1*L2,MPI_DOUBLE_PRECISION,bot,b1Tag,cartcomm,req(12),errcode)
		CALL MPI_iRECV(Buffer_z1,L1*L2,MPI_DOUBLE_PRECISION,bot,t1Tag,cartcomm,req(11),errcode)
	ELSE
		!CALL external_BCs_Zero_Derivative(phi1,.FALSE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)
		Buffer_z1 = phi1(S1:S1+L1-1,S2:S2+L2-1,S3+bb)
	ENDIF
	
	IF ( Id3 .NE. num_subdomain_z-1 ) call MPI_waitall(2,req(9:10),MPI_STATUS_IGNORE,errcode)
	IF ( Id3 .NE. 0 ) call MPI_waitall(2,req(11:12),MPI_STATUS_IGNORE,errcode)
	phi1(S1:S1+L1-1,S2:S2+L2-1,S3)=Buffer_z1
	phi1(S1:S1+L1-1,S2:S2+L2-1,S3+L3-1)=Buffer_z4
	
ENDIF

!! External BCs in the z-direction
!CALL external_BCs_Dirichlet(phi1,U_BC_top,.TRUE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)
!CALL external_BCs_Zero_Derivative(phi1,.FALSE.,3,Id1,Id2,Id3,L1-2,L2-2,L3-2)

RETURN
END SUBROUTINE applyBCs

!! ===================================================================================
!! External BC calculation subroutines
!! ===================================================================================

SUBROUTINE external_BCs_Dirichlet(array,BC_value,is_top,dimension_num,Id1,Id2,Id3,px,py,pz)
USE parameters
IMPLICIT NONE

INTEGER :: px, py, pz
REAL(KIND=8), DIMENSION(0:px+1,0:py+1,0:pz+1) :: array
REAL(KIND=8) :: BC_value
INTEGER :: dimension_num,Id1,Id2,Id3
LOGICAL :: is_top

SELECT CASE (dimension_num)

CASE (1)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id1 .eq. num_subdomain_x-1) THEN
			array(px+1,0:py+1,0:pz+1) = BC_value
		END IF
	ELSE
		IF (Id1 .eq. 0) THEN
			array(0,0:py+1,0:pz+1) = BC_value
		END IF
	END IF
	
CASE (2)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id2 .eq. num_subdomain_y-1) THEN
			array(0:px+1,py+1,0:pz+1) = BC_value
		END IF
	ELSE
		IF (Id2 .eq. 0) THEN
			array(0:px+1,0,0:pz+1) = BC_value
		END IF
	END IF
	
CASE (3)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id3 .eq. num_subdomain_z-1) THEN
			array(0:px+1,0:py+1,pz+1) = BC_value
		END IF
	ELSE
		IF (Id3 .eq. 0) THEN
			array(0:px+1,0:py+1,0) = BC_value
		END IF
	END IF
	
END SELECT

!PRINT *, 'in BC subroutine', array(40, 10, pz+1),  array(40, 10, 0), top_val, bottom_val

END SUBROUTINE external_BCs_Dirichlet

!! ===================================================================================

SUBROUTINE external_BCs_Zero_Derivative(array,is_top,dimension_num,Id1,Id2,Id3,px,py,pz)
USE parameters
IMPLICIT NONE

INTEGER :: px, py, pz
REAL(KIND=8), DIMENSION(0:px+1,0:py+1,0:pz+1) :: array
REAL(KIND=8) :: BC_value
INTEGER :: dimension_num,Id1,Id2,Id3
LOGICAL :: is_top

SELECT CASE (dimension_num)

CASE (1)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id1 .eq. num_subdomain_x-1) THEN
			array(px+1,0:py+1,0:pz+1) = array(px,0:py+1,0:pz+1)
		END IF
	ELSE
		IF (Id1 .eq. 0) THEN
			array(0,0:py+1,0:pz+1) = array(1,0:py+1,0:pz+1)
		END IF
	END IF
	
CASE (2)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id2 .eq. num_subdomain_y-1) THEN
			array(0:px+1,py+1,0:pz+1) = array(0:px+1,py,0:pz+1)
		END IF
	ELSE
		IF (Id2 .eq. 0) THEN
			array(0:px+1,0,0:pz+1) = array(0:px+1,1,0:pz+1)
		END IF
	END IF
	
CASE (3)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id3 .eq. num_subdomain_z-1) THEN
			array(0:px+1,0:py+1,pz+1) = array(0:px+1,0:py+1,pz)
		END IF
	ELSE
		IF (Id3 .eq. 0) THEN
			array(0:px+1,0:py+1,0) = array(0:px+1,0:py+1,1)
		END IF
	END IF
	
END SELECT

!PRINT *, 'in BC subroutine', array(40, 10, pz+1),  array(40, 10, 0), top_val, bottom_val

END SUBROUTINE external_BCs_Zero_Derivative

!! ===================================================================================

SUBROUTINE external_BCs_Neumann(array,derivative,is_top,dimension_num,Id1,Id2,Id3,px,py,pz)
USE parameters
IMPLICIT NONE

INTEGER :: px, py, pz
REAL(KIND=8), DIMENSION(0:px+1,0:py+1,0:pz+1) :: array
REAL(KIND=8) :: derivative
INTEGER :: dimension_num,Id1,Id2,Id3
LOGICAL :: is_top

SELECT CASE (dimension_num)

CASE (1)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id1 .eq. num_subdomain_x-1) THEN
			array(px+1,0:py+1,0:pz+1) = array(px,0:py+1,0:pz+1) + derivative
		END IF
	ELSE
		IF (Id1 .eq. 0) THEN
			array(0,0:py+1,0:pz+1) = array(1,0:py+1,0:pz+1) - derivative
		END IF
	END IF
	
CASE (2)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id2 .eq. num_subdomain_y-1) THEN
			array(0:px+1,py+1,0:pz+1) = array(0:px+1,py,0:pz+1) + derivative
		END IF
	ELSE
		IF (Id2 .eq. 0) THEN
			array(0:px+1,0,0:pz+1) = array(0:px+1,1,0:pz+1) - derivative
		END IF
	END IF
	
CASE (3)
	IF (is_top .eqv. .TRUE.) THEN
		IF (Id3 .eq. num_subdomain_z-1) THEN
			array(0:px+1,0:py+1,pz+1) = array(0:px+1,0:py+1,pz) + derivative
		END IF
	ELSE
		IF (Id3 .eq. 0) THEN
			array(0:px+1,0:py+1,0) = array(0:px+1,0:py+1,1) - derivative
		END IF
	END IF
	
END SELECT

END SUBROUTINE external_BCs_Neumann

!! ===================================================================================

SUBROUTINE external_BCs_Periodic(array,dimension_num,Id1,Id2,Id3,px,py,pz,left,right,front,back,top,bot,cartcomm)
USE parameters
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER :: px, py, pz
REAL(KIND=8), DIMENSION(0:px+1,0:py+1,0:pz+1) :: array
INTEGER :: dimension_num,Id1,Id2,Id3
INTEGER :: left,right,front,back,top,bot,cartcomm
INTEGER :: errcode
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

SELECT CASE (dimension_num)

CASE (1)
	IF ( num_subdomain_x .eq. 1) THEN
		array(0,0:py+1,0:pz+1) = array(px,0:py+1,0:pz+1)
		array(px+1,0:py+1,0:pz+1) = array(1,0:py+1,0:pz+1)
	ELSE		
		IF (Id1 .eq. num_subdomain_x-1) &
		call MPI_SEND(array(px,0:py+1,0:pz+1),(py+2)*(pz+2),MPI_DOUBLE_PRECISION,right,r1Tag+100,cartcomm,errcode)
	
		IF (Id1 .eq. 0) &
		call MPI_RECV(array(0,0:py+1,0:pz+1),(py+2)*(pz+2),MPI_DOUBLE_PRECISION,left,r1Tag+100,cartcomm,nstatus,errcode)
	
		IF (Id1 .eq. 0) &
		call MPI_SEND(array(1,0:py+1,0:pz+1),(py+2)*(pz+2),MPI_DOUBLE_PRECISION,left,l1Tag+100,cartcomm,errcode)
	
		IF (Id1 .eq. num_subdomain_x-1) &
		call MPI_RECV(array(px+1,0:py+1,0:pz+1),(py+2)*(pz+2),MPI_DOUBLE_PRECISION,right,l1Tag+100,cartcomm,nstatus,errcode)
	END IF
			
CASE (2)
	IF ( num_subdomain_y .eq. 1) THEN
		array(0:px+1,0,0:pz+1) = array(0:px+1,py,0:pz+1)
		array(0:px+1,py+1,0:pz+1) = array(0:px+1,1,0:pz+1)
	ELSE	
		IF (Id2 .eq. num_subdomain_y-1) THEN
			call MPI_SEND(array(0:px+1,py,0:pz+1),(px+2)*(pz+2),MPI_DOUBLE_PRECISION,back,f1Tag+100,cartcomm,errcode)
		END IF

		IF (Id2 .eq. 0) THEN
			call MPI_RECV(array(0:px+1,0,0:pz+1),(px+2)*(pz+2),MPI_DOUBLE_PRECISION,front,f1Tag+100,cartcomm,nstatus,errcode)
		END IF
	
		IF (Id2 .eq. 0) &
		call MPI_SEND(array(0:px+1,1,0:pz+1),(px+2)*(pz+2),MPI_DOUBLE_PRECISION,front,b1Tag+100,cartcomm,errcode)
	
		IF (Id2 .eq. num_subdomain_y-1) &
		call MPI_RECV(array(0:px+1,py+1,0:pz+1),(px+2)*(pz+2),MPI_DOUBLE_PRECISION,back,b1Tag+100,cartcomm,nstatus,errcode)
	END IF
CASE (3)	
	IF ( num_subdomain_z .eq. 1) THEN
		array(0:px+1,0:py+1,0) = array(0:px+1,0:py+1,pz)
		array(0:px+1,0:py+1,pz+1) = array(0:px+1,0:py+1,1)
	ELSE
		IF (Id3 .eq. num_subdomain_z-1) &
		call MPI_SEND(array(0:px+1,0:py+1,pz),(px+2)*(py+2),MPI_DOUBLE_PRECISION,top,t1Tag+100,cartcomm,errcode)
	
		IF (Id3 .eq. 0)	&
		call MPI_RECV(array(0:px+1,0:py+1,0),(px+2)*(py+2),MPI_DOUBLE_PRECISION,bot,t1Tag+100,cartcomm,nstatus,errcode)
	
		IF (Id3 .eq. 0) &
		call MPI_SEND(array(0:px+1,0:py+1,1),(px+2)*(py+2),MPI_DOUBLE_PRECISION,bot,d1Tag+100,cartcomm,errcode)
	
		IF (Id3 .eq. num_subdomain_z-1) &
		call MPI_RECV(array(0:px+1,0:py+1,pz+1),(px+2)*(py+2),MPI_DOUBLE_PRECISION,top,d1Tag+100,cartcomm,nstatus,errcode)
	END IF
		
END SELECT

END SUBROUTINE external_BCs_Periodic

!!===================================================================================
!! Internal BC Communication Subroutine
!!===================================================================================

SUBROUTINE commuBC(Id1,Id2,Id3,L1,L2,L3,S1,S2,S3,phi1,left,right,front,back,top,bot,cartcomm)
use parameters
implicit none
INCLUDE 'mpif.h'

INTEGER S1, S2, S3, Id1, Id2, Id3,L1,L2,L3
INTEGER :: left,right,front,back,top,bot,cartcomm !ranks of processors to shift to and cartcomm
INTEGER :: i,j,k
INTEGER :: errcode
REAL(KIND=8) :: phi1(S1:S1+L1-1,S2:S2+L2-1,S3:S3+L3-1)
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

! Interior boundary condition
! communication in n1 direction
	
If(Id1 .NE. num_subdomain_x-1) call MPI_SEND(phi1(S1+L1-1-bb,S2:S2+L2-1,S3:S3+L3-1), &
												bb*L2*L3,MPI_DOUBLE_PRECISION,right,r1Tag,cartcomm,errcode)    !Send real grid  
												
If(Id1 .NE. 0) call MPI_RECV(phi1(S1,S2:S2+L2-1,S3:S3+L3-1), &
												bb*L2*L3,MPI_DOUBLE_PRECISION,left,r1Tag,cartcomm,nstatus,errcode)          !receive by ghost grid

If(Id1 .NE. 0) call MPI_SEND(phi1(S1+bb,S2:S2+L2-1,S3:S3+L3-1), &
												bb*L2*L3,MPI_DOUBLE_PRECISION,left,l1Tag,cartcomm,errcode)

If(Id1 .NE. num_subdomain_x-1) call MPI_RECV(phi1(S1+L1-1,S2:S2+L2-1,S3:S3+L3-1), &
												bb*L2*L3,MPI_DOUBLE_PRECISION,right,l1Tag,cartcomm,nstatus,errcode)

!  communication in n2 direction
   
If(Id2 .NE. num_subdomain_y-1) call MPI_SEND(phi1(S1:S1+L1-1,S2+L2-1-bb,S3:S3+L3-1), &
												L1*bb*L3,MPI_DOUBLE_PRECISION,back,f1Tag,cartcomm,errcode)

If(Id2 .NE. 0) call MPI_RECV(phi1(S1:S1+L1-1,S2,S3:S3+L3-1), &
												L1*bb*L3,MPI_DOUBLE_PRECISION,front,f1Tag,cartcomm,nstatus,errcode)

If(Id2 .NE. 0) call MPI_SEND(phi1(S1:S1+L1-1,S2+bb,S3:S3+L3-1), &
												L1*bb*L3,MPI_DOUBLE_PRECISION,front,b1Tag,cartcomm,errcode)

If(Id2 .NE. num_subdomain_y-1) call MPI_RECV(phi1(S1:S1+L1-1,S2+L2-1,S3:S3+L3-1), &
												L1*bb*L3,MPI_DOUBLE_PRECISION,back,b1Tag,cartcomm,nstatus,errcode)

! communication in n3 direction

  
If(Id3 .NE. num_subdomain_z-1) call MPI_SEND(phi1(S1:S1+L1-1,S2:S2+L2-1,S3+L3-1-bb), &
												L1*L2*bb,MPI_DOUBLE_PRECISION,top,t1Tag,cartcomm,errcode)

If(Id3 .NE. 0) call MPI_RECV(phi1(S1:S1+L1-1,S2:S2+L2-1,S3), &
												L1*L2*bb,MPI_DOUBLE_PRECISION,bot,t1Tag,cartcomm,nstatus,errcode)

If(Id3 .NE. 0) call MPI_SEND(phi1(S1:S1+L1-1,S2:S2+L2-1,S3+bb), &
												L1*L2*bb,MPI_DOUBLE_PRECISION,bot,d1Tag,cartcomm,errcode)

If(Id3 .NE. num_subdomain_z-1) call MPI_RECV(phi1(S1:S1+L1-1,S2:S2+L2-1,S3+L3-1), &
												L1*L2*bb,MPI_DOUBLE_PRECISION,top,d1Tag,cartcomm,nstatus,errcode)

RETURN
END SUBROUTINE 

!! ===================================================================================
!! Input subroutine
!! ===================================================================================

!SUBROUTINE loadLotsFiles(arrIn,name_string,Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3)
!! Purpose:
!!   Write arrIn to an individual file for each core.
!!	File has Id1,Id2,Id3 in filename as _XXYYZZ.frt
!
!USE parameters
!IMPLICIT NONE
!
!! Declare variables:--------------------------------------------------------+
!INTEGER, INTENT(in) :: Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3
!REAL(KIND=8), DIMENSION(Is1-bb:Ie1+bb,Is2-bb:Ie2+bb,Is3-bb:Ie3+bb) :: arrIn     ! Array to load
!CHARACTER(LEN=29) :: filename                                   ! Filename of input or output
!CHARACTER(LEN=5) :: tstring, name_string
!CHARACTER(LEN=2) :: Id1string, Id2string, Id3string
!CHARACTER(LEN=3) :: counter
!INTEGER :: status
!INTEGER :: tpStep, i, outFreq
!!---------------------------------------------------------------------------+
!
!WRITE(counter,'(I3)') 100+load_from_output
!
!write(Id1string,'(i2.2)')Id1
!write(Id2string,'(i2.2)')Id2
!write(Id3string,'(i2.2)')Id3
!
!filename=name_string//load_from_run//'_'//counter(2:3)//'_'//Id1string//Id2string//Id3string//'.frt'
!!PRINT *,'Filename:', filename
!OPEN (FILE=output_path//filename, UNIT=300, STATUS='UNKNOWN', FORM='UNFORMATTED')
!READ (300) arrIn(Is1:Ie1,Is2:Ie2,Is3:Ie3)
!CLOSE (300)
!
!RETURN
!
!END SUBROUTINE loadLotsFiles


!! ===================================================================================
!! Output subroutine
!! ===================================================================================

SUBROUTINE writeLotsFiles(arrIn,name_string,step,Is1,Ie1,Is2,Ie2,Is3,Ie3,Id1,Id2,Id3)
! Purpose:
!   Write arrIn to an individual file for each core.
!	File has Id1,Id2,Id3 in filename as _XXYYZZ.frt

USE parameters
IMPLICIT NONE

! Declare variables:--------------------------------------------------------+
INTEGER, INTENT(in) :: Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3
REAL(KIND=8), DIMENSION(Is1-bb:Ie1+bb,Is2-bb:Ie2+bb,Is3-bb:Ie3+bb), INTENT(IN) :: arrIn     ! Array to write
CHARACTER(LEN=29) :: filename                                   ! Filename of input or output
CHARACTER(LEN=5) :: tstring, name_string
CHARACTER(LEN=2) :: Id1string, Id2string, Id3string
CHARACTER(LEN=3) :: counter
INTEGER, INTENT(in) :: step
INTEGER :: tpStep, i, outFreq
!---------------------------------------------------------------------------+

WRITE(counter,'(I3)') 100+step

write(Id1string,'(i2.2)')Id1
write(Id2string,'(i2.2)')Id2
write(Id3string,'(i2.2)')Id3

filename=name_string//run_num//'_'//counter(2:3)//'_'//Id1string//Id2string//Id3string//'.frt'
OPEN (FILE=output_path//filename, UNIT=300, STATUS='UNKNOWN', FORM='UNFORMATTED')
WRITE (300) arrIn(Is1:Ie1,Is2:Ie2,Is3:Ie3)
CLOSE (300)

RETURN

END SUBROUTINE writeLotsFiles

!! ===================================================================================

SUBROUTINE writeLotsFiles2(arrIn,name_string,step,Is1,Ie1,Is2,Ie2,Is3,Ie3,Id1,Id2,Id3)
! Purpose:
!   Write arrIn to an individual file for each core.
!	File has Id1,Id2,Id3 in filename as _XXYYZZ.frt

USE parameters
IMPLICIT NONE

! Declare variables:--------------------------------------------------------+
INTEGER, INTENT(in) :: Is1,Ie1,Is2,Ie2,Is3,Ie3, Id1,Id2,Id3
REAL(KIND=8), DIMENSION(Is1:Ie1,Is2:Ie2,Is3:Ie3), INTENT(IN) :: arrIn     ! Array to write
CHARACTER(LEN=29) :: filename                                   ! Filename of input or output
CHARACTER(LEN=5) :: tstring, name_string
CHARACTER(LEN=2) :: Id1string, Id2string, Id3string
CHARACTER(LEN=3) :: counter
INTEGER, INTENT(in) :: step
INTEGER :: tpStep, i, outFreq
!---------------------------------------------------------------------------+

WRITE(counter,'(I3)') 100+step

write(Id1string,'(i2.2)')Id1
write(Id2string,'(i2.2)')Id2
write(Id3string,'(i2.2)')Id3

filename=name_string//run_num//'_'//counter(2:3)//'_'//Id1string//Id2string//Id3string//'.frt'
OPEN (FILE=output_path//filename, UNIT=300, STATUS='UNKNOWN', FORM='UNFORMATTED')
WRITE (300) arrIn(Is1:Ie1,Is2:Ie2,Is3:Ie3)
CLOSE (300)

RETURN

END SUBROUTINE writeLotsFiles2

!! =================================================================================================

SUBROUTINE outputVector(vector,outName,length)
USE parameters
IMPLICIT NONE
INTEGER :: length, outNum
REAL(KIND=8),DIMENSION(1:length) :: vector
CHARACTER(LEN=8) :: outName
CHARACTER(LEN=12) :: filename

!! This subroutine outputs arrays of length "length" to file.

WRITE(filename, '(A8,A4)') outName,'.dat'
OPEN(UNIT=300,FILE=output_path//filename,FORM='FORMATTED',ACTION='WRITE')
WRITE(300,'(E25.12)') vector
CLOSE(300)

END SUBROUTINE outputVector

!! =================================================================================================

!SUBROUTINE loadVector(vector,outName,length)
!USE parameters
!IMPLICIT NONE
!INTEGER :: length, outNum
!REAL(KIND=8),DIMENSION(1:length) :: vector
!CHARACTER(LEN=8) :: outName
!CHARACTER(LEN=12) :: filename
!
!!! This subroutine outputs arrays of length "length" to file.
!
!WRITE(filename, '(A8,A4)') outName,'.dat'
!OPEN(UNIT=300,FILE=output_path//filename,FORM='FORMATTED',ACTION='READ')
!READ(300,'(E25.12)') vector
!CLOSE(300)
!
!END SUBROUTINE loadVector

!! =================================================================================================

SUBROUTINE outputDomainInfo()
USE parameters
IMPLICIT NONE

CHARACTER(LEN=3) :: counter
CHARACTER(LEN=18) :: filename

!! This subroutine outputs the dimensions of the domain and the number of subdomains

WRITE(filename, '(A11,A3,A4)') 'domain_info',run_num,'.frt'
OPEN(UNIT=20,FILE=output_path//filename,ACTION='WRITE')
WRITE(20,*) xg
WRITE(20,*) yg
WRITE(20,*) zg
WRITE(20,*) num_subdomain_x
WRITE(20,*) num_subdomain_y
WRITE(20,*) num_subdomain_z
CLOSE(20)

END SUBROUTINE
	

