
	subroutine mc_rho_A(x, y, z, atom)	! generate nucleon position (x, y , z) according to rho_A(r) distribution (nuclear density function)
	
	implicit none

	include "parameters_common_blocks.inc"
	
	double precision x, y, z
	double precision cos_theta, sin_theta, phi
	integer atom

	double precision rho_A, r
	
	double precision rmax_ws

	double precision rdn1, rdn2, rdn3, rdn4
	double precision accept_prob 


	r = 0.d0
	rmax_ws = 2.d0 * radius

100	rdn1 = ran2(idum)
	r = rmax_ws * rdn1**(1.d0/3.d0)
c       what does this mean? why rdn1**(1.0/3)?????

	rdn2 = ran2(idum)
	accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)

	if (rdn2 .GT. accept_prob) go to 100

!	do 
!		rdn1 = ran2(idum)
!		r = rmax_ws * rdn1**(1.d0/3.d0)
!		rdn2 = ran2(idum)
!		accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)
!		if (rdn2 .LE. accept_prob) exit
!	end do

!	rdn2 = 1.d0
!	accept_prob = 0.d0
!
!	do while (rdn2 .GT. accept_prob) 
!		rdn1 = ran2(idum)
!		r = rmax_ws * rdn1**(1.d0/3.d0)
!		rdn2 = ran2(idum)
!		accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)
!	end do

	rdn3 = ran2(idum)
	cos_theta = 1.d0 - 2.d0 * rdn3
	sin_theta = sqrt(1.d0 - cos_theta**2)

	rdn4 = ran2(idum)
	phi = 2.d0 * PI * rdn4

	x = r * sin_theta * cos(phi)
	y = r * sin_theta * sin(phi)
	z = r * cos_theta
	
	end subroutine



