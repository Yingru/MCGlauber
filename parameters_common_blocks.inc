!	constant parameters
	
	double precision PI	
	parameter (PI = 3.14159265358979324d0)

	double precision ran2

	integer idum
	common /idum/ idum

	integer N_evt
	parameter (N_evt = 3000) 

!	nucleon size
	double precision nucleon_dr
	parameter (nucleon_dr = 0.2d0)

	integer atom_A, atom_B
!	parameter (atom_A = 207, atom_B = 207)
	parameter (atom_A = 197, atom_B = 197)
!	parameter (atom_A = 63, atom_B = 63)

	double precision radius
!	parameter (radius = 6.62d0) ! for Pb
	parameter (radius = 6.38d0) ! for Au
!	parameter (radius = 4.20641d0) ! for Cu

!	sigma_NN = 42mb for sqrt(s) = 200GeV
!	sigma_NN = 41mb for sqrt(s) = 130GeV
!	sigma_NN = 36mb for sqrt(s) = 62.4GeV
!	sigma_NN = 32mb for sqrt(s) = 19.6GeV

	double precision sqrt_sNN, sigma_NN
!	parameter (sqrt_sNN = 2760.d0, sigma_NN = 6.4d0)
!	parameter (sqrt_sNN = 200.d0, sigma_NN = 4.2d0)
	parameter (sqrt_sNN = 62.4d0, sigma_NN = 3.6d0)

	double precision frac_coll
!	parameter (frac_coll = 0.11d0)
	parameter (frac_coll = 0.13d0)

!	the parameters in negative binomial distribution for pp multiplicity
	double precision mean_npp, k_pp
	parameter (mean_npp = 2.35d0, k_pp = 1.9d0)	! My favorite one
!	parameter (mean_npp = 2.34d0, k_pp = 2.d0)	! STAR: Jon Gan's thesis
!	parameter (mean_npp = 2.55d0, k_pp = 1.86d0)	! A.K. Dash & B. Mohanty: arxi:0908.0888
!	parameter (mean_npp = 2.48d0, k_pp = 2.3d0)	! UA5 p-pbar analysis
!	parameter (mean_npp = 0.32d0, k_pp = 1.88d0)	! PHENIX paper: PRC 78, 044902(2008), with limited acceptance

!	LHC: 7TeV
!	parameter (mean_npp = 6.0d0, k_pp = 1.1d0)	! My favorite one: rapidity cut
!	LHC: 2.36TeV
!	parameter (mean_npp = 4.47d0, k_pp = 1.3d0)	! My favorite one: rapidity cut
!	LHC: 0.9TeV
!	parameter (mean_npp = 3.6d0, k_pp = 1.7d0)	! My favorite one: rapidity cut


!	the flag of Glauber model
	integer flag_glauber
	common /flag_glauber/ flag_glauber

!	we first get the non-contracted profile, then transform to pancake-like profile by shrinking z-direction
	double precision gamma_pancake
	parameter (gamma_pancake = 1380.d0)
!	parameter (gamma_pancake = 100.d0)
!	parameter (gamma_pancake = 31.2d0)
!	parameter (gamma_pancake = 1.d0)


!!!!!!!!!!!!!!!!!!!!!!!

! 	commmon blocks

	double precision x_n1, y_n1, z_n1, x_n2, y_n2, z_n2
	common /xyz_part/ x_n1(atom_A), y_n1(atom_A), z_n1(atom_A), x_n2(atom_B), y_n2(atom_B), z_n2(atom_B)
	integer part_n1, part_n2
	common /if_part/ part_n1(atom_A), part_n2(atom_B)
	double precision nch_n1, nch_n2
	common /nch_part/ nch_n1(atom_A), nch_n2(atom_B)

!	double precision xxx_n1, yyy_n1, zzz_n1, xxx_n2, yyy_n2, zzz_n2	
!	common /participant_positions/ xxx_n1(n_evt, atom_A), yyy_n1(n_evt, atom_A), zzz_n1(n_evt, atom_A), xxx_n2(n_evt, atom_B), yyy_n2(n_evt, atom_B), zzz_n2(n_evt, atom_B)
!	integer participant_n1, participant_n2
!	common /if_participant/ participant_n1(n_evt, atom_A), participant_n2(n_evt, atom_B)

	double precision x_coll, y_coll, z_coll
	common /xyz_coll/ x_coll(atom_A, atom_B), y_coll(atom_A, atom_B), z_coll(atom_A, atom_B)
	integer coll_n1_n2
	common /if_coll/ coll_n1_n2(atom_A, atom_B)
	double precision nch_n1_n2
	common /nch_coll/ nch_n1_n2(atom_A, atom_B)

!	double precision xxx_coll, yyy_coll, zzz_coll
!	common /collision_positions/ xxx_coll(n_evt, atom_A, atom_B), yyy_coll(n_evt, atom_A, atom_B), zzz_coll(n_evt, atom_A, atom_B)
!	integer collision_n1_n2
!	common /if_collision/ collision_n1_n2(n_evt, atom_A, atom_B)

