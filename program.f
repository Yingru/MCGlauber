
        program main

        implicit none

        include "parameters_common_blocks.inc"

        double precision temp 
        double precision impact_b
        double precision x, y, z
        double precision r, phi
        double precision rho_A
        integer atom
        double precision dr
        integer i_atom, i_n1, i_n2 
        double precision N_coll, N_part, N_charge
        double precision Nsq_coll, Nsq_part, Nsq_charge
        double precision sigma_Ncoll, sigma_Npart, sigma_Ncharge
        double precision Ncoll_evt, Npart_evt, Ncharge_evt
        integer i_evt, ntr_evt, ntot_evt, ntr_ptc_evt
        double precision sigma_xsq, sigma_ysq, sigma_xy
        double precision ecc_std, ecc_pp
        double precision eccsq_std, eccsq_rp
        double precision sigma_eccstd, sigma_eccpp
        double precision eccstd_evt, eccpp_evt
        double precision dsigma_db
        double precision dsigmadb_evt
        double precision dsigma_db2
        double precision mean_x, mean_y
        double precision meanx_evt, meany_evt
        double precision b_min, b_max
        integer iosend,count_pt
        character (len = 100) filename

        count_pt=0

!	qin
9999    open (99999,  file="/dev/urandom", 
     $       access="stream", form="unformatted", iostat = iosend)

        if (iosend .eq. 0) then
                read (99999) idum
                write (*,*) "random seed from noise!"
        else
                idum = time()
                write (*,*) "random seed from time!"
        end if
        if (idum .eq. 0) go to 9999
        if (idum .gt. 0) idum = -idum
        write (*,*) idum
!	idum = -123456789

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (9, file="dt_impact_b.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (10, file="dt_error_messages.dat")

!	output nuclei profiles

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (11, file="dt_nucleus1.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (12, file="dt_part_nucleus1.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (13, file="dt_nucleus2.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (14, file="dt_part_nucleus2.dat")

        open (16, file="dt_binary_collisions.dat")
        open (17, file="dt_all_binary_collisions.dat")
        open (18, file="dt_all_binary_collisions_rphi.dat")

!	output npart, ncoll, ncharge, eccentricity

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
!	open (20, file=filename)
        open (20, file="dt_ncoll_npart_ncharge.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
        open (21, file="dt_dsigma_db.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
!	open (22, file=filename)
        open (22, file="dt_ncharge_per_pair.dat")

        write (*,*) "please give the filename for data output: " 
! 	read (*,*) filename
!	open (23, file=filename)
        open (23, file="dt_eccentricity.dat")

!	beginning of code

        dr = nucleon_dr         ! nucleon size

!	initialize the nuclear density distribution
        atom = atom_A
        r = 0.d0
        temp = rho_A(r, atom)
!	write (*,*) "nuclear distribution function:", atom, r, temp

!	read (*,*) flag_glauber
        flag_glauber = 21
        write (*,*) "model:", flag_glauber

        b_min = 0.d0
        b_max = 6.56

!	do loop in impact parameters
!	impact_b = 0.d0
!	do while (impact_b .le. 2.d0*radius)

!		initialization before calling glauber

                ntr_evt = 0
                ntot_evt = 0

                N_coll = 0.d0
                N_part = 0.d0
                N_charge = 0.d0

                Nsq_coll = 0.d0
                Nsq_part = 0.d0
                Nsq_charge = 0.d0

                ecc_std = 0.d0
                ecc_pp = 0.d0

                eccsq_std = 0.d0
                eccsq_rp = 0.d0

                dsigma_db = 0.d0

!		do loop in events
                do i_evt = 1, n_evt

                        if (b_min .eq. b_max) then
                                impact_b = b_min
                        else
                                temp = ran2(idum)
                                impact_b = sqrt(b_min**2 + (b_max**2 - b_min**2) * temp)
c                               dont understand why here there's nothing to do ith radius? but a number in the range of (0,5)
                        end if

!                        write (*,*) "one glauber event:", i_evt

!			write (*,*) temp, impact_b
                        write (9,*) temp, impact_b

                        call single_evt_glauber(impact_b, Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccpp_evt, dsigmadb_evt)

!			write (*,*) Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccpp_evt

!			output nucleon positions and if they participate the collision
                        if (i_evt .eq. 1) then
                                do i_atom = 1, atom_A

                                        write (11, *) x_n1(i_atom), y_n1(i_atom)
                                        if (part_n1(i_atom) .eq. 1) write (12, *) x_n1(i_atom), y_n1(i_atom)
c                                       Is this recording the position of before the collision??? ( I don't think it is the positiona after the collision, whcih should be x(n1,n2)

                                end do

                                do i_atom = 1, atom_B

                                        write (13, *) x_n2(i_atom), y_n2(i_atom)
                                        if (part_n2(i_atom) .eq. 1) write (14, *) x_n2(i_atom), y_n2(i_atom)

                                end do
                        end if

!			count binary collisions and calculate moments of binary collision distribution
                        do i_n1 = 1, atom_A
                                do i_n2 = 1, atom_B

                                        x = x_coll(i_n1, i_n2)
                                        y = y_coll(i_n1, i_n2)
                                        z = z_coll(i_n1, i_n2)

                                        if (coll_n1_n2(i_n1, i_n2) .EQ. 1) then

                                                r = sqrt(x*x + y*y)
                                                phi = atan2(y, x)
                                                if (phi .lt. 0.d0) phi = phi + 2.d0*PI

                                                if (i_evt .eq. 1) write (16,*) x, y, r, phi

                                                write (17,33) x, y
                                                count_pt=count_pt+1
                                                write (18,*) r, phi
                                        end if

                                end do
                        end do

33                      format(E12.6,2X,E12.6)

!			record every event for later use
!			do i_atom = 1, atom_A
!
!				xxx_n1(i_evt, i_atom) = x_n1(i_atom)
!				yyy_n1(i_evt, i_atom) = y_n1(i_atom)
!				zzz_n1(i_evt, i_atom) = z_n1(i_atom)
!
!				participant_n1(i_evt, i_atom) = part_n1(i_atom)
!
!			end do
!
!			do i_atom = 1, atom_B
!
!				xxx_n2(i_evt, i_atom) = x_n2(i_atom)
!				yyy_n2(i_evt, i_atom) = y_n2(i_atom)
!				zzz_n2(i_evt, i_atom) = z_n2(i_atom)
!
!				participant_n2(i_evt, i_atom) = part_n2(i_atom)
!
!			end do
!
!			do i_n1 = 1, atom_A
!
!				do i_n2 = 1, atom_B
!
!				xxx_coll(i_evt, i_n1, i_n2) = x_coll(i_n1, i_n2)
!				yyy_coll(i_evt, i_n1, i_n2) = y_coll(i_n1, i_n2)
!				zzz_coll(i_evt, i_n1, i_n2) = z_coll(i_n1, i_n2)
!
!				collision_n1_n2(i_evt, i_n1, i_n2) = coll_n1_n2(i_n1, i_n2)
!
!				end do
!			end do

                        N_coll = N_coll + Ncoll_evt
                        N_part = N_part + Npart_evt
                        N_charge = N_charge + Ncharge_evt

                        Nsq_coll = Nsq_coll + Ncoll_evt**2
                        Nsq_part = Nsq_part + Npart_evt**2
                        Nsq_charge = Nsq_charge + Ncharge_evt**2

                        dsigma_db = dsigma_db + dsigmadb_evt

                        if (Npart_evt .GT. 0.d0 .AND. Ncharge_evt .GT. 0.d0) then

                                ntr_evt = ntr_evt + 1
                                ecc_std = ecc_std + eccstd_evt
                                ecc_pp = ecc_pp + eccpp_evt

                                eccsq_std = eccsq_std + eccstd_evt**2
                                eccsq_rp = eccsq_rp + eccpp_evt**2

                        end if

                        ntot_evt = ntot_evt + 1

                end do

!		average over events

                N_coll = N_coll / ntot_evt
                N_part = N_part / ntot_evt
                N_charge = N_charge / ntot_evt

                Nsq_coll = Nsq_coll / ntot_evt
                Nsq_part = Nsq_part / ntot_evt
                Nsq_charge = Nsq_charge / ntot_evt

                sigma_Ncoll = sqrt(Nsq_coll - N_coll**2)
                sigma_Npart = sqrt(Nsq_part - N_part**2)
                sigma_Ncharge = sqrt(Nsq_charge - N_charge**2)

!		the differential cross section: Monte Carlo result & the optical limit
                dsigma_db = dsigma_db / ntot_evt
                dsigma_db2 = 2.d0*PI*impact_b * (1.d0 - (1.d0 - N_coll / (atom_A*atom_B))**(atom_A*atom_B))

!		for eccenetricity we need to average over true events (if no collision happens, no medium, eccentricity is undetermined)
                ecc_std = ecc_std / ntr_evt
                ecc_pp = ecc_pp / ntr_evt

                eccsq_std = eccsq_std / ntr_evt
                eccsq_rp = eccsq_rp / ntr_evt

                sigma_eccstd = sqrt(eccsq_std - ecc_std**2)
                sigma_eccpp = sqrt(eccsq_rp - ecc_pp**2)

                write (*,*) impact_b, ntr_evt, N_coll, N_part, N_charge
!		write (*,*) impact_b, dsigma_db, dsigma_db2 
                write (*,*) impact_b, N_part, ecc_std, ecc_pp

                write (20,*) impact_b, N_coll, sigma_Ncoll, N_part, sigma_Npart, N_charge, sigma_Ncharge
                write (21,*) impact_b, dsigma_db, dsigma_db2 
                write (22,*) impact_b, N_part, N_charge/(N_part/2.d0)
                write (23,*) impact_b, N_part, ecc_std, sigma_eccstd, ecc_pp, sigma_eccpp 

                write(6,*) "count_pt: ",count_pt,count_pt*1d0/n_evt
!		impact_b = impact_b + 5.d0*dr

!	end do

        end program
