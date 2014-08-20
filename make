
cat mc_nuclear_density_fun.f nuclear_density_function.f program.f RAN2.f single_evt_glauber.f > all.f
gfortran --fixed-line-length-none all.f
