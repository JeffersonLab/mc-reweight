#/usr/bin/tcsh
gfortran -c callmod.f model_new.f f1f209.f rc_mod.f
gfortran -o callmod callmod.o model_new.o f1f209.o rc_mod.o
#
