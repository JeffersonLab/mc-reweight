#/usr/bin/tcsh

# Link required Fortran files
echo
echo "Linking callmod and recon_mc from source code"
unlink callmod; unlink recon_mc
ln -fs SRC/callmod; ln -fs SRC/recon_mc
echo "Done."

# Define the Monte Carlo (mc) and radiative correction (rc) input files
#set mc_input_file = input/recon-mc/input_c1_30_1.1.dat
#set rc_input_file = input/rad-corr/rc94_c1_2.3.dat
set mc_input_file = input/recon-mc/kpp_shms_488_input.dat
set rc_input_file = input/rad-corr/kpp_shms_488_rc.dat

echo 
echo "Linking $mc_input_file as mc_input.dat"
unlink mc_input.dat
ln -fs $mc_input_file mc_input.dat
echo "Done."

echo
set  mc_file = `head -1 mc_input.dat`
echo "Linking input/monte-carlo/$mc_file as $mc_file" 
unlink $mc_file
ln -fs input/monte-carlo/$mc_file
echo "Done"

echo 
echo "Linking $rc_input_file as both rad_corr.dat & paw/macros/rad_corr.dat"
unlink rad_corr.dat
unlink paw/macros/rad_corr.dat
ln -fs $rc_input_file rad_corr.dat
ln -fs ../../$rc_input_file paw/macros/rad_corr.dat
echo "Done."

echo
echo "Deleting log file output/logs/$1.log"
rm output/logs/$1.log
echo "Done."

#echo
#echo "Deleting input/recon-mc/reconmc.in and creating input/recon-mc/reconmc.in"
#unlink input/recon-mc/reconmc.in
#./get_scaler_info.prl $1 > input/recon-mc/reconmc.in
#ln -fs input/recon-mc/kpp_shms_488_recon.in input/recon-mc/reconmc.in
#echo "Done."

echo
echo "Executing ./recon_mc ..."
echo
echo "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
echo
echo "New re-weighted ntuple is: output/mc-ntuples/mc$1.rzdat"
echo
echo "Created new log file: output/logs/$1.log"
echo
./recon_mc<<+ > output/logs/$1.log
$1
