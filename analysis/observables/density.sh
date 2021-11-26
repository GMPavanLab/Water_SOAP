
rm density_T.dat
for t in 273 283 293 298 303 313 323 333 343 353 363 373; do
   echo ${t} > ${t}_dens.dat;	
   cd ${t}K; 
   for model in SPC SPCE SPCEb TIP3P TIP3P-FB OPC3 TIP4P TIP4P-2005 TIP4P-FB TIP4P-ICE TIP4P-eps TIP4PEW OPC TIP5P TIP5PEW TIP5P-2018; do 
      cd ${model}/run; 
      plumed driver --mf_xtc run.xtc --plumed ../../../plumed-density.dat; 
      rm bck.*; 
      rm avg_std_density.dat; 
      python3 ../../../avg_std.py density.dat >> avg_std_density.dat; 
      cat avg_std_density.dat >> ../../../${t}_dens.dat;
      cd ../../; 
   done; 
   cd ..;
   cat ${t}_dens.dat | tr '\n' ' ' >>density_T.dat
   echo "" >>density_T.dat
done
