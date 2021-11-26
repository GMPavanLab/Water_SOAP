for t in 273 283 293 298 303 313 323 333 343 353 363 373; do
   echo ${t} > ${t}_eps.dat;	
   echo ${t} > ${t}_mu.dat;	
   cd ${t}K; 
   for model in SPC SPCE SPCEb TIP3P TIP3P-FB OPC3 TIP4P TIP4P-2005 TIP4P-FB TIP4P-ICE TIP4P-eps TIP4PEW OPC TIP5P TIP5PEW TIP5P-2018; do 
      cd ${model}/run; 
      rm \#*; 
      echo 0 | gmx_mpi dipoles -s run.tpr -f run.xtc -temp ${t} > log_epsilon  2>&1
      cat log_epsilon | grep "Epsilon =" | awk '{print $3}' >> ../../../${t}_eps.dat;
      cat log_epsilon | grep "Average  =" | awk '{print $3,$7}' >> ../../../${t}_mu.dat;
      cd ../../; 
   done; 
   cd ..;
   cat ${t}_eps.dat | tr '\n' ' ' >>epsilon_T.dat
   cat ${t}_mu.dat | tr '\n' ' ' >>mu_T.dat
   echo "" >>epsilon_T.dat
   echo "" >>mu_T.dat
done
