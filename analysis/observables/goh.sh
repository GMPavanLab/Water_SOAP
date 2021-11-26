
for model in SPC SPCE SPCEb TIP3P TIP3P-FB OPC3 TIP4P TIP4P-2005 TIP4P-FB TIP4P-ICE TIP4P-eps TIP4PEW OPC TIP5P TIP5PEW TIP5P-2018; do 
   cd ${model}/run; 
   rm index.ndx
   echo "a O*" > input_makendx
   echo "a H*" >> input_makendx
   echo "quit" >> input_makendx
   cat input_makendx | gmx_mpi make_ndx -f run.gro -o index.ndx
   rm input_makendx
   echo -e "3\n4\n" | gmx_mpi rdf -f run -s run -n
   cat rdf.xvg | grep -v \# | grep -v @ > ../../plots/rdf_oh/${model}_rdf.dat
   cd ../../; 
done; 
