# Loop over input files
N=9
#for i in $(seq 0 $N); do
#    # run ais++
#    ais++ -i input-files/input_initial_${i}.aisi -o output-files/data_initial_${i}.h5
#    ais++ -i input-files/input_final_${i}.aisi -o output-files/data_final_${i}.h5
#done
ais++ -i input-files/input_initial_0.aisi -o output-files/data_initial_0.h5
ais++ -i input-files/input_final_0.aisi -o output-files/data_final_0.h5

# run analysis
python analysis_test3.py

# Clean up
rm output-files/*.h5
rm output-files/*.txt
