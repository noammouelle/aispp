# Loop over input files
N=9
for i in $(seq 0 $N); do
    # run ais++
    ais++ -i input-files/input_initial_${i}.aisi -o output-files/data_initial_${i}.h5
    ais++ -i input-files/input_final_${i}.aisi -o output-files/data_final_${i}.h5
done

# run analysis
python analysis_test4.py

# Clean up
rm output-files/*.h5
rm output-files/*.txt
