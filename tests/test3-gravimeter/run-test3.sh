# run ais++
ais++ -i input-files/input_initial.aisi -o output-files/data_initial.h5
ais++ -i input-files/input_final.aisi -o output-files/data_final.h5

# run analysis
python analysis_test3.py

# Clean up
rm output-files/*.h5
rm output-files/*.txt