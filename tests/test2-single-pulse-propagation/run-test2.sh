# Run the tests using ais++
ais++ -i input-files/input_ground_initial.aisi -o output-files/ground_initial.h5
ais++ -i input-files/input_ground_final.aisi -o output-files/ground_final.h5
ais++ -i input-files/input_excited_initial.aisi -o output-files/excited_initial.h5
ais++ -i input-files/input_excited_final.aisi -o output-files/excited_final.h5

# Compare analytical and numerical results with python script
python analysis-test2-ground.py
python analysis-test2-excited.py

# Clean up
rm output-files/*.h5
rm output-files/*.txt