# Run the tests using ais++
ais++ -i input-files/input_initial.aisi -o output-files/initial.h5
ais++ -i input-files/input_final.aisi -o output-files/final.h5

# Compare analytical and numerical results with python script
python analysis-test2.py