# Run the tests using ais++
ais++ -i input-files/input_zero_pot_initial.aisi -o output-files/zero_pot_initial.h5
ais++ -i input-files/input_zero_pot_final.aisi -o output-files/zero_pot_final.h5
ais++ -i input-files/input_linear_pot_initial.aisi -o output-files/linear_pot_initial.h5
ais++ -i input-files/input_linear_pot_final.aisi -o output-files/linear_pot_final.h5
ais++ -i input-files/input_quadratic_pot_initial.aisi -o output-files/quadratic_pot_initial.h5
ais++ -i input-files/input_quadratic_pot_final.aisi -o output-files/quadratic_pot_final.h5

# Compare analytical and numerical results with python script
python analysis-test1.py