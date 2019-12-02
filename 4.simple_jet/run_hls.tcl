# open the project, don't forget to reset
open_project -reset proj
set_top algo_main
add_files src/algo.cpp
add_files -tb algo_test.cpp 
add_files -tb algo_ref.cpp
add_files -tb pfcands_ttbar.txt

# reset the solution
open_solution -reset "solution1"
###  MP7 (Virtex-7 690T)
# set_part {xc7vx690tffg1927-2}
##   VCU118 dev kit (VU9P), 320 MHz
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 3.125

# just check that the C++ compiles
csim_design

# synthethize the algorithm
csynth_design

exit
