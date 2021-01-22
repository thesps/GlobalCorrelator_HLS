# open the project, don't forget to reset
open_project -reset JetLoop
set_top jet_loop
add_files seedcone.cpp -cflags -std=c++0x
add_files -tb algo_test.cpp 
add_files -tb algo_ref.cpp
add_files -tb pfcands_ttbar.txt

# reset the solution
open_solution -reset "solution1"
###  MP7 (Virtex-7 690T)
# set_part {xc7vx690tffg1927-2}
##   VCU118 dev kit (VU9P), 320 MHz
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 2.5

# just check that the C++ compiles
#csim_design

# synthethize the algorithm
csynth_design

exit
