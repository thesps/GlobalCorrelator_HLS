# get the configuration
source config_hls_fullpfalgo_mp7.tcl

# open the project, don't forget to reset
open_project -reset "l1pfpuppi-ku15p-II2-nomux-csim"
set_top ${l1pfTopFunc}
#set_top pfalgo3_full
add_files firmware/simple_fullpfalgo.cpp -cflags "-DTESTKU15P -DHLS_pipeline_II=2"
add_files puppi/firmware/simple_puppi.cpp -cflags "-DTESTKU15P -DHLS_pipeline_II=2"
add_files -tb fullpfalgo_test.cpp  -cflags "-DTESTKU15P -DHLS_pipeline_II=2 -DMP7_TOP_FUNC=${l1pfTopFunc} -DMP7_REF_FUNC=${l1pfRefFunc} -DMP7_VALIDATE=${l1pfValidate}"
add_files -tb simple_fullpfalgo_ref.cpp -cflags "-DTESTKU15P"
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTKU15P"
add_files -tb utils/test_utils.cpp -cflags "-DTESTKU15P"
add_files -tb DiscretePFInputs.h    -cflags "-DTESTKU15P -std=c++0x"
add_files -tb utils/DiscretePFInputs_IO.h -cflags "-DTESTKU15P -std=c++0x"
add_files -tb data/regions_TTbar_PU140.dump
add_files -tb puppi/simple_puppi_ref.cpp -cflags "-DTESTKU15P"
add_files -tb vertexing/simple_vtx_ref.cpp -cflags "-DTESTKU15P"

# reset the solution
open_solution -reset "solution"
#set_part {xc7vx690tffg1927-2}
#set_part {xcvu9p-flgb2104-2-i}
set_part {xcku15p-ffva1760-2-e}
#set_part {xcku115-flvd1517-2-i}
#create_clock -period 3.125 -name default
create_clock -period 2.7 -name default
set_clock_uncertainty 1.5

config_interface -trim_dangling_port
# do stuff
csim_design
#csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog -vendor "cern-cms" -version ${l1pfIPVersion} -description "${l1pfTopFunc}"

# exit Vivado HLS
exit
