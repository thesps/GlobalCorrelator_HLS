
# Structure of this example

* `run_hls.tcl`: configuration & startup file for Vivado HLS defining the project, the input files, etc
* `src` directory with the header file and implementation for the synthesis (i.e. to to be compiled into firmware)
   * `data.h`: define the dataformats
   * `algo.h`, `algo.cpp`: header and source file for the code to be synthethised
* `algo_ref.cpp`: a clean and readable reference implementation of the algorithm, for validation
* `algo_test.cpp`: C++ testbench, that compares the implementation to be synthetised to the reference one to make sure they're bitwise identical (i.e. that any re-writing of the C++ code to make HLS work better didn't change the results)

# Running the example
## Batch mode using Tcl script
`vivado_hls -f run_hls.tcl`

## To open the project in the GUI
After creating the project using the tcl script (default name from script is `proj`)

`vivado_hls -p proj`
