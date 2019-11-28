# First tutorial example

An algorithm that computes particle-level jets


# Structure of this example

* `run_hls.tcl`: configuration & startup file for Vivado HLS defining the project, the input files, etc
* `src` directory with the header file and implementation for the synthesis (i.e. to to be compiled into firmware)
   * `data.h`: define the dataformats
   * `algo.h`, `algo.cpp`: header and source file for the code to be synthethised
* `algo_ref.cpp`: a clean and readable reference implementation of the algorithm, for validation
* `algo_test.cpp`: C++ testbench, that runs the algorithm. 
* `pfcands_ttbar.txt`: an input file of Puppi candidates and ak4 Puppi jets dumped from CMSSW, used as input for the test bench


# Running the example
## Batch mode using Tcl script
`vivado_hls -f run_hls.tcl`

## Synthesis results (Vivado 2018.3)

For the first implementation at II=9, latency 156 clock cycles, 356 DSPs, 78k FFs, 67k LUTs.

## Changing the clock

We can change the clock to 320 MHz (3.125 ns) from the TCL, recreating the project, or in the GUI from the solution settings (golden gear button in the toolbar) and going in the 'Synthesis' page

After this change, the algorithm now takes 2 clock cycles, 227 FFs, 1195 LUTs: the amount of computation resources used is unchanged, but more registers are needed to hold the result in memory during the processing.

One can go to an even faster clock of 400 MHz (2.5 ns), and you'll see a further increase in FFs to 407.



