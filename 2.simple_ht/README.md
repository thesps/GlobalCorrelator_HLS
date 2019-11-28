# First tutorial example

Second tutorial example: an algorithm that computes the scalar sum p<sub>T</sub> of all the objects with |&eta| &lt; 2.4
 * define a simple structure to hold a "particle" object with a p<sub>T</sub> and an &eta; value, stored as integer
 * define a reference implementation
 * try different C++ implementatons for synthesis that result in different performances


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

In this case, the Tcl script is configured to also run the synthesis. 
The report of the synthesis is saved in `proj/solution1/syn/report/algo_main_csynth.rpt`

## Synthesis results (Vivado 2018.3)

For the first implementation, it has a latency of 10 clock cycles, using 568 FFs and 1551 LUTs.
Vivado is not understanding that it can parallelize the sum, apparently because of the if inside the loop

Changing the loop code so that the sum is always executed, but we add zero in some entries, makes HLS understand that the code can be parallelized.
Now the loop takes 1 clock, 82 FFs, 1195 LUTs.

## Changing the clock

We can change the clock to 320 MHz (3.125 ns) from the TCL, recreating the project, or in the GUI from the solution settings (golden gear button in the toolbar) and going in the 'Synthesis' page

After this change, the algorithm now takes 2 clock cycles, 227 FFs, 1195 LUTs: the amount of computation resources used is unchanged, but more registers are needed to hold the result in memory during the processing.

One can go to an even faster clock of 400 MHz (2.5 ns), and you'll see a further increase in FFs to 407.



