This package contains the source files to perform an optimization for the study 
"High-fidelity Musculoskeletal Modeling Reveals a Motor Planning Contribution to the Speed-Accuracy Tradeoff".

The code was developed on a Unix environment with OpenSim 3.3. 
I included a CMakeLists.txt file as an example, but the user should follow the instructions on how to compile
example C++ code from the OpenSim website (e.g., https://simtk-confluence.stanford.edu/display/OpenSim33/API+Examples), and use
the CMakeLists.txt from there.

The main file is ArmController.cpp.
When performing an optimization, the user must specify in the main file the following:
the initial coordinates of the arm (initial_coords), the target hand position (end_eff_target),
and the duration of the optimization (topt_full_dim.set_opt_time).  In general, the duration should increase for targets that are "further away".
However, increasing the duration makes the optimization slower and more prone to fall in local minima.  

The parameters maxIter and lambda determine the number of iterations and the population size in the Covariance Matrix Adaptation optimization.
Our default parameters (650 and 32, respectively) takes about 4 hours when running our code in parallel with 16 threads.  

The cost function is specified in the file OptSystem.cpp (function: evaluate_cost_function).   
The default value for the variable task_term_accuracy is set to 0.0025, which means that the final hand position
is only encouraged to be within 5 cm of the target (i.e., sqrt(0.0025) meters).  To encourage finding a closer solution to the target,
the term should be decreased.  

Our code uses the model "ue_rigid.osim" from the study "Benchmarking of dynamic simulation predictions in two software platforms using an upper limb musculoskeletal model". 
 
When the optmization is done, it produces a file "ue_full_bestSoFar_states.sto".  The movement can be visualized and analyzed in
the OpenSim 3.3 GUI.

  

