
## Compilation

There are three different executable files, one for each optimization algorithm implemented. 

Compile the with `METHOD=-DAMOSA` to get `fly_amosa` executable which optimizes the problem using **Archived Multi-objective Simulated Annealing** algorithm.

`make veryclean; make METHOD=-DAMOSA deps; make METHOD=-DAMOSA`

Compilation with the flag `METHOD=-DNSGA2` produces the `fly_nsga2` executable which optimize the problem using **Non-dominated Sorting Genetic Algorithm**. 

`make veryclean; make METHOD=-DNSGA2 deps; make METHOD=-DNSGA2`

And finally compilation with the `METHOD=-DSS` generated the `fly_ss` executable which optimize the problem using **Scatter Search** algorithm.

`make veryclean; make METHOD=-DSS deps; make METHOD=-DSS`

## Input file

Almost all the information needed for the simulator and the optimizer are included in the input file. Three different sections: `$nsga`, `$amosa`, and `$ss` are containing optimization parameters for NSGA2, AMOSA, and SS respectively.

## Execution

There are several options available via command line to run the optimizer. Almost all the command line parameters are also available on the input file as well. In order to get the list of parameters run: `./fly_[yourmethod] -h` to get:

    Usage: flyMOP [options] <datafile>

    Argument:
      <datafile>          input data file

    Options:
      -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers
      -b <bkup_freq>      write state file every <bkup_freq> * tau moves
      -B                  run in benchmark mode (only do fixed initial steps)
      -D                  debugging mode, prints all kinds of debugging info
      -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>
      -E                  run in equilibration mode
      -f <param_prec>     float precision of parameters is <param_prec>
      -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh
      -h                  prints this help message
      -i <stepsize>       sets ODE solver stepsize (in minutes)
      -l                  echo log to the terminal
      -m <score_method>   w = wls, o=ols score calculation method
      -n                  nofile: don't print .log or .state files
      -N                  generates landscape to .landscape file in equilibrate mode 
      -o                  use oldstyle cell division times (3 div only)
      -p                  prints move acceptance stats to .prolix file
      -Q                  quenchit mode, T is lowered immediately to zero
      -s <solver>         choose ODE solver
      -t                  write timing information to .times file
      -v                  print version and compilation date
      -w <out_file>       write output to <out_file> instead of <datafile>
      -y <log_freq>       write log every <log_freq> * tau moves

    Please report bugs to <yoginho@usa.net>. Thank you!

Sample run command would be like:

`./fly/fly_ss -s rck -t -i 4.0 -a 0.001 input/sample_input.inp`

