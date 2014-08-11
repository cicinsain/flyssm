clear mex;
clear all;
close all;

%========================= PROBLEM SPECIFICATIONS ===========================
small_value = 1E-32;

problem.f='score';           
%problem.vtr = 0;

%initial values


    %promoter_strengths:
problem.x_0 = [  25.19644112   24.11349283  24.84952218  25.64152879 ... 
    ...%genetic_interconnect_matrix:
                 -0.01392977   -0.01937628  -0.01572716  -0.01162660 ...  
                 -0.01339641    0.01874027   0.03135462   0.02433955 ...  
                 -0.00763001    0.01777379   0.02498920   0.01621901 ...  
                  0.00083402    0.00136981   0.02535139  -0.02114692 ...
    ...%external_input_strengths:
                  0.02471927   0.03145242   0.01257521  -0.04127696 ... 
                  0.01640896  -0.05376330  -0.03293237   0.00000000 ... 
                  0.03466005  -0.04684228  -0.03388611   0.00000000 ...  
                 -0.02398744   0.03099123   0.04594317   0.00000000 ... 
    ...%maternal_connection_strengths:
                  0.000000e+00  0.000000e+00 0.000000e+00 0.000000e+00 ...
    ...%promoter_thresholds:
                 -2.500000e+00 -2.500000e+00 -2.500000e+00 -2.500000e+00 ...
    ...%diffusion_parameter(s):
                  0.23700000    0.29999999   0.29999999   0.11500000 ...
    ...%protein_half_lives:
                 14.91865154    5.19861693   6.29368303   5.87602929 ...           
    ...%translational_transcriptional_delays:
                  6.00000000    6.00000000   6.00000000   6.00000000];

              
%lower_bounds:              
problem.x_L = [ 10.0    10.0    10.0    10.0 ...
...%genetic_interconnect_matrix:
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ... %N/A was -0.6
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
...%external_input_strengths:
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
           -double(intmax)    -double(intmax)    -double(intmax)    -double(intmax)  ...
...%maternal_connection_strengths:
                 0.0   0.0   0.0   0.0  ...
...%promoter_thresholds:
                    2    2    2    2  ...
...%diffusion_parameter(s):
                    0       0       0       0 ...
...%protein_half_lives:
                    5       5       5       5 ...
...%translational_transcriptional_delays:
                    2       2       2       2];


%upper_bounds: 
problem.x_U = [ 30.0    30.0    30.0    30.0 ...
...%genetic_interconnect_matrix:
                double(intmax) double(intmax) double(intmax) double(intmax) ... %this should be double(intmax) %n/a WAS +0.6
                double(intmax) double(intmax) double(intmax) double(intmax) ...
                double(intmax) double(intmax) double(intmax) double(intmax) ...
                double(intmax) double(intmax) double(intmax) double(intmax) ...
...%external_input_strengths:
                double(intmax) double(intmax) double(intmax) double(intmax) ...
                double(intmax) double(intmax) double(intmax) double(intmax) ...
                double(intmax) double(intmax) double(intmax) double(intmax) ...
                double(intmax) double(intmax) double(intmax) double(intmax) ...
...%maternal_connection_strengths: 
                0.00000001    0.00000001   0.00000001   0.00000001 ...
...%promoter_thresholds:
                3 3 3 3 ...
...%diffusion_parameter(s):
                0.30000     0.30000     0.30000     0.30000 ...
...%protein_half_lives:
                20      20      20      20 ...
...%translational_transcriptional_delays:
                10      10      10      10];
          
            
opts.maxtime=120000;        
%opts.maxtime=600;        
opts.maxeval=0.5e7; 
%opts.ndiverse=25; %for multistart
%opts.diverse_criteria=2;
%opts.local.bestx=1;
%number of func evals before applying LS for the 1st time
%opts.local.n1=300;
%number of func evals between two LSs 
%opts.local.n2=0;
%opts.weight=1e9;

%this doesn't work with solnp solver - outputs some error that sounds
% - dutch?! Update: it's not dutch; it's just that the error message should
%be written in two rows, but he prints it columnwise, character by
%character. So if you read first the characters in odd position (first row of 
%the error message), and than the ones in even position (second row of the 
%error message), you will get the error message transleted from 'dutch' to 
%english.

%opts.log_var=[1:27, 29:31, 33:35, 49:52];
%opts.local.solver='solnp'; %works fine
%opts.local.solver='dhc'; %Warning: The P-code file dhc.p was generated prior to MATLAB version 7.5 (R2007b) and will not be supported in a future release.
opts.local.solver='fminsearch';%works fine
%opts.local.solver='dn2fb'; % undefined variable TY 
%opts.local.solver='n2fb'; % undefined variable TY 
%opts.local.solver='lsqnonlin'; %The LevenbergMarquardt option is no longer valid. Set the Algorithm option instead.

%opts.local.solver=0;
%opts.local.finish='n2fb';
%opts.local.finish='fminsearch';
%opts.local.finish='dhc';
%opts.local.finish='nomadm'; %best one for n2fb %needs optimization toolbox
opts.local.finish='solnp'; %the only one that works fine


%========================= END OF PROBLEM SPECIFICATIONS =====================
%time intervals

 
    

Results=ssm_kernel(problem,opts);
