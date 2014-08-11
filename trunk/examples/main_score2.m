clear mex;
clear all;
close all;

%========================= PROBLEM SPECIFICATIONS ===========================
small_value = 1E-32;

problem.f='score2';           
%problem.vtr = 0;

%initial values


    %promoter_strengths:
problem.x_0 = [  12.48991375   26.64956536  24.21400803  19.14696636 ... 
    ...%genetic_interconnect_matrix:
                 -0.02875641   0.03773355  -0.08696411   0.02085833 ...  
                  0.07752513  -0.05127861   0.08351538  -0.03693460 ...  
                 -0.02272315  -0.05630814   0.03840626   0.07150631 ...  
                  0.02276078  -0.02350780   0.01856495  -0.06116937 ...
    ...%external_input_strengths:
                 -0.14654266   0.01651517   0.06442419  -0.05467585 ...
                 -0.12247616  -0.10550997   0.10687539   ...
                  0.02896215  -0.07365022  -0.06193835   ...
                  0.02710513  -0.04491538   0.10493181   ... 
    ...%maternal_connection_strengths:
                  ...
    ...%promoter_thresholds:
                  ...
    ...%diffusion_parameter(s):
                  ...
    ...%protein_half_lives:
                 13.76469068    7.27890037  11.63317492  12.03105457 ...           
    ...%translational_transcriptional_delays:
                  ];

%disp(problem.x_0)             
              
%lower_bounds:              
problem.x_L = [ 10.0    10.0    10.0    10.0 ...
...%genetic_interconnect_matrix:
           -5    -5    -5    -5  ... %N/A was -0.6
           -5    -5    -5    -5  ...
           -5    -5    -5    -5  ...
           -5    -5    -5    -5  ...
...%external_input_strengths:
           -5    -5    -5    -5  ...
           -5    -5    -5  ...
           -5    -5    -5  ...
           -5    -5    -5  ...
...%maternal_connection_strengths:
                 ...
...%promoter_thresholds:
                 ...
...%diffusion_parameter(s):
                 ...
...%protein_half_lives:
                    5       5       5       5 ...
...%translational_transcriptional_delays:
                 ];


%upper_bounds: 
problem.x_U = [ 30.0    30.0    30.0    30.0 ...
...%genetic_interconnect_matrix:
                5 5 5 5 ... %this should be double(intmax) %n/a WAS +0.6
                5 5 5 5 ...
                5 5 5 5 ...
                5 5 5 5 ...
...%external_input_strengths:
                5 5 5 5 ...
                5 5 5 ...
                5 5 5 ...
                5 5 5 ...
...%maternal_connection_strengths: 
                 ...
...%promoter_thresholds:
                 ...
...%diffusion_parameter(s):
                 ...
...%protein_half_lives:
                20      20      20      20 ...
...%translational_transcriptional_delays:
                ];
          
            
opts.maxtime=60000;        
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
