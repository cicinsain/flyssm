clear mex;
clear all;
close all;


% global optimum
%
%x*=[2.23607, 0, 1, 0];
%    
%f(x*)=-40.9575;

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex4';                                          
problem.x_L=[0 0 0 0];                                          
problem.x_U=[10 10 10 10];                                    
problem.x_0=[3 4 5 1];                                               
problem.int_var=3;
problem.c_L=[-inf -inf -inf ];                               
problem.c_U=[8 10 5];                     

opts.local.solver='misqp';
opts.local.finish='misqp';

opts.maxeval=1e7;
opts.maxtime=5;
%========================= END OF PROBLEM SPECIFICATIONS =====================

Results=ssm_kernel(problem,opts);
