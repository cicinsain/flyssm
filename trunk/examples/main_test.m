clear mex
clear all
close all


nproblem=5;
noptim=1;                                               %Number of optimization per problem

pnames={'ex1', 'ex2', 'ex3', 'ex4','ex5'};                    %Name of all problems

% % or
% pnames={};
% for i=1:4
%     pnames{i}=strcat('ex',num2str(i));
% end


lb={-1*ones(1,2), [0 0], [0 0 0 0 0 0], [0 0 0 0], zeros(1,5)};     %Lower bounds for all problems
ub={ones(1,2), [3 4], [1 1 1 1 16 16], [10 10 10 10], ones(1,5)};  %Upper bounds for all problems

%Specific problem settings

%Problem 1
problem(1).vtr=-1.031628;                               %Value to reach for problem 1

%Problem 2
problem(2).c_L=[-inf -inf];                             %Lower bounds for problem 2 nonlinear inequality constraints
problem(2).c_U=[2 36];                                  %Upper bounds for problem 2 nonlinear inequality constraints

%Problem 3
problem(3).neq=4;                                       %Number of nonlinear equality constraints in problem 3                
problem(3).c_L=-inf;                                    %Lower bounds for problem 3 nonlinear inequality constraints
problem(3).c_U=4;                                       %Upper bounds for problem 3 nonlinear inequality constraints
problem(3).vtr=-0.3888;                                 %Value to reach for problem 3

%Problem 4
problem(4).x_0=[3 4 5 1];                               %Initial point for problem 4
problem(4).int_var=3;                                   %Number of integer variables in problem 4
problem(4).c_L=[-inf -inf -inf ];                       %Lower bounds for problem 4 nonlinear inequality constraints
problem(4).c_U=[8 10 5];                                %Upper bounds for problem 4 nonlinear inequality constraints

%Problem 5
problem(5).vtr=19.8722;                                 %Value to reach for problem 5



%Options for all problems
opts.maxeval=1e5;                                       %Maximum number of function evaluations for all problems
opts.maxtime=3;                                         %Maximum computation time for all problems


%Specific options for some problems
test(3).local.solver='solnp';                           %Specific local solver for problem 3
test(4).local.solver='misqp';                           %Specific local solver for problem 4
test(5).maxtime=100;                                    %Increase the optimization time for problem 5
test(5).log_var=[1:5];                                  %Declare all variables as log_var for problem 5
test(5).local.solver='n2fb';                            %Specific local solver for problem 5

%Extra input parameters for some problems
%Problem 3
k1=0.09755988;
k3=0.0391908;
k2=0.99*k1;
k4=0.9*k3;
param{3}={k1,k2,k3,k4};

%Problem5
%time intervals

t=[0.0 1230.0	3060.0 4920.0 7800.0 10680.0 15030.0 22620.0 36420.0];
 
% Distribution of species concentration
 %	     y(1)    y(2)    y(3)    y(4)    y(5)
 
 yexp=[ 100.0	 0.0	 0.0	 0.0     0.0
    	88.35    7.3     2.3     0.4     1.75
        76.4    15.6     4.5     0.7     2.8
        65.1    23.1     5.3     1.1     5.8
        50.4    32.9     6.0     1.5     9.3
        37.5    42.7     6.0     1.9    12.0
        25.9    49.1     5.9     2.2    17.0
        14.0    57.4     5.1     2.6    21.0
         4.5    63.1     3.8     2.9    25.7 ];

param{5}={t,yexp};
     
         
 
%Call testssm
ssm_test(nproblem,noptim,pnames,lb,ub,problem,opts,test,param);


