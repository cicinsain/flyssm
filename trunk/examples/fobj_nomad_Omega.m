function [A,l,u]= fobj_nomad_Omega (n)
A = [eye(n)];
l= [ 10; 10; 10; 10; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -1; -0.6; -0.6; -0.6; -0.6; -10; -10; -10; -10; 0; 0; 0; 0; 5; 5; 5; 5; 2; 2; 2; 2; ];
u= [ 30; 30; 30; 30; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0.6; 0.6; 0.6; 0.6; 10; 10; 10; 10; 0.3; 0.3; 0.3; 0.3; 20; 20; 20; 20; 10; 10; 10; 10; ];
return
