function [R]= objf_n2fb (N,P,x,NF,R,LTY,TY)
global n_fun_eval 
[f,g,R]= score(x);
n_fun_eval=n_fun_eval+1; 
return
