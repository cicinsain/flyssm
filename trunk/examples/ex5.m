%
%-------------------------------------------------------------------
% 						ISOMERIZACION del ALPHA-PINENO
%
% Problema obtenido de COPS(1999)
%-------------------------------------------------------------------
% Parametros a determinar: Coeficientes de reaccion (theta)
% Variables:
%					Concentraciones de las especies (x)
% 					Variacion de la concentracion con el tiempo (y)
%-------------------------------------------------------------------
% Variables: x=theta (model parameters)
%
% CGM 16-05-00


function [f,gg,R] = ex5(x,t,yexp)


%Inicializaci�n necesaria en Matlab 5.0
	f=[];
   

% INTEGRACION MEDIANTE LSODES
% Vectores de entrada para el M-file

   vx=[100. 0. 0. 0. 0.];
   nvar=5;
   xinfo(1)=x(1);
   xinfo(2)=x(2);
   xinfo(3)=x(3);
   xinfo(4)=x(4);
   xinfo(5)=x(5);
   
   
[f,xinfoconstf,xinfopath,xinfomatrix, yt] = pineno(vx,nvar,xinfo); 

if (xinfopath(1)< 0.0),
   disp('error in istate:')
   xinfopath
   pause 
end

%fprintf('istate:%d\n',xinfopath(1))
%fprintf('n_steps(5000):%d\n',xinfopath(2))
%fprintf('rwork(2000):%d\n',xinfopath(3))
%fprintf('iwork(1000):%d\n',xinfopath(4))
%fprintf('\n')
%pause

yteor(1,1:5)=vx;
yteor(2:9,1:5)=xinfomatrix(1:8,1:5);
%yteor
%pause
   
% Objective function:   
%x
R=(yteor-yexp);
R=reshape(R,numel(R),1);



f = sum(sum((yteor-yexp).^2));
gg=0;



   
     
   
 
 
 