function [dda12,p] = d2ressum_robust(pars,y,t_sampled,options,epsilon,a_delta)
%		Second derivative approximation for residual sum of squares

a1_=pars(1);
a2_=pars(2);
[~,ressum2,p] = constrImpLS(pars,y,t_sampled,options,-1,epsilon);
[~,ressum1]=constrImpLS_weighted([a1_-a_delta,a2_],y,t_sampled,options,-1,p);
[~,ressum3]=constrImpLS_weighted([a1_+a_delta,a2_],y,t_sampled,options,-1,p);
dda12=(ressum1 - 2*ressum2 + ressum3)/a_delta^2;