function [ressumgrad,ressumDer] = ressumgrad_robust(pars,y,t_sampled,options,epsilon,a_delta,ressumeps)
%		Finite difference approximation of derivative of residual sum of squares for robust least squares

if nargin<7
    ressumeps=0;
end
a1_=pars(1);
a2_=pars(2);
[~,ressum2,p] = constrImpLS(pars,y,t_sampled,options,-1,epsilon);
[~,ressum1]=constrImpLS_weighted([a1_-a_delta,a2_],y,t_sampled,options,-1,p);
[~,ressum3]=constrImpLS_weighted([a1_+a_delta,a2_],y,t_sampled,options,-1,p);
ressumDer=(ressum3-ressum1)/(2*a_delta);
if ressumDer>0
    ressumDer=NaN;
end
ressumgrad=-(ressum2+ressumeps)/ressumDer;
end