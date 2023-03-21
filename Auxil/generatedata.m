function [y,t_sampled,x,t,dSeq] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,delta_t_data,delta_t_sampled,firstorder)
%		Synthetic data generation function

if nargin<10
    firstorder=false;
end
A=[-a1 0; 1 -a2];
impulsediffs = (Tmax-Tmin)*rand(1,n_impulses)+Tmin;
impulsetimes = [0 cumsum(impulsediffs)];
impulseweights = (dmax-dmin)*rand(1,n_impulses+1)+dmin;
dSeq = [impulsetimes; impulseweights];
t_end=dSeq(1,end)+rand(1,1)*(Tmax-Tmin)+Tmin;
t=0:delta_t_data:t_end;
x = xSol(A, dSeq,t,[0;0]);
Tstart = impulsetimes(2)/2;%rand(1,1)*(Tmax-Tmin)+Tmin;
t_sampled=round(Tstart/delta_t_sampled)*delta_t_sampled:delta_t_sampled:t(end);
if firstorder
    x_sampled=interp1(t,x(1,:),t_sampled);
else
    x_sampled=interp1(t,x(2,:),t_sampled);
end

dSeq = dSeq(:,2:end);
dSeq(1,:)=dSeq(1,:)-t_sampled(1);
y = x_sampled';
t = t-t_sampled(1);
t_sampled=t_sampled-t_sampled(1);
x(:,t<0)=[];
t(t<0)=[];
end