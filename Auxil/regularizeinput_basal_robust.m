function [ressum,betasum,tkcell] = regularizeinput_basal_robust(a1est,m1val,a2,l_range,y,t_sampled,a_delta,options,ressumeps,epsilon)
%		Same as regularizeinput_basal but with robust least squares.
a1estadj=a1est+m1val;
c2 = zeros(1,numel(l_range)-2);
p = cell(1,numel(l_range)-2);
for k3=1:numel(l_range)-2
    l = l_range(k3+1);
    pars=[a1est(k3) a2];
    [dda12,p{k3}] = d2ressum_robust(pars,y-l,t_sampled,options,epsilon,a_delta);
    c2(k3) = dda12/2;
end
c0 = m1val.^2.*c2;

betasum=zeros(1,numel(a1estadj));
ressum=zeros(1,numel(a1estadj));
tkcell=cell(1,numel(a1estadj));
parfor k1=1:numel(a1estadj)
    if a1estadj(k1)>0
        [ressum(k1),betasum(k1),tkcell{k1}] = sparseinput_robust(a1estadj(k1),a2,l_range(k1),c0(k1)+ressumeps,y,t_sampled,options,p{k1});
    else
        betasum(k1)=NaN;
        ressum(k1)=NaN;
    end
end

end