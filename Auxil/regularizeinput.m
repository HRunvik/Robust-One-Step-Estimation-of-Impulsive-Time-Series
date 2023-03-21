function [ressum,betasum,tkcell] = regularizeinput(a1est,m1val,a2_range,y,t_sampled,a_delta,options,ressumeps)
%		Along gamma curve (a1est+m1val,a2_range), the set of input impulses is regularized through comparison
%		with estimated residual sum of squares based on quadratic approximation from the curve of preliminary
%		estimates (a1est,a2_range).
a1estadj=a1est+m1val;
c2 = zeros(1,numel(a2_range)-2);
for k3=1:numel(a2_range)-2
    ressumG=zeros(1,3);
    a2_ = a2_range(k3+1);
    for k1=1:3
        a1_ = a1est(k3)+(k1-2)*a_delta;
        pars=[a1_ a2_];
        [~,ressumG(k1)] = constrImpLS(pars,y,t_sampled,options);
    end
    dda12 = (ressumG(1) - 2*ressumG(2) + ressumG(3))/a_delta^2;
    c2(k3) = dda12/2;
end
c0 = m1val.^2.*c2;

betasum=zeros(1,numel(a1estadj));
ressum=zeros(1,numel(a1estadj));
tkcell=cell(1,numel(a1estadj));
parfor k1=1:numel(a1estadj)
    if a1estadj(k1)>0
        [ressum(k1),betasum(k1),tkcell{k1}] = sparseinput(a1estadj(k1),a2_range(k1+1),0,c0(k1)+ressumeps,y,t_sampled,options);
    else
        betasum(k1)=NaN;
        ressum(k1)=NaN;
    end
end

end