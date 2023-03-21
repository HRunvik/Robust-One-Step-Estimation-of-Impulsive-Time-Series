function [a1est,m1val] = gammaCurve_robust(a1_range,a2_range,y,t_sampled,options,bigstep,epsilon,ressumeps)
%		Same as gammaCurve, but with robust least squares
if nargin<8
    ressumeps=0;
end
a_delta=a1_range(2)-a1_range(1);
a_delta_small = a_delta*0.1;

a1est=zeros(1,numel(a2_range));
m1val=zeros(1,numel(a2_range));
for k2=1:numel(a2_range)
    a2_=a2_range(k2);
    if mod(k2,bigstep)==1
        if k2>1
            ressumDer=zeros(1,numel(a1grid));
            ressumgrad=zeros(1,numel(a1grid));
            for k1=1:numel(a1grid)
                a1_=a1grid(k1);
                pars=[a1_,a2_];
                [ressumgrad(k1),ressumDer(k1)]=ressumgrad_robust(pars,y,t_sampled,options,epsilon,a_delta_small,ressumeps);
            end
            [m1val(k2),k1r]=min(ressumgrad);
            a1est_=a1grid(k1r);
        end
        ressumDer=zeros(1,numel(a1_range));
        ressumgrad=zeros(1,numel(a1_range));
        parfor k1=1:numel(a1_range)
            a1_=a1_range(k1);
            pars=[a1_,a2_];
            [ressumgrad(k1),ressumDer(k1)]=ressumgrad_robust(pars,y,t_sampled,options,epsilon,a_delta_small,ressumeps);
        end
        A = find(ressumDer>0);
        ressumgrad(A:end) = NaN;
        [m1val(k2),k1r]=min(ressumgrad);

        a1est(k2)=a1_range(k1r);
        a1e = max(a1est(k2)-4*a_delta,a1_range(1))+4*a_delta;
        a1e = min(a1e+4*a_delta,a1_range(end))-4*a_delta;
        a1grid = a1e-4*a_delta:a_delta:a1e+4*a_delta;
        k2_=k2;
        if k2>1 && abs(a1est_-a1est(k2))>a_delta/2
            runback=true;
            a1grid_b=a1grid;
        else
            runback=false;
        end
        while runback
            k2_=k2_-1;
            a2_=a2_range(k2_);
            ressumDer=zeros(1,numel(a1grid));
            ressumgrad=zeros(1,numel(a1grid));
            for k1=1:numel(a1grid)
                a1_=a1grid(k1);
                pars=[a1_,a2_];
                [ressumgrad(k1),ressumDer(k1)]=ressumgrad_robust(pars,y,t_sampled,options,epsilon,a_delta_small,ressumeps);
            end
            [m1val_,k1r]=min(ressumgrad);
            if m1val_<m1val(k2_)&&k2_>k2-bigstep
                a1est(k2_)=a1grid(k1r);
                m1val(k2_)=m1val_;
                a1e = max(a1est(k2_)-4*a_delta,a1_range(1))+4*a_delta;
                a1e = min(a1e+4*a_delta,a1_range(end))-4*a_delta;
                a1grid = a1e-4*a_delta:a_delta:a1e+4*a_delta;
            else
                runback=false;
                a1grid=a1grid_b;
            end
        end
    else
        ressumDer=zeros(1,numel(a1grid));
        ressumgrad=zeros(1,numel(a1grid));
        for k1=1:numel(a1grid)
            a1_=a1grid(k1);
            pars=[a1_,a2_];
            [ressumgrad(k1),ressumDer(k1)]=ressumgrad_robust(pars,y,t_sampled,options,epsilon,a_delta_small,ressumeps);
        end
        [m1val(k2),k1r]=min(ressumgrad);
        a1est(k2)=a1grid(k1r);
        if mod(k2,bigstep)==2
            a1e = max(a1est(k2)-4*a_delta,a1_range(1))+4*a_delta;
            a1e = min(a1e+4*a_delta,a1_range(end))-4*a_delta;
        else
            a1e = 2*a1est(k2)-a1est(k2-1);
            a1e = max(a1e-4*a_delta,a1_range(1))+4*a_delta;
            a1e = min(a1e+4*a_delta,a1_range(end))-4*a_delta;
        end
        a1grid = a1e-4*a_delta:a_delta:a1e+4*a_delta;
    end

end
end