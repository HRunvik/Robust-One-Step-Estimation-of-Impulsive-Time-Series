function [B,fval,d_opt,tk_opt] = mergeimpulsetrain(d,t_sampled,b1,b2)
%		Finds impulse merging configuration that minimizes total sum of
%		impulse weights through linear boolean programming

n=numel(d);
d_merged = (d(1:end-1)+d(2:end).*exp(diff(t_sampled)*b2)).^(-b1/(b2-b1))...
    .*(d(1:end-1)+d(2:end).*exp(diff(t_sampled)*b1)).^(b2/(b2-b1));
d_merged(isnan(d_merged))=0;
d_merged(isinf(d_merged))=0;
D = [d d_merged]';
Aeq = [eye(n) eye(n,n-1)];
Aeq(2:end,n+1:end) = Aeq(2:end,n+1:end) + eye(n-1);
beq=ones(1,n);
intcon = 1:2*n-1;
lb=zeros(1,2*n-1);
ub=ones(1,2*n-1);
opts=optimoptions('intlinprog','Display','off');
[B,fval] = intlinprog(D,intcon,[],[],Aeq,beq,lb,ub,[],opts);
dix_1 = find(B(1:n));
dix_2 = find(B(n+1:end));
dix_2_=[dix_2(1); dix_2(1)+cumsum(diff(dix_2)-1)];
d_opt=zeros(numel(dix_1)+numel(dix_2),1);
d_opt(dix_2_) = D(dix_2+n);
if numel(dix_1)>0
    d_opt(d_opt==0) = D(dix_1(1));
end

tk_opt = -1*ones(2,numel(d_opt));
tk_opt(1,dix_2_)=t_sampled(dix_2);
tk_opt(2,dix_2_)=t_sampled(dix_2+1);
tk_opt(1,tk_opt(2,:)==-1)=t_sampled(dix_1);
end