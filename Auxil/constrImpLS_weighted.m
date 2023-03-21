function [beta,ressum] = constrImpLS_weighted(pars,y,t,options,tk,p)
%		Same as constrImpLS, but weighted least squares
m=length(y);

a1=pars(1);
a2=pars(2);
if nargin<5 || max(tk)==-1
    tk = t(2:end-1);
end
n=numel(tk);

A=[-a1 0; 1 -a2];
B=[1; 0];
C=[0 1];
[Td,Ad] = eig(A);
Tdi = Td^-1;
Z = zeros(m,n+2);
Z(:, 1) = exp(-a2*t);
Z(:, 2) = C*expmtx(Ad, Td, Tdi, t,B);
for i=3:n+2
    tarr = t - tk(i-2);
    firstp = find(tarr>=0,1,'first');
    Z(:, i) = [zeros(firstp-1,1); transpose(C*expmtx(Ad,Td,Tdi,...
        tarr(firstp:end),B))];
end

%[beta,ressum]=lsqlin((Z'*(p.*Z)),(Z'*(p.*y)),[],[],[],[],zeros(n+2,1),...
%    ones(n+2,1)*10000,[],options);
[beta,ressum]=lsqlin(sqrt(p).*Z,sqrt(p).*y,[],[],[],[],zeros(n+2,1),...
    ones(n+2,1)*10000,[],options);
end

