function [ressum,impsum,tk] = sparseinput(a1,a2,l,c0,y,t_sampled,options)
%		Computes sparse input by removing impulses until residual sum is
%		approximately equal to c0

[betab0,resb0]=constrImpLS([a1 a2],y-l,t_sampled,options);
[~,~,d_opt,tk_opt]=mergeimpulsetrain(betab0(3:end)',t_sampled(2:end-1),a1,a2);

I0 = numel(t_sampled);
resb = resb0;
tk=NaN;
if resb>=c0
    impsum = NaN;
    ressum = NaN;
    return
end
while resb<c0
    resbold = resb;
    tkold=tk;
    I0=I0-1;
    if I0==0
        impsum=NaN;
        ressum=NaN;
        return
    end
    [~,I] = maxk(d_opt,I0);
    tk=[tk_opt(1,I) tk_opt(2,I)];
    tk(tk==-1) = [];
    tk=sort(tk);
    %tk = [0 tk];
    [~,resb]=constrImpLS([a1 a2],y-l,t_sampled,options,tk);
end
alpha = (c0-resbold)/(resb-resbold);
%impsum = I0 + 1-alpha;
impsum=I0+round(1-alpha);
if alpha<0.5
    tk=tkold;
end
%ressum = [resbold;resb];
ressum=resbold*(impsum-I0)+resb*(1-(impsum-I0));


