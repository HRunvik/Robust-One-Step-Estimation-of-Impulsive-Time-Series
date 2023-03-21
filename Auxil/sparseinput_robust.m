function [ressum,impsum,tk] = sparseinput_robust(a1,a2,l,c0,y,t_sampled,options,p)
%		Same as sparseinput but with robust least squares

[betab0,resb0]=constrImpLS_weighted([a1 a2],y-l,t_sampled,options,-1,p);
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
    [~,resb]=constrImpLS_weighted([a1 a2],y-l,t_sampled,options,tk,p);
end
alpha = (c0-resbold)/(resb-resbold);
%impsum = I0 + 1-alpha;
impsum=I0+round(1-alpha);
if alpha<0.5
    tk=tkold;
end
%ressum = [resbold;resb];
ressum=resbold*(impsum-I0)+resb*(1-(impsum-I0));


