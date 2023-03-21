% Parameter and input estimation for second order system
% Experiments with estimates with large errors

n_impulses=3; % number of impulses
dmin=2; % lowest impulse amplitude
dmax = 7; % Impulse amplitudes uniformly distributed between dmin
% and dmax

Tmax = 5;
Tmin = 2; % Distance between impulses uniformly distributed between
% Tmin and timescale

min_a1=0.4; % lowest a1 value. a1 uniformly distributed between min_a1 and
% min_a1+1
min_adiff = 4; % lowest difference between a1 and a2. a2 uniformly
% distributed between a1+min_adiff and a1+min_adiff+1

delta_t_sampled=0.5; % sampling time

mnoisestd=0.008; % noise std

plotandpause = false; % look at each run separately
options = optimset('lsqlin');
options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);

a2vals = [5.030457470200914 7.144203939055338];
for k=[26 43]
    rng((k-1)*5)
    a1= rand(1)*0.6+min_a1;
    a2= a1+min_adiff+rand(1);

    [y_nonoise,t_sampled,x,t,dSeq] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,0.02,delta_t_sampled);
    gaussiannoise = mnoisestd*randn(size(t_sampled));
    y = y_nonoise+gaussiannoise';

    tk = t_sampled(2:end-1);
    m = length(y);
    n = length(tk);

    a1_min=a1*0.5;
    a1_max=a1*2;
    a1_delta=0.002;
    a1_range=a1_min:a1_delta:a1_max;

    % Basal level
    lmin = -0.15;
    lmax = 0.08;
    level = linspace(lmin,lmax,116);


    a2val = a2vals((k>35)+1);
    [a1est,m1val]=gammaCurve_basal(a1_range,a2val,level,y,t_sampled,options,15,mnoisestd^2);
    a_delta_ = 0.001;
    [ressumb,betasumb,tk]=regularizeinput_basal(a1est,m1val,a2val,level,y,t_sampled,a_delta_,options,mnoisestd^2);

    bicAll = 2*betasumb*log(m) + m*log(ressumb);
    bicAll(betasumb>m/4) = NaN;
    [~,minbicix] = min(bicAll);
    lval=level(minbicix);
    a1val=a1est(minbicix)+m1val(minbicix);

    A=[-a1val 0; 1 -a2val];
    [beta,ressumest] = constrImpLS([a1val a2val],y-lval,t_sampled,options,tk{minbicix});

    dSeqest = [tk{minbicix};beta(3:end)'];
    dSeqest = mergeImpulses(dSeqest,a1val,a2val,t_sampled);
    xest = xSol(A,dSeqest,t,[beta(2);beta(1)]);
    %%

    figure(1)
    subplot(2,1,1)
    plot(t,x(2,:),'k')
    hold on
    plot(t_sampled,y,'bo')
    plot(t,xest(2,:)+lval,'r--')
    yline(0,'-.')
    yline(lval,'r-.')
    if k<35
        xlim([0 10])
        ylim([-0.2 1])
    else
        ylim([-0.2 1.68])
        xlim([0 12])
    end
    hold off
    ylabel('y','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    legend({'Synthetic data' 'Gaussian noise', '$\hat y$', '$y_0$', '$\hat y_0$'}, 'interpreter','latex','Location','north','Orientation','horizontal')

    subplot(2,1,2)
    stem(dSeq(1,:),dSeq(2,:),'k')
    hold on
    stem(dSeqest(1,:),dSeqest(2,:),'r--')
    hold off
    ylabel('Impulse weight','Interpreter','latex')
    xlabel('Time','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')


    figure(2)
    plot(level,a1est+m1val)
    hold on
    plot(0,a1,'r*')
    hold off

    %pause
end
