% Parameter and input estimation for second order system
% Monte Carlo runs with synthetic data with added noise
% No outliers, unknown basal level

n_runs=200; % number of runs

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

options = optimset('lsqlin');
options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);

truea=zeros(n_runs,2);
levelestvect=zeros(1,n_runs);
b1estvect=zeros(1,n_runs);
b2estvect=zeros(1,n_runs);
miny=zeros(1,n_runs);

for k=1:n_runs
    disp(['Run ' num2str(k)])
    rng((k-1)*5)
    a1= rand(1)*0.6+min_a1;
    a2= a1+min_adiff+rand(1);
    truea(k,:)=[a1 a2];
    
    [y_nonoise,t_sampled,x,t] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,0.02,delta_t_sampled);
    gaussiannoise = mnoisestd*randn(size(t_sampled));
    y = y_nonoise+gaussiannoise';
    miny(k)=min(y);
    
    tk = t_sampled(2:end-1);
    m = length(y);
    n = length(tk);

    a1_min=a1*0.5;
    a1_max=a1*2;%(a1+a2)/2;
    a1_delta=0.002;
    a2_min=a1_max+a1_delta;
    a2_max=a2*1.5;
    a1_range=a1_min:a1_delta:a1_max;
    a2_range=linspace(a2_min,a2_max,20); 

    %% Basal level
    lmin = -0.15;
    lmax = 0.08;
    level = linspace(lmin,lmax,116);

    levelest=zeros(size(a2_range));
    bictest=zeros(size(a2_range));
    a1estvect=zeros(size(a2_range));
    tic
    for k1=1:numel(a2_range)
        a2val = a2_range(k1);
        [a1est,m1val]=gammaCurve_basal(a1_range,a2val,level,y,t_sampled,options,15,mnoisestd^2);
        
        a_delta_ = 0.001;
        [ressumb,betasumb]=regularizeinput_basal(a1est,m1val,a2val,level,y,t_sampled,a_delta_,options,mnoisestd^2);

        bicAll = 2*betasumb*log(m) + m*log(ressumb);
        bicAll(betasumb>m/4) = NaN;
        [bictest(k1),minbicix] = min(bicAll);
        levelest(k1)=level(minbicix);
        a1estvect(k1) = a1est(minbicix)+m1val(minbicix);
    end
    toc
    [~,minbicix]=min(bictest);
    levelestvect(k)=levelest(minbicix);
    b1estvect(k)=a1estvect(minbicix);
    b2estvect(k)=a2_range(minbicix);

end
%%
figure(1)
histogram(levelestvect)
hold on
histogram(miny,16)
hold off
legend('One-step estimation','$\min y_k$','Interpreter','Latex')
xlabel('$y_0$ estimation error','Interpreter','Latex')
ylabel('Count','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
%%
save('basalests.mat','miny','levelestvect','b1estvect','b2estvect','truea');