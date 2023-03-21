% Parameter and input estimation for second order system
% Monte Carlo runs with synthetic data with added noise
% Added outliers, fixed basal level

n_runs=50; % number of runs

n_impulses=3; % number of impulses
dmin=0.1*4; % lowest impulse amplitude
dmax = 1*4; % Impulse amplitudes uniformly distributed between dmin
% and dmax

Tmax = 5;
Tmin = 2; % Distance between impulses uniformly distributed between
% Tmin and Tmax

min_a1=0.4; % lowest a1 value. a1 uniformly distributed between min_a1 and
% min_a1+1
min_adiff = 0.3; % lowest difference between a1 and a2. a2 uniformly
% distributed between a1+min_adiff and a1+min_adiff+1

delta_t_sampled=0.5; % sampling time

options = optimset('lsqlin');
options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);

truea=zeros(n_runs,2);

minerrorgamma_iid=zeros(n_runs,1);
minerrorgamma_mixed=zeros(n_runs,1);
minerrorgamma_robust=zeros(n_runs,1);

a1estvect_iid = zeros(n_runs,1);
a2estvect_iid = zeros(n_runs,1);
a1estvect_robust = zeros(n_runs,1);
a2estvect_robust = zeros(n_runs,1);

mnoisestd= 0.006; % noise std
n_outliers = 2;
outlierrange = dmax*0.125;

for k=1:n_runs
    disp(['Run ' num2str(k)])
    rng((k-1)*5)
    a1= rand(1)+min_a1;
    a2= a1+min_adiff+rand(1);
    truea(k,:)=[a1 a2];

    [y_nonoise,t_sampled,x,t] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,0.02,delta_t_sampled);

    gaussiannoise = mnoisestd*randn(size(t_sampled));
    y_iidnoise = y_nonoise+gaussiannoise';

    outlierix = randperm(numel(y_nonoise));
    outlierix = outlierix(1:n_outliers);
    outliernoise = gaussiannoise;
    outliernoise(outlierix) = (rand(1,2)-0.5)*outlierrange;
    y_outliers = y_nonoise+outliernoise';

    %% Prepare stuff
    tk = t_sampled(2:end-1);
    m = length(y_iidnoise);
    n = length(tk);

    a1_min=a1*0.5;
    a1_max=(a1+a2)/2;
    a_delta=0.002;%0.01
    a2_min=a1_max+a_delta;
    a2_max=a2*1.5;
    a1_range=a1_min:a_delta:a1_max;
    a2_range=a2_min:a_delta:a2_max;

    %% a1-a2
    tic
    [a1est_iid,m1val_iid] = gammaCurve(a1_range,a2_range,y_iidnoise,t_sampled,options,15,mnoisestd^2);
    a1est_iid=a1est_iid(2:end-1);
    m1val_iid=m1val_iid(2:end-1);

    a1estadj_iid=a1est_iid+m1val_iid;

    [a1est_outliers,m1val_outliers] = gammaCurve(a1_range,a2_range,y_outliers,t_sampled,options,15,mnoisestd^2);
    a1est_outliers=a1est_outliers(2:end-1);
    m1val_outliers=m1val_outliers(2:end-1);

    a1estadj_outliers=a1est_outliers+m1val_outliers;

    [a1est_robust,m1val_robust] = gammaCurve_robust(a1_range,a2_range,y_outliers,t_sampled,options,15,0.1,mnoisestd^2/numel(y_outliers));
    a1est_robust=a1est_robust(2:end-1);
    m1val_robust=m1val_robust(2:end-1);
    a1estadj_robust = a1est_robust+m1val_robust;

    toc

    errordist = ((a1estadj_iid-a1).^2+(a2_range(2:end-1)-a2).^2).^0.5;
    minerrorgamma_iid(k) = min(errordist);
    minerrorgamma_mixed(k) = min(((a1estadj_outliers-a1).^2+(a2_range(2:end-1)-a2).^2).^0.5);
    minerrorgamma_robust(k) = min(((a1estadj_robust-a1).^2+(a2_range(2:end-1)-a2).^2).^0.5);

    a_delta_ = 0.001;
    [ressumb_iid,betasumb_iid,tkcell_iid] = regularizeinput(a1est_iid,m1val_iid,a2_range,y_iidnoise,t_sampled,a_delta_,options,mnoisestd^2);
    [ressumb_robust,betasumb_robust,tkcell_robust] = regularizeinput_robust(a1est_robust,m1val_robust,a2_range,y_outliers,t_sampled,a_delta_,options,mnoisestd^2/numel(y_outliers),0.1);

    bic_iid=[];
    ressumbic_iid=[];
    minix_iid=zeros(1,max(betasumb_iid)-min(betasumb_iid));
    for k1 = min(betasumb_iid):max(betasumb_iid)
        ressumtemp = ressumb_iid;
        ressumtemp(betasumb_iid>k1) = inf;
        [~,minix_] = min(ressumtemp);
        if min(ressumtemp)<inf && betasumb_iid(minix_)==k1
            bic_iid = [bic_iid; [2*k1*log(m) + m*log(min(ressumtemp)) k1]];
            minix_iid(k1-min(betasumb_iid)+1)=minix_;
            ressumbic_iid = [ressumbic_iid min(ressumtemp)];
        end
    end
    minix_iid(minix_iid==0)=[];
    bic_iid(bic_iid(:,2)>m/4,1) = NaN;
    [~,minbicix_iid] = min(bic_iid(:,1));

    a2estvect_iid(k) = a2_range(minix_iid(minbicix_iid)+1);
    a1estvect_iid(k) = a1estadj_iid(minix_iid(minbicix_iid));

    bic_robust=[];
    ressumbic_robust=[];
    minix_robust=zeros(1,max(betasumb_robust)-min(betasumb_robust));
    for k1 = min(betasumb_robust):max(betasumb_robust)
        ressumtemp = ressumb_robust;
        ressumtemp(betasumb_robust>k1) = inf;
        [~,minix_] = min(ressumtemp);
        if min(ressumtemp)<inf && betasumb_robust(minix_)==k1
            bic_robust = [bic_robust; [2*k1*log(m) + m*log(min(ressumtemp)) k1]];
            minix_robust(k1-min(betasumb_robust)+1)=minix_;
            ressumbic_robust = [ressumbic_robust min(ressumtemp)];
        end
    end
    minix_robust(minix_robust==0)=[];
    bic_robust(bic_robust(:,2)>m/4,1) = NaN;
    [~,minbicix_robust] = min(bic_robust(:,1));

    a2estvect_robust(k) = a2_range(minix_robust(minbicix_robust)+1);
    a1estvect_robust(k) = a1estadj_robust(minix_robust(minbicix_robust));

end
%%
errora1_iid = truea(:,1)-a1estvect_iid;
errora2_iid = truea(:,2)-a2estvect_iid;
erroreuclid_iid = (errora1_iid.^2+errora2_iid.^2).^0.5;
errora1_robust = truea(:,1)-a1estvect_robust;
errora2_robust = truea(:,2)-a2estvect_robust;
erroreuclid_robust = (errora1_robust.^2+errora2_robust.^2).^0.5;