% Parameter and input estimation for second order system
% Monte Carlo runs with synthetic data with added noise
% No outliers, fixed basal level

n_runs=100; % number of runs

n_impulses=3; % number of impulses
dmin=0.1*4; % lowest impulse amplitude

Tmax = 5;
Tmin = 2; % Distance between impulses uniformly distributed between
% Tmin and Tmax

min_a1=0.4; % lowest a1 value. a1 uniformly distributed between min_a1 and
% min_a1+1
min_adiff = 0.3; % lowest difference between a1 and a2. a2 uniformly
% distributed between a1+min_adiff and a1+min_adiff+1
dmax = 1*4; % Impulse amplitudes uniformly distributed between dmin
% and dmax
endtime = 3; % fixed time added after last impulse
delta_t_sampled=0.5; % sampling time

k0max=6;
meanminsteperror=zeros(1,k0max);
meanBICminerror=zeros(1,k0max);
meanminpointerror=zeros(1,k0max);
meanimps=zeros(1,k0max);
meanmingammadist=zeros(1,k0max);

options = optimset('lsqlin');
options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);

for k0=1:k0max
    truea=zeros(n_runs,2);
    a1estcell = cell(n_runs,1);
    a2estcell = cell(n_runs,1);
    a1estvect = zeros(n_runs,1);
    a2estvect = zeros(n_runs,1);
    n_impest=zeros(n_runs,1);
    n_impmin=zeros(n_runs,1);
    a1estmin=zeros(n_runs,1);
    a2estmin=zeros(n_runs,1);
    minerrorgamma=zeros(n_runs,1);
    meanerrorgamma=zeros(n_runs,1);
    minsteperror=zeros(n_runs,1);

    mnoisestd= k0*0.002+0; % noise std

    for k=1:n_runs
        disp(['Run ' num2str(k)])
        rng((k-1)*5)
        a1= rand(1)+min_a1;
        a2= a1+min_adiff+rand(1);
        truea(k,:)=[a1 a2];
        
        [y_nonoise,t_sampled,x,t] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,0.02,delta_t_sampled);
        gaussiannoise = mnoisestd*randn(size(t_sampled));
        y = y_nonoise+gaussiannoise';

        tk = t_sampled(2:end-1);
        m = length(y);
        n = length(tk);
        
        a1_min=a1*0.5;
        a1_max=(a1+a2)/2;
        a_delta=0.002;%0.01
        a2_min=a1_max+a_delta;
        a2_max=a2*1.5;
        a1_range=a1_min:a_delta:a1_max;
        a2_range=a2_min:a_delta:a2_max;

        ressum=zeros(numel(a1_range),numel(a2_range));

        %% a1-a2
        tic
        [a1est,m1val] = gammaCurve(a1_range,a2_range,y,t_sampled,options,15,mnoisestd^2);
        a1est=a1est(2:end-1);
        m1val=m1val(2:end-1);
        toc
        a1estadj=a1est+m1val;

        a_delta_ = 0.001;
        [ressumb,betasumb,tkcell] = regularizeinput(a1est,m1val,a2_range,y,t_sampled,a_delta_,options,mnoisestd^2);

        bic=[];
        minix=zeros(1,max(betasumb)-min(betasumb));
        for k1 = min(betasumb):max(betasumb)
            ressumtemp = ressumb;
            ressumtemp(betasumb>k1) = inf;
            [~,minix_] = min(ressumtemp);
            if min(ressumtemp)<inf && betasumb(minix_)==k1
                bic = [bic; [2*k1*log(m) + m*log(min(ressumtemp)) k1]];
                minix(k1-min(betasumb)+1)=minix_;
            end
            if k1==min(betasumb)
                n_impmin(k)=k1;
                a1estmin(k)=a1estadj(minix_);
                a2estmin(k)=a2_range(minix_+1);
            end
        end
        minix(minix==0)=[];
        bic(bic(:,2)>m/4,1) = NaN;
        [~,minbicix] = min(bic(:,1));
        a2estcell{k} = a2_range(minix+1);
        a1estcell{k} = a1estadj(minix);
        a2estvect(k) = a2_range(minix(minbicix)+1);
        a1estvect(k) = a1estadj(minix(minbicix));
        n_impest(k) = bic(minbicix,2);

        errordist = ((a1estadj-a1).^2+(a2_range(2:end-1)-a2).^2).^0.5;
        minerrorgamma(k) = min(errordist);
        meanerrorgamma(k) = mean(errordist);
        [~,minstepix]=min(m1val);
        minsteperror(k)=((a1estadj(minstepix)-a1)^2+(a2_range(minstepix+1)-a2)^2)^0.5;

    end
    %%
%     figure(11)
%     histogram(n_impest)
%     hold on
%     histogram(n_impmin)
%     hold off

    errors = ((a1estvect-truea(:,1)).^2+(a2estvect-truea(:,2)).^2).^0.5;
    errormin = ((a1estmin-truea(:,1)).^2+(a2estvect-truea(:,2)).^2).^0.5;
    errorclosest = zeros(1,n_runs);
    for k=1:n_runs
        errorclosest(k) = (min((a1estcell{k}-truea(k,1)).^2 + (a2estcell{k}-truea(k,2)).^2))^0.5;
    end
    % figure(12)
    % histogram(errors)
    % hold on
    % histogram(errorclosest)
    % histogram(minsteperror)
    % hold off

    meanminsteperror(k0)=mean(minsteperror);
    meanBICminerror(k0)=mean(errors);
    meanminpointerror(k0)=mean(errorclosest);
    meanimps(k0)=mean(n_impest);
    meanmingammadist(k0)=mean(minerrorgamma);


end

%%
save('pointests.mat',"meanmingammadist","meanimps","meanminsteperror","meanBICminerror","meanminpointerror")