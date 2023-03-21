% Parameter and input estimation for first order system
% Monte Carlo runs with synthetic data with added noise

rng(12)

n_runs=400; % number of runs

n_impulses=3;
dmin=0.1; % lowest impulse amplitude
dmax=1;

Tmax = 5; 
Tmin = 2; % Distance between impulses uniformly distributed between
                % min_time and timescale
min_a1=0.4; % lowest a1 value. a1 uniformly distributed between min_a1 and
            % min_a1+1
delta_t_sampled=0.5; % sampling time

noisestd=0.01; % noise std

deltabase=zeros(1,n_runs); % error when using min f/(df/db)
deltaadj=zeros(1,n_runs); % error when correcting with newton step
delta2=zeros(1,n_runs);

ressumest=zeros(1,n_runs);
nstepmin_= zeros(1,n_runs);
ressum_min = zeros(1,n_runs);
ressum_minDelta = zeros(1,n_runs);
delta_est = zeros(1,n_runs);

for k=1:n_runs
    %% Base parameters
    a1= rand(1)+min_a1;
    a2= 2*a1; %not used
        
    [y_nonoise,t_sampled,x,t,dSeq] = generatedata(a1,a2,n_impulses,dmin,dmax,Tmin,Tmax,0.02,delta_t_sampled,true);
    gaussiannoise = noisestd*randn(size(t_sampled));
    y = y_nonoise+gaussiannoise';
    tk = t_sampled(2:end-1);
    m = length(y);
    n = length(tk);
    options = optimset('lsqlin');
    options = optimset(options,'Display','off','TolFun',1e-10,'TolCon',1e-10);%'Algorithm','interior-point',
    
    %%
    a1_min=0.01;
    a1_max=(a1+a2)/2;
    a_delta=0.001;
    
    a1_range=a1_min:a_delta:a1_max;

    ressum=zeros(numel(a1_range),1);
    ressum_=zeros(numel(a1_range),1);

    beta=zeros(numel(a1_range),n+1);
    beta0=[];
    beta0_=[];
    
    parfor k1=1:numel(a1_range)
        a1_=a1_range(k1);
            Z = zeros(m,n+1); % Unkown impulse times
            Z(:, 1) = exp(-a1_*t_sampled);
            for i=2:n+1
                tarr = t_sampled - tk(i-1);
                firstp = find(tarr>=0,1,'first');
                Z(:, i) = [zeros(firstp-1,1); exp(-a1_*tarr(firstp:end))'];
            end
            Z_ = zeros(m,2); % known impulse times
            for i=1:n_impulses
                tarr = t_sampled - dSeq(1,i);
                firstp = find(tarr>=0,1,'first');
                Z_(:, i) = [zeros(firstp-1,1); exp(-a1_*tarr(firstp:end))'];
            end
            
            [beta01,ressum(k1)]=lsqlin(Z,y,[],[],[],[],zeros(n+1,1),...
                ones(n+1,1)*10000,beta0,options); % Unkown impulse times
            [beta01_,ressum_(k1)]=lsqlin(Z_,y,[],[],[],[],zeros(n_impulses,1),...
                ones(n_impulses,1)*10000,beta0_,options); % known impulse times

            beta(k1,:)=beta01;
    end

    ressumDer = (ressum(3:end)-ressum(1:end-2))/(2*a_delta);
    A = find(ressumDer>0,1);
    ressum(A:end) = NaN;
    %ressum(sum(sparsity,2)>=0.5*numel(t_sampled))=NaN;

    ressum_min(k) = min(ressum_); % impulse time is known
    ressum_minDelta(k) =a1 - a1_range(ressum_== ressum_min(k)); % estimation error when impulse time is known
    
    ntstep=-(ressum(2:end-1)+noisestd^2)./(ressum(3:end)-ressum(1:end-2))*2*a_delta;
    ntstep_=-ressum_(2:end-1)./(ressum_(3:end)-ressum_(1:end-2))*2*a_delta;
    ntstep_(ntstep_<0)=NaN;

    minidx=find(ntstep==min(ntstep));
    minidx_=find(ntstep_==min(ntstep_));
    nstepmin=ntstep(minidx);
    nstepmin_(k)=ntstep_(minidx_);
    minidx=minidx+1;
    
    deltabase(k)=a1-a1_range(minidx);
    deltaadj(k)=a1-a1_range(minidx)-nstepmin; % error with unknown impulse times
    delta_est(k)=a1-(a1_range(minidx_)+nstepmin_(k)); % error with one step and known impulse time. for some reason this estimate is a lot worse
    
end
%%
figure(5)
histogram(deltabase,25)%,'Normalization','probability'
hold on
histogram(deltaadj,15)%,'Normalization','probability'
hold off

legend('$\bar \omega$','$\hat \omega$','Interpreter','Latex')%,'Fitted normal distribution')
xlabel('Estimation error','Interpreter','Latex')
ylabel('Count','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

figure(6)
plot(deltaadj, ressum_minDelta,'k.')
xlabel('One-step estimation error','Interpreter','Latex')
ylabel('Estimation error with $f^\dagger(\omega)$ ','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')