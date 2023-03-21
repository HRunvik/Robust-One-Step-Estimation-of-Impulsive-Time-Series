function [theta_new, p_new, Theta, P, ressum] = robust_constrained_linReg(x_aug,y_aug,theta_old,epsilon,lsopts)
%INPUT
%x_aug : n x d matrix
%y_aug : n x 1 vector
%epsilon : scalar

%Output:
%theta_new : d x 1 vector of robust estimate

n = size(x_aug,1);
n2=size(x_aug,2);

cnt = 1;
P = [1/n.*ones(n,1)]; Theta = [theta_old];

% optimize
while (cnt==1)
    alpha = (y_aug-x_aug*theta_old).^2;
    
    %p_new = optimize_probability(alpha,epsilon,n);
    
    p_new = my_barrier_opt_prb(alpha,epsilon,n);
    
    %p_new = alternate_method_prb_opt(alpha, epsilon);
    
    P = [P p_new];
        
    [theta_new,ressum] = lsqlin(sqrt(p_new).*x_aug,sqrt(p_new).*y_aug,...
        [],[],[],[],zeros(n2,1),ones(n2,1)*10000,[],lsopts);

    Theta = [Theta theta_new];
    
    if norm(theta_old-theta_new)/norm(theta_old)>1e-3&&size(Theta,2)<100
        cnt = 1;
        theta_old = theta_new;
    else
        cnt = 0;
    end
end
