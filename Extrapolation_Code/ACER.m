function [barrier_levels, eps_hat_mean, CI, Akj, Bkj] = ...
    ACER(process, ML, k_memory, N_barrier, conf_level, flagACER, flagCI)

% -------------------------------------------------------------------- %
%                                                                      %
% Calculating the ACER functions for the given array of numbers        %
%                                                                      %
%  Inputs:                                                             %
%   process: Extracted peaks from time series realizations of the      %
%            processes                                                 %
%   ML(ii): Length of peaks of realization "ii"                        %
%   k_memory: The number of previous non-exceedances-1                 %
%   N_barrier: Number of barrier levels                                %
%   conf_level: The confidence level e.g. 0.95 for 95%                 %
%               confidence intervals in papers                         %
%   flagACER: 1 or 2 - defines the way of calculating ACER                 %
%                                                                      %
%  Outputs:                                                            %
%   Barrier_levels: Equidistance barrier levels                        %
%   eps_hat_mean: Mean of the ACER functions; \epsilo_{k}^{hat}        %
%   CI: Confidence Intervals of the ACER functions                     %
%                                                                      %
%                                                                      %
Dim = size(process);  % [(Number of time series simulated)x(maximum number of extracted peaks)]
peaks = process;      % Peaks extracted from time series; the ingredients for estimation of ACER functions

%% Making equidistance barrier levels
barrier_levels = linspace(min(min(peaks))-eps,max(max(peaks))+eps,N_barrier);
D_eta = barrier_levels(2)-barrier_levels(1);   % Delta barrier level

%% Calculating ACER functions for each of the realizations
% tic;

Akj = zeros(Dim(1),N_barrier); % added 12.09.2011
Bkj = zeros(Dim(1),N_barrier); % added 14.09.2011

eps_hat = zeros(Dim(1),N_barrier);
eps_hat_ = zeros(Dim(1),N_barrier);
for jj=1:Dim(1)
    N_D_eta = (peaks(jj,:)-barrier_levels(1))/D_eta+1;   % Number of barrier levels each peak contains
    N_D_eta_floor = floor(N_D_eta);                      % Number of barrier levels each peak crosses
    N_D_eta_ceil = ceil(N_D_eta);                        % Number of barrier levels not met by each peak
    
    if k_memory==1
        num = zeros(size(barrier_levels));
        for ii=1:ML(jj);                                 % Attention: This counts on number of peaks
            counter = 1:N_D_eta_floor(ii);
            num(counter) = num(counter)+1;
        end
        den = ML(jj)*ones(size(num));                    % Attention: Should be devided by the number of peaks
    else
        num = zeros(size(barrier_levels));
        den = zeros(size(barrier_levels));
        kk = k_memory-1;
        for ii=kk+1:ML(jj)
            counter_den = max(N_D_eta_ceil(ii-kk:ii-1));
            den(counter_den:N_barrier) = den(counter_den:N_barrier)+1;
            counter_den = counter_den:N_D_eta_floor(ii);
            if ~isempty(counter_den)
                num(counter_den) = num(counter_den)+1;
            end
        end
        den(den==0) = eps;
    end
    
    Akj(jj,:) = num; % added 12.09.2011
    Bkj(jj,:) = den; % added 14.09.2011
    
    eps_hat(jj,:) = num./den;
    eps_hat_(jj,:) = num./((ML(jj)-k_memory+1)*ones(size(num)));
end
% toc;
clear ii jj kk;

%% Calculating Mean & Confidence Intervals of ACER functions
if (Dim(1)~=1)
    
    if (flagACER == 1)
        eps_hat_mean = mean(eps_hat);
        eps_hat_std = std(eps_hat);
        
%         logACER = log(eps_hat);
%         logACER(logACER == (-Inf)) = 0;%0; %NaN; %log(eps);
%         logACERvar = var(logACER);
%         for jj = 1:size(logACER,2)
%             dummy = logACER(:,jj);
%             dummy(isnan(dummy)) = [];
%             logACERvar(jj) = var(dummy);
%             coef(jj) = mean(dummy)./std(dummy);
% %             logACERvar(jj) = sum( ( dummy - sum(dummy)/Dim(1) ).^2 )/(Dim(1)-1);
%         end
        
    else % <=> flagACER == 2
        eps_hat_mean = mean(eps_hat_);
        eps_hat_std = std(eps_hat_);
        
%         logACER = log(eps_hat_);
%         logACER(logACER == (-Inf)) = 0;%0; %NaN; %log(eps);
%         logACERvar = var(logACER);
    end
    
    if conf_level>1 && conf_level<100
        conf_level = conf_level/100;
    elseif conf_level>100
        conf_level = 0.95;
    end
    %     conf_coef = norminv((1+conf_level)/2);       % Calculating confidence interval coefficient i.e. 1.96 for 95% confidence interval
    
    if (flagCI == 1)
        conf_coef = tinv((1+conf_level)/2,Dim(1)-1);
        CI = conf_coef*eps_hat_std/sqrt(Dim(1));     % Confidence interval
    else % <=> flagCI == 2
        conf_coef = norminv((1+conf_level)/2);       % Calculating confidence interval coefficient i.e. 1.96 for 95% confidence interval
        CI = conf_coef*sqrt(eps_hat_mean)/sqrt(sum(ML)-k_memory+1);
    end
    
else
    
    if (flagACER == 1)
        eps_hat_mean = eps_hat;
    else % <=> flagACER == 2
        eps_hat_mean = eps_hat_;
    end
    %         eps_hat_mean = counter0/(sum(ML)-k_memory+1);
  
    if conf_level>1 && conf_level<100
        conf_level = conf_level/100;
    elseif conf_level>100
        %         fprintf('Confidence level not in range [0,100].\n')
    end
    
    conf_coef = norminv((1+conf_level)/2);       % Calculating confidence interval coefficient i.e. 1.96 for 95% confidence interval
    CI = conf_coef*sqrt(eps_hat_mean)/sqrt(sum(ML)-k_memory+1);
end

end
