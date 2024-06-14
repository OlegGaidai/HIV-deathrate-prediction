% -------------------------------------------------------- %
%                                                          %
%   Main program which does the following tasks:           %
%                                                          %
%   1. Simulates time series of a Gaussian stochastic      %
%      process with a given power spectrum                 %
%   2. Calculates the ACER functions for it                %
%   3. Calculates the coefficients of the fit              %
%   4. Plots the data versus fit                           %
%                                                          %
%%
close all; clear all; clc
%% -------------------------------------------------------- %
% ---       Input data, constants and variables        --- %
% -------------------------------------------------------- %

x_star = @(level,x) x(2)+(1/x(3)*(log(x(1))-log(level)))^(1/x(4));
epsilon = @(eta,x) x(1)*exp(-x(3)*(eta-x(2)).^x(4));


[FileName, PathName] = uigetfile({'*.txt';'*.dat';'*.mat';'*.xls'},...
    'Open File','Multiselect','off');

if ~isempty(strfind(FileName, '.xls'))
    x = xlsread(strcat(PathName,FileName));
else
    x = importdata(strcat(PathName,FileName));
end

if  size(x,1)> size(x,2)
    x = transpose(x);
end
Dim = size(x);

choice = questdlg('Extract peaks from time series?','Choosing data','No');

% choice = 'No';

switch choice
    case 'Yes'
        % --- Extracting peaks from time series --- %
        process = zeros(Dim);
        ML = zeros(Dim(1),1);   % Maximum length of the peaks extracted from time series
        h = waitbar(0,'Extracting peaks, please wait');
        for jj=1:Dim(1)
            clear dummy
            dummy = peaks_from_timeseries(x(jj,:));
            ML(jj) = length(dummy);
            process(jj,1:ML(jj)) = dummy;
            waitbar(jj/Dim(1));
        end
        close(h);
        clear h;

        process = process(:,1:max(ML));
        helpdlg(horzcat(['N =',' ',num2str(sum(ML),'%10.4e\n')]),'Number of peaks');
    case {'No', 'Cancel'}
        process = x;
        ML = zeros(Dim(1),1);
        for jj=1:Dim(1)
            ML(jj) = sum(isfinite(process(jj,:))); 
        end
        helpdlg(horzcat(['N =',' ',num2str(sum(ML),'%10.4e\n')]),'Number of points');        
        
end

% --- input for the k_memory --- %
text = {'Define vector of k - subindexes of the desired ACER functioms \epsilon_k(\eta):'};
title = 'Input for the k values';
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
k_memory = inputdlg(text,title,1,{'1:6'},options);
k_memory = str2num(k_memory{1});                % Memory of the process to be taken into account
clear options text title

% k_memory = 1:6;

%% ------------------------------------------------------ %
% --- Estimation of theACER functions with "k" memory --- %
% ------------------------------------------------------- %
flagACERstr = questdlg('Is the time series stationary?','Stationarity','No');
if strcmp(flagACERstr, 'Yes')
    flagACER = 1;
else
    flagACER = 2;
end
h = waitbar(0,'Calculating ACER functions, please wait');
for kk=1:length(k_memory)
        [bar_lev{kk}, acer_hat{kk}, CI{kk}, Akj{kk}, Bkj{kk}] = ...
            ACER(process, ML, k_memory(kk), 750, 0.95, flagACER, 1);
        waitbar(kk/length(k_memory));
end
close(h);
clear kk h;

if Dim(1)~=1
    helpdlg(horzcat(['There are',' ',num2str(Dim(1)),' ','realizations. ',...
        'Empirical CI is estimated.']),'Note');
else
    helpdlg(horzcat(['There is only one realization provided. ',...
        'CI is estimated asymptotically under Poisson assumption (See User guide, p.7)']),'Note');
end

process = process(:)';
process(isnan(process)) = [];
 
x = x(:)';
x(isnan(x)) = [];
    
sigma = std(x);
%% ----------------------------------------------------- %
% --- Plotting of the ACER functions with "k" memory --- %
% ------------------------------------------------------ %
MrKr = {'o','>','-.',':','s','-','+','d','p','x','v','^','h','s'};
if length(k_memory)>length(MrKr)
    [MK{1:length(k_memory)-length(MrKr)}] = deal('-');
    MrKr = [MrKr MK];
end
clear MK

MrkSz = [3 3 3 3 2.5 3 3 3 3 3 3 3 3 3];
if length(k_memory)>length(MrkSz)
    [MS(1:length(k_memory)-length(MrkSz))] = deal(3);
    MrkSz = [MrkSz MS];
end
clear MS

line = [0.5 0.5 1 1.5 0.5 1 1 1 1 1 1 1 1 1];
if length(k_memory)>length(line)
    [ln(1:length(k_memory)-length(line))] = deal(1);
    line = [line ln];
end
clear ln
% 
% figure
% clf
% for ii=1:length(k_memory)
%     
%     %     condition = bar_lev{ii}>bar_lev{ii}(15);
%     %     CI{ii} = CI{ii}(condition);
%     
%     semilogy(bar_lev{ii}(bar_lev{ii}>bar_lev{ii}(20)),...
%         acer_hat{ii}(bar_lev{ii}>bar_lev{ii}(20)),...
%         MrKr{ii},'Color','k','MarkerSize',MrkSz(ii),'Linewidth',line(ii))
%     hold on
%     leg{ii} = ['k=' num2str(k_memory(ii))];
%     
%     eps_menu_options{ii} = ['ACER(k=' num2str(k_memory(ii)) ')'];
%     
% end
% axis tight
% xlabel('\eta')
% ylabel('ACER_k(\eta)')
% legend(leg)
% editplot;
% ax = axes('Units', 'Normalized', 'Position', [0.13 0.13 0.25 0.25], ...
%     'XTick', [], 'YTick', [], 'Box', 'on');
% hold on;
% for ii=1:length(k_memory)
%     
%     %     condition = bar_lev{ii}>bar_lev{ii}(15);
%     %     CI{ii} = CI{ii}(condition);
%     
%     semilogy(bar_lev{ii}(bar_lev{ii}>bar_lev{ii}(20)),...
%         log10(acer_hat{ii}(bar_lev{ii}>bar_lev{ii}(20))),...
%         MrKr{ii},'Color','k','MarkerSize',MrkSz(ii),'Linewidth',line(ii))
%     hold on
%     
% end


figure(1)
clf
for ii=1:length(k_memory)
    
    %     condition = bar_lev{ii}>bar_lev{ii}(15);
    %     CI{ii} = CI{ii}(condition);
    
    semilogy(bar_lev{ii}(bar_lev{ii}>bar_lev{ii}(20))/sigma,...
        acer_hat{ii}(bar_lev{ii}>bar_lev{ii}(20)),...
        MrKr{ii},'Color','k','MarkerSize',MrkSz(ii),'Linewidth',line(ii))
    hold on
    leg{ii} = ['k=' num2str(k_memory(ii))];
    
    eps_menu_options{ii} = ['ACER(k=' num2str(k_memory(ii)) ')'];
    
end
axis tight
xlabel('\eta/\sigma')
ylabel('ACER_k(\eta)')
legend(leg)
editplot;
ax = axes('Units', 'Normalized', 'Position', [0.13 0.13 0.25 0.25], ...
    'XTick', [], 'YTick', [], 'Box', 'on');
hold on;
for ii=1:length(k_memory)
    
    %     condition = bar_lev{ii}>bar_lev{ii}(15);
    %     CI{ii} = CI{ii}(condition);
    
    semilogy(bar_lev{ii}(bar_lev{ii}>bar_lev{ii}(20))/sigma,...
        log10(acer_hat{ii}(bar_lev{ii}>bar_lev{ii}(20))),...
        MrKr{ii},'Color','k','MarkerSize',MrkSz(ii),'Linewidth',line(ii))
    hold on
    
end


clear Dim rnd dummy kk ii jj leg MrKr
%% ------------------------------------------------- %
% ---    Choosing the desired ACER function     --- %
% ------------------------------------------------- %

choice_of_ACER = menu('Choose ACER function to be analysed',eps_menu_options);

bar_lev = bar_lev{choice_of_ACER};
acer_hat = acer_hat{choice_of_ACER};
CI = CI{choice_of_ACER};

Akj = Akj{choice_of_ACER};
Bkj = Bkj{choice_of_ACER};
eta = bar_lev;
[~, position] = max(Akj(:,find(eta>=mean(x),1,'first'):end),[],2);
ind = ceil(mean(position)+5*std(position));

clear eps_menu_options Akj Bkj position;
%% ------------------------------------------------- %
% ---    Cutting data from the tail marker      --- %
% ---                \eta_1                     --- %
% ------------------------------------------------- %

figure(2)
clf
semilogy(bar_lev(bar_lev>bar_lev(15)),acer_hat(bar_lev>bar_lev(15)),'o','Color','k','MarkerSize',3)
axis tight
xlabel('\eta')
ylabel(['ACER_{' num2str(k_memory(choice_of_ACER)) '}(\eta)'])
legend(['ACER_k(\eta) for k=' num2str(k_memory(choice_of_ACER))])
editplot;

text = {'Enter value of the tail marker \eta_1:'};
title = 'Input for the tail marker';
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
if k_memory(choice_of_ACER)~=1
    eta_1 = eta(eta>=mean(x));
    if numel(eta_1)>=ind
        eta_1 = inputdlg(text,title,1,{num2str(eta_1(ind))},options);
    else
        eta_1 = inputdlg(text,title,1,{''},options);
    end
else
    eta_1 = inputdlg(text,title,1,{''},options);
end

eta_1 = str2num(eta_1{1});

condition = bar_lev>eta_1;               % Finding the index of the tail marker

bar_lev = bar_lev(condition);     % Refining barrier levels
CI = CI(condition);
acer_hat = acer_hat(condition);         % Refining ACER function

clear condition options text title
%% ------------------------------------------------- %
% ---    Cutting from consideration data        --- %
% ---         with high uncertainty             --- %
% ------------------------------------------------- %

text = {'Define delta :'};
title = 'Input for the \delta value';
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
delta = inputdlg(text,title,1,{'1'},options);
delta = str2num(delta{1});                % Memory of the process to be taken into account
clear options text title


% delta = 1;                                   % parametr from formula 27 "Estimation of extreme values 
                                             % from sampled time series. Structural Safety 31, 325-334."
condition = ( CI < (acer_hat)*delta );
condition = updatecond(condition);
bar_lev = bar_lev(condition);
acer_hat = acer_hat(condition);
CI = CI(condition);

CI_plus = acer_hat + CI;
CI_minus = acer_hat - CI;

clear condition text title

%%
text = {'Enter value of the desired level you want ACER function to be extrapolated to:'};
title = 'Input for the value of the target level';
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
level_of_interest = inputdlg(text,title,1,{'1e-6'},options);
level_of_interest = str2num(level_of_interest{1});                % Memory of the process to be taken into account

% level_of_interest = 1e-6;

clear options text title
%% ------------------------------------------------- %
% ---              optimization                  --- %
% -------------------------------------------------- %
penalty = questdlg('Use penalized objective function:','Objective function','No');
pause(2);
h = waitbar(0,'Optimization routine has started, please wait','Name','Optimal curve fitting...');
pause(1);

W = (log(CI_plus)-log(CI_minus)).^(-2);
% Checking for all weights to be "positive" and "real"
pos_index = W>0;
real_index = imag(W)~=0; 
if ~isempty(pos_index)
    %fprintf('Non-positive weights of the objective function changed to positive.\n')
    W = abs(W); 
elseif ~isempty(real_index)
    %fprintf('Imaginary weights of the objective function changed to real.\n')
    W = abs(W); 
end
% Normalizing weights to have sum(w)=1
W = W/sum(W);
clear pos_index real_index;
%------------------------------------------------------
bmin = min(min(x));
%------------------------------------------------------
if exist('x','var') ~= 0
    [q0, b0, a0, c0, R] = guess(bar_lev, eta_1, x, acer_hat);
    [~, ind] = min(abs(bar_lev - mean2(x)));
else
    [q0, b0, a0, c0, R] = guess(bar_lev, eta_1, process, acer_hat);
    [~, ind] = min(abs(bar_lev - mean2(process)));
end
qMedian = acer_hat(ind); % ACER(mean of process)

%------------------------------------------------------
% qMin = exp(max(log(acer_hat))); % minimal q value
% q_init = max(qMin, q0);
% qMax =  qMedian - (qMin - qMedian);
waitbar(0.15,h,'Step 1 of 3...');
[fin_sol, fmin, pos, sol] = Optimization(bar_lev, ...
    acer_hat, eta_1, W, q0, b0, a0, c0, qMedian, bmin, 0.05, penalty);
acer_hat_fit = epsilon(bar_lev,fin_sol);
x_star_main = x_star(level_of_interest,fin_sol);

for ss=1:size(sol,1)
    x_eps(ss) = x_star(level_of_interest,sol(ss,1:4));
end

clear ss q0 b0 a0 c0 qMedian R;
%% ------------------------------------------------- %
% ---                Reanchoring                 --- %
% -------------------------------------------------- %
CI_plus_reanch = acer_hat_fit + CI;
CI_minus_reanch = acer_hat_fit - CI;

%% -------------------------------------------------- %
% ---      optimization for CI_plus_reanch        --- %
% --------------------------------------------------- % 
if exist('x','var') ~= 0
    [q0_pl, b0_pl, a0_pl, c0_pl, R_pl] = guess(bar_lev, eta_1, x, CI_plus_reanch);
    [~, ind] = min(abs(bar_lev - mean2(x)));
else
    [q0_pl, b0_pl, a0_pl, c0_pl, R_pl] = guess(bar_lev, eta_1, process, CI_plus_reanch);
    [~, ind] = min(abs(bar_lev - mean2(process)));
end
qMedian_pl = CI_plus_reanch(ind); % ACER(mean of process)

waitbar(1/3,h,'Step 2 of 3...');
[fin_sol_pl, fmin_pl, pos_pl, sol_pl] = Optimization(bar_lev, ...
    CI_plus_reanch, eta_1, W, q0_pl, b0_pl, a0_pl, c0_pl, qMedian_pl, bmin, 0.05, penalty);
x_CI_pl = x_star(level_of_interest,fin_sol_pl);
clear q0_pl b0_pl a0_pl c0_pl qMedian_pl R_pl;

%%
condition = (CI_minus_reanch > 0); %& (CI_plus_reanch > 0);
condition = updatecond(condition);

CI_minus_reanch = CI_minus_reanch(condition);
CI_plus_reanch = CI_plus_reanch(condition);

bar_lev = bar_lev(condition);
W = W(condition);
acer_hat = acer_hat(condition);
acer_hat_fit = acer_hat_fit(condition);
CI_plus = CI_plus(condition);
CI_minus = CI_minus(condition);
clear condition;
%% -------------------------------------------------- %
% ---      optimization for CI_minus_reanch       --- %
% --------------------------------------------------- % 
if exist('x','var') ~= 0
    [q0_mn, b0_mn, a0_mn, c0_mn, R_mn] = guess(bar_lev, eta_1, x, CI_minus_reanch);
    [~, ind] = min(abs(bar_lev - mean2(x)));
else
    [q0_mn, b0_mn, a0_mn, c0_mn, R_mn] = guess(bar_lev, eta_1, process, CI_minus_reanch);
    [~, ind] = min(abs(bar_lev - mean2(process)));
end
qMedian_mn = CI_minus_reanch(ind); % ACER(mean of process)

waitbar(2/3,h,'Step 3 of 3...');
[fin_sol_mn, fmin_mn, pos_mn, sol_mn] = Optimization(bar_lev, ...
    CI_minus_reanch, eta_1, W, q0_mn, b0_mn, a0_mn, c0_mn, qMedian_mn, bmin, 0.05, penalty);
x_CI_mn = x_star(level_of_interest,fin_sol_mn);
clear q0_mn b0_mn a0_mn c0_mn qMedian_mn R_mn;

waitbar(3/3,h,'Finished!');
pause(2);
close(h);
clear h;
%% ------------------------------------------------- %
% --- Plotting the fit with extrapolated        --- %
% ---         confidence intervals              --- %
% ------------------------------------------------- %
if strcmp(penalty, 'Yes')
    fprintf('General solution, penalized: \r\n');
else
    fprintf('General solution: \r\n');
end

finalfigure(level_of_interest, fin_sol, fin_sol_pl, fin_sol_mn,...
    bar_lev, acer_hat, CI_plus_reanch, CI_minus_reanch, k_memory, choice_of_ACER, sigma)

% fprintf('General solution, penalized: \r\n');
% finalfigure(level_of_interest, fin_sol_penal, fin_sol_pl_penal, fin_sol_mn_penal,...
%     bar_lev, acer_hat, CI_plus_reanch, CI_minus_reanch, sigma)

fid = fopen(strcat(PathName,FileName(1:end-4),...
    '_ACER_k',num2str(k_memory(choice_of_ACER)),'_results.txt'), 'w');
fprintf(fid,'----------------------------------------------------------------------------\r\n');
fprintf(fid,'                                WORK STATEMENT                              \r\n');
fprintf(fid,'----------------------------------------------------------------------------\r\n');
fprintf(fid,' \r\n');
fprintf(fid,'Input data:\r\n');
fprintf(fid,' \r\n');

fprintf(fid,'Time series loaded from:                              %s\r\n', strcat(PathName,FileName));
fprintf(fid,'Extraction of peaks:                                  %s\r\n', choice);
fprintf(fid,'Stationarity of the loaded time series:               %s\r\n', flagACERstr);
fprintf(fid,' \r\n');

fprintf(fid,'Vector of k:                                          [%s]\r\n', num2str(k_memory));
fprintf(fid,'Analysis was made for:                                %s\r\n', strcat('ACER(k=',num2str(k_memory(choice_of_ACER)),')'));
fprintf(fid,'Tail marker:                                          %4.3f\r\n', eta_1);
fprintf(fid,'Level of cutting uncertain data:                      %4.3f\r\n', delta);
fprintf(fid,'Level of interest:                                    %2.3e\r\n',level_of_interest);
fprintf(fid,'Power of weights (1 or 2):                            %i\r\n', 2);
fprintf(fid,'Use the penalized objective function:                 %s\r\n', penalty);
fprintf(fid,' \r\n');
fprintf(fid,'----------------------------------------------------------------------------\r\n');
fprintf(fid,' \r\n');
fprintf(fid,'Output results:\r\n');
fprintf(fid,' \r\n');
fprintf(fid,'Max. value of the loaded process:                     %g\r\n', max(max(x)));
fprintf(fid,'Min. value of the process:                            %g\r\n', min(min(x)));
fprintf(fid,'Mean value of the process:                            %g\r\n', mean2(x));
fprintf(fid,'Standard deviation:                                   %g\r\n', sigma);
fprintf(fid,' \r\n');
fprintf(fid,'Predicted T-years return level estimate is:           %g\r\n', x_star_main);
fprintf(fid,'Predicted confidence interval:                  CI_ = %g\r\n', x_CI_mn);
fprintf(fid,'                                                CI+ = %g\r\n', x_CI_pl);
fprintf(fid,' \r\n');
fprintf(fid,'Parameters of optimal curve are:                  q = %g\r\n', fin_sol(1));
fprintf(fid,'                                                  b = %g\r\n', fin_sol(2));
fprintf(fid,'                                                  a = %g\r\n', fin_sol(3));
fprintf(fid,'                                                  c = %g\r\n', fin_sol(4));
fprintf(fid,' \r\n');
fprintf(fid,' \r\n');
fprintf(fid,datestr(now,21));
fclose(fid);
