function [fin_sol, fmin, pos, sol] = Optimization(eta, ACER, eta_1, W, q0, b0, a0, c0, qMed, bmin, alpha, penalty)

sol = 100*ones(6,5);
warning off all;

%% 1) Global 4 param LOG-level constraint optimization
if strcmp(penalty, 'Yes')
    F = @(x) exp(alpha.*(x(4).^sign(x(4) - 1) - 1)) * sum(W.*((log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4))).^2); % x = [q b a c]
else
    F = @(x) sum(W.*((log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4))).^2); % x = [q b a c]
end

nlnineq = @(x) [log(x(1)) - x(3).*(eta - x(2)).^x(4)];
nlneq = @(x) [];
nonlinfcn = @(x)deal(nlnineq(x),nlneq(x));

optst_logGlobal = optimset('Display','off',...
    'Algorithm', 'interior-point', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol(1,1:4), sol(1,5)] = fmincon(F,[q0 b0 a0 c0],[],[],[],[],...
    [eps bmin eps 1],...
    [Inf eta_1 Inf Inf],...
    nonlinfcn,optst_logGlobal);
if  any(imag(sum(sol(1,:))))
    sol(1,:) = 100*ones(1,5);
end
clear F optst_logGlobal nlnineq nlneq nonlinfcn;

%% 2) Global 4 param LOG-level LS constraint optimization = 5

if strcmp(penalty, 'Yes')
    F = @(x) sqrt(exp(alpha.*(x(4).^sign(x(4) - 1) - 1))) * sqrt(W).*(log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4)); % x = [q b a c]
else
    F = @(x) sqrt(W).*(log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4)); % x = [q b a c]
end

optst_log_all = optimset('Display','off',...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol(2,1:4), sol(2,5)] = lsqnonlin(F, [q0 b0 a0 c0], ...
    [eps bmin eps 1],...
    [Inf eta_1 Inf Inf],...
    optst_log_all);
if  any(imag(sum(sol(2,:))))
    sol(2,:) = 100*ones(1,5);
end
clear F optst_log_all;

%% 3) Global 4 param LOG-level LS-optimization Levenberg-Marquard
if strcmp(penalty, 'Yes')
    F = @(x) sqrt(exp(alpha.*(x(4).^sign(x(4) - 1) - 1))) * sqrt(W).*(log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4)); % x = [q b a c]
else
    F = @(x) sqrt(W).*(log(ACER)-log(x(1))+...
        x(3)*(eta - x(2)).^x(4)); % x = [q b a c]
end

optst_logLM_all = optimset('Display','off',...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol(3,1:4), sol(3,5)] = lsqnonlin(F, [q0 b0 a0 c0], [], [], optst_logLM_all);
if  any(imag(sum(sol(3,:))))
    sol(3,:) = 100*ones(1,5);
end

clear F optst_logLM_all;

%% 4) 2 param LOG-level constraint optimization

outp = log(ACER);
% x = [b c]
if strcmp(penalty, 'Yes')
    F = @(x) exp(alpha.*(x(2).^sign(x(2) - 1) - 1)) * sum(W.*(outp.*outp)) - (sum(W.*outp))^2 - ...
        (sum(W.*((eta - x(1)).^x(2).*outp)) - ...
        sum(W.*(eta - x(1)).^x(2))*sum(W.*outp))^2/ ...
        (sum(W.*((eta - x(1)).^x(2).*(eta - x(1)).^x(2))) - ...
        (sum(W.*(eta - x(1)).^x(2)))^2);
else
    F = @(x) sum(W.*(outp.*outp)) - (sum(W.*outp))^2 - ...
        (sum(W.*((eta - x(1)).^x(2).*outp)) - ...
        sum(W.*(eta - x(1)).^x(2))*sum(W.*outp))^2/ ...
        (sum(W.*((eta - x(1)).^x(2).*(eta - x(1)).^x(2))) - ...
        (sum(W.*(eta - x(1)).^x(2)))^2);
end

ObFun_logNLINFIT = optimset('Display','off',...
    'Algorithm', 'interior-point', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol4, sol(4,5)] = fmincon(F,[b0 c0],[],[],[],[],...
    [bmin 1], [eta_1 Inf], [], ObFun_logNLINFIT);

inp = (eta - sol4(1)).^sol4(2);
A = (sum(W.*(inp.*outp)) - sum(W.*inp)*sum(W.*outp))/(sum(W.*(inp.*inp)) - (sum(W.*inp))^2);
B = sum(W.*outp) - A*sum(W.*inp);
sol(4,1:4) = [exp(B) sol4(1) -A sol4(2)];
if  any(imag(sum(sol(4,:))))
    sol(4,:) = 100*ones(1,5);
end
clear F inp A B sol4 ObFun_logNLINFIT;

%% 5) fixed q, 2 param LOG-level LS constraint optimization

if strcmp(penalty, 'Yes')
    F = @(x) sqrt(exp(alpha.*(x(2).^sign(x(2) - 1) - 1))) * ...
        sqrt(W).*(outp - log(qMed) - ((eta - x(1)).^x(2))*...
        sum(W.*((eta - x(1)).^x(2)).*(outp - log(qMed)))/...
        sum(W.*((eta - x(1)).^x(2)).*((eta - x(1)).^x(2))));
else
    F = @(x) sqrt(W).*(outp - log(qMed) - ((eta - x(1)).^x(2))*...
        sum(W.*((eta - x(1)).^x(2)).*(outp - log(qMed)))/...
        sum(W.*((eta - x(1)).^x(2)).*((eta - x(1)).^x(2))));
end

ObFun_log_LS_NLINFIT = optimset('Display','off',...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol5, sol(5,5)] = lsqnonlin(F, [b0 c0], ...
    [bmin 1], [eta_1 Inf], ObFun_log_LS_NLINFIT);
inp = (eta - sol5(1)).^sol5(2);
A = sum(W.*inp.*(outp - log(qMed)))/sum(W.*(inp.^2));
sol(5,1:4) = [qMed sol5(1) -A sol5(2)];
if  any(imag(sum(sol(5,:))))
    sol(5,:) = 100*ones(1,5);
end
clear F inp A B sol5 ObFun_log_LS_NLINFIT;

%% 6) 2 param LOG-level LS LM optimization

if strcmp(penalty, 'Yes')
    F = @(x) sqrt(exp(alpha.*(x(2).^sign(x(2) - 1) - 1))) *...
        sqrt(W).*(outp - sum(W.*outp) - ...
        ((sum(W.*((eta - x(1)).^x(2).*outp)) - sum(W.*(eta - x(1)).^x(2))*sum(W.*outp))/...
        (sum(W.*((eta - x(1)).^x(2).*(eta - x(1)).^x(2))) - (sum(W.*(eta - x(1)).^x(2)))^2))*...
        ((eta - x(1)).^x(2) - sum(W.*(eta - x(1)).^x(2))));
else
    F = @(x) sqrt(W).*(outp - sum(W.*outp) - ...
        ((sum(W.*((eta - x(1)).^x(2).*outp)) - sum(W.*(eta - x(1)).^x(2))*sum(W.*outp))/...
        (sum(W.*((eta - x(1)).^x(2).*(eta - x(1)).^x(2))) - (sum(W.*(eta - x(1)).^x(2)))^2))*...
        ((eta - x(1)).^x(2) - sum(W.*(eta - x(1)).^x(2))));
end

ObFun_log_LS_NLINFIT_LM = optimset('Display','off',...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 10000, 'MaxIter', 10000,...
    'TolX', 1e-12, 'TolFun', 1e-12);
[sol6, sol(6,5)] = lsqnonlin(F, [b0 c0], [], [], ObFun_log_LS_NLINFIT_LM);
inp = (eta - sol6(1)).^sol6(2);
A = (sum(W.*(inp.*outp)) - sum(W.*inp)*sum(W.*outp))/(sum(W.*(inp.*inp)) - (sum(W.*inp))^2);
B = sum(W.*outp) - A*sum(W.*inp);
sol(6,1:4) = [exp(B) sol6(1) -A sol6(2)];
if  any(imag(sum(sol(6,:))))
    sol(6,:) = 100*ones(1,5);
end

clear F inp outp A B sol6 ObFun_log_LS_NLINFIT_LM;

%%
% dim_fails = size(sol(sol(:,4)>=5,:));
% sol(sol(:,4)>=5,:) = 100*ones(dim_fails);
good_sol = sol(sol(:,4)<5,:);
[fmin, pos] = min(good_sol(:,5));
fin_sol = good_sol(pos,1:4);
%%
super_sol = sol(sol(:,4)<=3.5,:);
super_fmin = min(super_sol(:,5));
if ~isempty(super_fmin)
    if (super_fmin - fmin)/super_fmin<0.5
        pos = (sol(:,5) == super_fmin);
        [~, pos] = max(pos);
        fin_sol = sol(pos,1:4);
        fmin = super_fmin;
    end
end

if isempty(fin_sol)
    [fmin, pos] = min(sol(:,5));
    fin_sol = sol(pos,1:4);
end


end

