function finalfigure(level_of_interest, fin_sol, fin_sol_pl, fin_sol_mn,...
    bar_lev, acer_hat, CI_plus_reanch, CI_minus_reanch, k_memory, choice_of_ACER, sigma)

x_star = @(level,x) x(2)+(1/x(3)*(log(x(1))-log(level)))^(1/x(4));
epsilon = @(eta,x) x(1)*exp(-x(3)*(eta-x(2)).^x(4));

x_eps_star = x_star(level_of_interest,fin_sol);
x_eps_L = x_star(level_of_interest,fin_sol_mn);
x_eps_U = x_star(level_of_interest,fin_sol_pl);

fprintf('%4.4g & (%4.4g, %4.4g) \\\\ \r\n',...
    x_eps_star, x_eps_L, x_eps_U);
fprintf('[q b a c] = [%6.4g, %6.4g, %6.4g, %6.4g] \\\\ \r\n',...
    fin_sol);

x_eps = x_star(level_of_interest*1e-1,fin_sol);
x_CI_L = x_star(level_of_interest*1e-1,fin_sol_mn);
x_CI_U = x_star(level_of_interest*1e-1,fin_sol_pl);

BL_eps = linspace(bar_lev(1),x_eps,150);
BL_CI_L = linspace(bar_lev(1),x_CI_L,150);
BL_CI_U = linspace(bar_lev(1),x_CI_U,150);

figure
clf
% warning off all;
semilogy(bar_lev,(CI_plus_reanch),'--k')
hold on
semilogy(bar_lev,(CI_minus_reanch),'--k')
semilogy(bar_lev, acer_hat,'*k','MarkerSize',2)
semilogy(BL_eps,(epsilon(BL_eps,fin_sol)),'k')
semilogy(BL_CI_L,(epsilon(BL_CI_L, fin_sol_mn)),':k')
semilogy(BL_CI_U,(epsilon(BL_CI_U,fin_sol_pl)),':k')
semilogy(BL_CI_U,(level_of_interest)*ones(size(BL_CI_U)),':k')
semilogy(x_eps_star,level_of_interest,'*k', 'Markersize',7)

% semilogy(1/sigma*[x_eps_star x_eps_star],[level_of_interest level_of_interest*1e-1],'k')
%     txtar =  annotation('textarrow',[x_eps_star x_eps_star], [level_of_interest level_of_interest*1e-1],...
%         'string','We are here.');

%     legend('CI^{+}','CI^{-}','\epsilon_k(\eta)','\epsilon_k^{fit}','CI^{+}_{extr}','CI^{-}_{extr}')
%      xlabel('\eta^{100yr}','Units','normalized','Position',[0.98 -0.03],'FontWeight','bold','FontSize',12)
%      ylabel(['ACER_1(\eta)'],'Rotation',0,'Units','normalized','Position',[-0.04 0.95],'FontWeight','bold','FontSize',12)
xlabel('\eta')
ylabel(['ACER_{' num2str(k_memory(choice_of_ACER)) '}(\eta)'])
text(x_eps_star+0.04, level_of_interest*2e-1, horzcat(num2str(x_eps_star,'%4.3g\n')),...
    'HorizontalAlignment','Center','FontName', 'Times New Roman', 'FontSize', 16)
editplot

end
