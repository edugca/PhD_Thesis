% Brault, Joshua and Khan, Hashmat (2019). "The real interest rate channel
% is structural in contemporary New-Keynesian models: A Note". Journal of
% Money, Credit and Banking.

% Running the file will generate figures 1-4 in the paper. 
% KAPPA controls the CAC and OMEGA controls the IAC.
BETA = 0.99;
ETA = 1;
DELTA = 0.025;
ALPHA = 0.3;
NU = 1.5;
THETA = 0.7;
EPS = 0.83;
PSI = ((1-BETA*THETA)*(1-THETA)/THETA)*((1-ALPHA)/(1-ALPHA+(ALPHA/(1-EPS))));
RHOM = 0.95;
HABIT = 0;
KAPPA = 0;
OMEGA = 0;
save PARAMFILE BETA ETA NU DELTA ALPHA EPS THETA RHOM PSI OMEGA HABIT KAPPA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock persistence vector
rhomvec = [0.5; 0.95];
for i = 1:2
    OMEGA = 5.48; KAPPA = 0;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    cvec(:,i) = C_eps_m*100; invvec(:,i) = I_eps_m*100;
    rvec(:,i) = R_eps_m*100; yvec(:,i) = Y_eps_m*100;
end

figure(1);
subplot(2,2,1)
plot(cvec(:,1), 'b', 'LineWidth',2); hold on
plot(cvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Consumption'}); ylim([-1.75 0.25])
subplot(2,2,2)
plot(invvec(:,1), 'b', 'LineWidth',2); hold on
plot(invvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Investment'}); ylim([-4.30 0.25])
subplot(2,2,3)
plot(rvec(:,1), 'b', 'LineWidth',2); hold on
plot(rvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Real interest rate'}); ylim([-0.1 1])
subplot(2,2,4)
plot(yvec(:,1), 'b', 'LineWidth',2); hold on
plot(yvec(:,2), 'r--', 'LineWidth',2); hold on
legend('\rho = 0.5','\rho = 0.95', 'Location', 'southeast')
grid on; title({'Output'}); ylim([-2 0.5])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock persistence vector
rhomvec = [0.5; 0.95];
for i = 1:2
    OMEGA = 2.5; KAPPA = 0;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    cvec(:,i) = C_eps_m*100; invvec(:,i) = I_eps_m*100;
    rvec(:,i) = R_eps_m*100; yvec(:,i) = Y_eps_m*100;
end

figure(2);
subplot(2,2,1)
plot(cvec(:,1), 'b', 'LineWidth',2); hold on
plot(cvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Consumption'}); ylim([-1.75 0.25])
subplot(2,2,2)
plot(invvec(:,1), 'b', 'LineWidth',2); hold on
plot(invvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Investment'}); ylim([-5.25 0.25])
subplot(2,2,3)
plot(rvec(:,1), 'b', 'LineWidth',2); hold on
plot(rvec(:,2), 'r--', 'LineWidth',2); hold on
grid on; title({'Real interest rate'}); ylim([-0.1 1])
subplot(2,2,4)
plot(yvec(:,1), 'b', 'LineWidth',2); hold on
plot(yvec(:,2), 'r--', 'LineWidth',2); hold on
legend('\rho = 0.5','\rho = 0.95', 'Location', 'southeast')
grid on; title({'Output'}); ylim([-2 0.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock persistence vector
rhomvec = [0; 0.5; 0.95];
for i = 1:3
    OMEGA = 0; KAPPA = 0;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    cvec(:,i) = C_eps_m*100; invvec(:,i) = I_eps_m*100;
    rvec(:,i) = R_eps_m*100; yvec(:,i) = Y_eps_m*100;
end

figure(3);
subplot(2,2,1)
plot(cvec(:,1), 'b', 'LineWidth',2); hold on
plot(cvec(:,2), 'b--', 'LineWidth',2); hold on
plot(cvec(:,3), 'r-.', 'LineWidth',2); hold on
grid on; title({'Consumption'}); ylim([-1.75 0.75])
subplot(2,2,2)
plot(invvec(:,1), 'b', 'LineWidth',2); hold on
plot(invvec(:,2), 'b--', 'LineWidth',2); hold on
plot(invvec(:,3), 'r-.', 'LineWidth',2); hold on
grid on; title({'Investment'}); ylim([-60 5])
subplot(2,2,3)
plot(rvec(:,1), 'b', 'LineWidth',2); hold on
plot(rvec(:,2), 'b--', 'LineWidth',2); hold on
plot(rvec(:,3), 'r-.', 'LineWidth',2); hold on
grid on; title({'Real interest rate'}); ylim([-0.6 0.15])
subplot(2,2,4)
plot(yvec(:,1), 'b', 'LineWidth',2); hold on
plot(yvec(:,2), 'b--', 'LineWidth',2); hold on
plot(yvec(:,3), 'r-.', 'LineWidth',2); hold on
legend('\rho = 0', '\rho = 0.5','\rho = 0.95', 'Location', 'southeast')
grid on; title({'Output'}); ylim([-13 2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock persistence vector
rhomvec = [0.5; 0.85];
for i = 1:2
    OMEGA = 0; KAPPA = 0;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    LRNOA(:,i) = LRR_eps_m*100; invNOA(:,i) = I_eps_m*100;
end
for i = 1:2
    OMEGA = 2.5; KAPPA = 0;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    LRIAC(:,i) = LRR_eps_m*100; invIAC(:,i) = I_eps_m*100;
end
for i = 1:2
    OMEGA = 0; KAPPA = 0.5;
    RHOM = rhomvec(i);
    save PARAMFILE BETA ETA NU DELTA ALPHA RHOM PSI OMEGA HABIT KAPPA
    dynare Brault_Khan_JMCB_2019 noclearall
    LRCAC(:,i) = LRR_eps_m*100; invCAC(:,i) = I_eps_m*100;
end


figure(4);
subplot(2,2,1)
plot(invNOA(:,1), 'b', 'LineWidth',2); hold on
plot(invCAC(:,1), 'b--', 'LineWidth',2); hold on
plot(invIAC(:,1), 'r-.', 'LineWidth',2); hold on
grid on; title({'\fontsize{16}\rho = 0.5', '\fontsize{10}Investment'}); ylim([-55 5])
subplot(2,2,2)
plot(invNOA(:,2), 'b', 'LineWidth',2); hold on
plot(invCAC(:,2), 'b--', 'LineWidth',2); hold on
plot(invIAC(:,2), 'r-.', 'LineWidth',2); hold on
grid on; title({'\fontsize{16}\rho = 0.85', '\fontsize{10}Investment'})
legend('No adj. costs', 'CAC','IAC', 'Location', 'southeast'); ylim([-35 5])
subplot(2,2,3)
plot(LRNOA(:,1), 'b', 'LineWidth',2); hold on
plot(LRCAC(:,1), 'b--', 'LineWidth',2); hold on
plot(LRIAC(:,1), 'r-.', 'LineWidth',2); hold on
grid on; title({'\fontsize{10}Long run real rate'}); ylim([-0.5 1.5])
subplot(2,2,4)
plot(LRNOA(:,2), 'b', 'LineWidth',2); hold on
plot(LRCAC(:,2), 'b--', 'LineWidth',2); hold on
plot(LRIAC(:,2), 'r-.', 'LineWidth',2); hold on
grid on; title({'\fontsize{10}Long run real rate'}); ylim([-0.5 2])




% doing some cleaning...
!del Brault_Khan_JMCB_2019.log
!del Brault_Khan_JMCB_2019_static.m
!del Brault_Khan_JMCB_2019.m
!del Brault_Khan_JMCB_2019_set_auxiliary_variables.m
!del Brault_Khan_JMCB_2019_results.mat
!del Brault_Khan_JMCB_2019_IRF_eps_m.eps
!del Brault_Khan_JMCB_2019_dynamic.m
!del PARAMFILE.mat


% End of file