var Y C L W RK K I LAMBDA Q MC MT PI i R LRR;
varexo eps_m;
parameters  BETA ETA DELTA ALPHA NU PSI RHOM HABIT KAPPA OMEGA; 

load PARAMFILE;
set_param_value('BETA',BETA);
set_param_value('ETA',ETA);
set_param_value('NU',NU);
set_param_value('DELTA',DELTA);
set_param_value('ALPHA',ALPHA);
set_param_value('RHOM',RHOM);
set_param_value('PSI',PSI);
set_param_value('OMEGA',OMEGA);
set_param_value('HABIT',HABIT);
set_param_value('KAPPA',KAPPA);


model(linear);
// Some local variables
#YBAR = 1;
#KBAR = (ALPHA/((1/BETA)-1+DELTA));
#IBAR = DELTA*KBAR;
#CBAR = YBAR-IBAR;
#RKBAR = (1/BETA)-1+DELTA;

// Model equations, see Appendix C
(1-HABIT)*(1-BETA*HABIT)*LAMBDA = HABIT*C(-1) - (1+BETA*HABIT)*C + BETA*HABIT*C(+1);
LAMBDA + W = (ETA/(1-ALPHA))*Y - ((ALPHA*ETA)/(1-ALPHA))*K(-1);
LAMBDA = LAMBDA(+1) + i - PI(+1);
L = (1/(1-ALPHA))*Y - (ALPHA/(1-ALPHA))*K(-1);
RK = RKBAR*(W-K(-1)+L);
MC = W + ((ALPHA)/(1-ALPHA))*Y - ((ALPHA)/(1-ALPHA))*K(-1);
PI = BETA*PI(+1) + PSI*MC;
i = NU*PI + MT;
MT = RHOM*MT(-1) + eps_m;
Y = (CBAR/YBAR)*C + (IBAR/YBAR)*I;
K = (1-DELTA)*K(-1) + DELTA*I;
Q = OMEGA*(1+BETA)*I - OMEGA*I(-1) - BETA*OMEGA*I(+1);
Q + KAPPA*KBAR*(K-K(-1)) = LAMBDA(+1) - LAMBDA + RK(+1) + (1-DELTA)*Q(+1) + KAPPA*KBAR*(K(+1)-K);

// Define real interest rate
R = i - PI(+1);
// Long run real interest rate
-C = LRR;
end;


steady;
shocks;
var eps_m; stderr 0.01;
end;

stoch_simul(order=1, irf=40, nomoments, nograph) K Y C PI MC i R I MT LRR RK;