%% SOLVE UNDETERMINED COEFFICIENTS AND PLOT GRAPH

% assume(a0 >= 0 & a1 <= 0 ...
%         & b1 <=0 ...
%         & d0 >= 0 & d1 <= 0 ...
%         & f0 >= 0 & f1 <= 0);

optSymbolic = false;
assume([a0 a1 a2 b0 b1 b2 d0 d1 d2 f0 f1 f2 rrhoM], 'clear');

if optSymbolic

    assume(a0 > 0);
    assumeAlso(d0 > 0);
    if rrhoNr_val == 0
        assumeAlso(a2 == 0);
    else
        assumeAlso(a2 < 0);
        assumeAlso(b2 < 0);
        assumeAlso(d2 < 0);
    end
    eqns_num_0 = subs(eqns_0, [varVals{:,1}], vpa([varVals{:,2}]));
    eqns_num_2 = subs(eqns_2, [varVals{:,1}], vpa([varVals{:,2}]));
    soln_0 = solve([eqns_num_0; eqns_num_2], [a0, b0, d0, f0, a2, b2, d2, f2], 'ReturnConditions' , true, ...
        'IgnoreAnalyticConstraints', false);

    % Select the solution
    idxSoln_0 = 1;
    if length(soln_0.a0) == 1
        vec_soln_0 = vpa([soln_0.a0, soln_0.b0, soln_0.d0, soln_0.f0, soln_0.a2, soln_0.b2, soln_0.d2, soln_0.f2]);
    else
        vec_soln_0 = vpa([  soln_0.a0(idxSoln_0), soln_0.b0(idxSoln_0), soln_0.d0(idxSoln_0), soln_0.f0(idxSoln_0), ...
                            soln_0.a2(idxSoln_0), soln_0.b2(idxSoln_0), soln_0.d2(idxSoln_0), soln_0.f2(idxSoln_0)
                        ]);
    end

    eqns_num_1 = subs(eqns_1, [varVals{:,1}], vpa([varVals{:,2}]));
    eqns_num_1 = subs(eqns_num_1, [a0, b0, d0, f0, a2, b2, d2, f2], vec_soln_0);
    assumeAlso(rrhoM >= 0 & rrhoM <1);
    soln_1 = solve(eqns_num_1, [a1, b1, d1, f1], 'ReturnConditions' , true, ...
        'IgnoreAnalyticConstraints', false);

else

    eqns_num_0 = subs(eqns_0, [varVals{:,1}], vpa([varVals{:,2}]));
    eqns_num_2 = subs(eqns_2, [varVals{:,1}], vpa([varVals{:,2}]));
    %matlabFunction([eqns_num_0; eqns_num_2], 'File', 'f_eqns.m');
    
    % When you solve a system of equations with nonunique solutions,
    % the behavior of vpasolve depends on whether the system is polynomial
    % or nonpolynomial. If polynomial, vpasolve returns all solutions by
    % introducing an arbitrary parameter. If nonpolynomial, a single
    % numerical solution is returned, if it exists.
    limitNum = 1e5;
    limitMin = 1e-20;
    
    if isAlways(vpa(subs(rrhoNr, [varVals{:,1}], vpa([varVals{:,2}]))) == 0)
        valsRange = [   
            0 limitNum;
            -limitNum limitNum;
            -limitNum limitNum;
            -limitNum limitNum;
            0 limitMin;
            0 limitMin;
            0 limitMin;
            0 limitMin];
    else
        valsRange = [   
            0 limitNum;
            -limitNum limitNum;
            -limitNum limitNum;
            -limitNum limitNum;
            -limitNum 0;
            -limitNum limitNum;
            -limitNum 0;
            -limitNum 0];
    end
    
    
    soln_0 = vpasolve([eqns_num_0, eqns_num_2], ...
        [a0, b0, d0, f0, a2, b2, d2, f2], ...
        valsRange);
    
    vec_soln_0 = struct2array(soln_0);

    assumeAlso(rrhoM >= 0 & rrhoM <1);
    eqns_num_1 = subs(eqns_1, [varVals{:,1}], vpa([varVals{:,2}]));
    eqns_num_1 = subs(eqns_num_1, [a0, b0, d0, f0, a2, b2, d2, f2], vec_soln_0);
    assumeAlso(rrhoM >= 0 & rrhoM <1);
    soln_1 = solve(eqns_num_1, [a1, b1, d1, f1], 'ReturnConditions' , true, ...
        'IgnoreAnalyticConstraints', false);
    
end

% Select the solution
idxSoln_1 = 1;
if length(soln_1.a1) == 1
    vec_soln_1 = vpa([soln_1.a1, soln_1.b1, soln_1.d1, soln_1.f1]);
else
    vec_soln_1 = vpa([soln_1.a1(idxSoln_1), soln_1.b1(idxSoln_1), soln_1.d1(idxSoln_1), soln_1.f1(idxSoln_1)]);
end

% Plot a1
f = figure;
edu_GraphSetInterpreter('latex');
if length(soln_0.a0) == 1
    indEff = rrhoM*soln_1.a1 - soln_1.a1 + soln_0.a2*(1-rrhoNr)*nnu*soln_1.b1 + soln_0.a2;
    dirEff = soln_0.a0*soln_1.f1;
else
    indEff = rrhoM*soln_1.a1(idxSoln_1) - soln_1.a1(idxSoln_1) + soln_0.a2(idxSoln_0)*(1-rrhoNr)*nnu*soln_1.b1(idxSoln_1) + soln_0.a2(idxSoln_0);
    dirEff = soln_0.a0(idxSoln_0) * soln_1.f1(idxSoln_1) ;
end

indEff = subs(indEff, [varVals{:,1}], vpa([varVals{:,2}]));
dirEff = subs(dirEff, [varVals{:,1}], vpa([varVals{:,2}]));


p_1 = fplot(dirEff, [0,1], 'DisplayName', 'Direct effect');
p_1.Color = 'blue';
p_1.LineWidth = 1.5;
p_1.LineStyle = '-.';
hold on;
p_2 = fplot(indEff, [0,1], 'DisplayName', 'Indirect effect');
p_2.Color = 'red';
p_2.LineWidth = 1.5;
p_2.LineStyle = ':';
hold on;
p_3 = fplot(dirEff + indEff, [0,1], 'DisplayName', 'Total effect');
p_3.Color = 'black';
p_3.LineWidth = 1.5;
edu_GraphDrawZeroAxis(p_3);
title('$\hat{R}_t$: coefficient of $\xi^m_t$')
xlabel('$\rho^m$');
legend('Location', 'best');
set(gca, 'FontSize', 12);

nameGraph = ['R_coeff_xiM_' num2str(kkappa_val) '_' num2str(rrhoNr_val) '_' num2str(kkSS_val) '.png'];
set(f, 'Position',  [100, 0, 600, 400]); % resize figure
exportgraphics(f, [pathImages filesep nameGraph]);