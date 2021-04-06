syms c_0 c_1 y_0 y_1 kk_m1 kk_0 kk_1 kk_2 ppi_0 ppi_1 nr_m1 nr_0 nr_1 g_1 xxiM_0 ;

eqns = [
    - c_0 == - c_1 + nr_0 - ppi_1 ;
    - c_0 == - c_1 + g_1 + rSS*( c_1 + (1+eeta)/(1-aalpha)*y_1 - (1+aalpha*eeta)/(1-aalpha)*kk_1 ) ;
    ppi_0 == ppsi*( (eeta+aalpha)/(1-aalpha)*y_0 - aalpha*(1 + eeta)/(1 - aalpha)*kk_0 + c_0 ) + bbeta*ppi_1 ;
    y_0 == cSS/ySS*c_0 + kkSS/ySS*kk_1 - (1-ddelta)*kkSS/ySS*kk_0 ;
];

% Iterate N times
for i = 1:10
    eqns = subs(eqns, c_0, a0*kk_0 + a1*xxiM_0 + a2*nr_m1) ;
    eqns = subs(eqns, c_1, a0*kk_1 + a1*rrhoM*xxiM_0 + a2*nr_0) ;

    eqns = subs(eqns, ppi_0, b0*kk_0    + b1*xxiM_0          + b2*nr_m1) ;
    eqns = subs(eqns, ppi_1, b0*kk_1    + b1*rrhoM*xxiM_0    + b2*nr_0) ;

    eqns = subs(eqns, y_0, d0*kk_0    + d1*xxiM_0          + d2*nr_m1) ;
    eqns = subs(eqns, y_1, d0*kk_1    + d1*rrhoM*xxiM_0    + d2*nr_0) ;

    eqns = subs(eqns, kk_1, f0*kk_0    + f1*xxiM_0        + f2*nr_m1) ;
    eqns = subs(eqns, kk_2, f0*kk_1    + f1*rrhoM*xxiM_0  + f2*nr_0) ;

    eqns = subs(eqns, g_1, kkappaLine*(kk_2 - kk_1) - kkappaLine*(kk_1 - kk_0) ) ;
    eqns = subs(eqns, nr_0, rrhoNr*nr_m1 + (1-rrhoNr)*nnu*ppi_0 + xxiM_0) ;
end

% Isolate coefficients
eqns_0 = subs(eqns, [xxiM_0, nr_m1], [0, 0]);
eqns_1 = subs(eqns, [kk_0, nr_m1], [0, 0]);
eqns_2 = subs(eqns, [xxiM_0, kk_0], [0, 0]);

eqns_0 = expand( subs(eqns_0, kk_0, 1) );
eqns_1 = expand( subs(eqns_1, xxiM_0, 1) );
eqns_2 = expand( subs(eqns_2, nr_m1, 1) );

% eqns_0 = subs(eqns_0, b2^2, 0);
% eqns_1 = subs(eqns_1, b2^2, 0);
% eqns_2 = subs(eqns_2, b2^2, 0);