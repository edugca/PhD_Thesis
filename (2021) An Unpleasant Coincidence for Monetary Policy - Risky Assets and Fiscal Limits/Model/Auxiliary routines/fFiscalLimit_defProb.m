function prob = fFiscalLimit_defProb(bValue, x1,x2,x3,x4,x5)

prob = normcdf(bValue, ...
                    fFiscalLimit_mu(x1,x2,x3,x4,x5), ...
                    fFiscalLimit_std(x1,x2,x3,x4,x5));
                            
end