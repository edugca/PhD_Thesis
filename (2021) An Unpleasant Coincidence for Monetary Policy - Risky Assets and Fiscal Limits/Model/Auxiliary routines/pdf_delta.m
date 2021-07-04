function result = pdf_delta(epsA, epsR, epsD, epsG, epsBeta)

simStruct = struct();
simStruct.params  = evalin('base', 'paramsStruct');

fiscalLimit_mu  = evalin('base', 'fiscalLimit_mu');
fiscalLimit_std = evalin('base', 'fiscalLimit_std');

result = pickFiscalLimit(simStruct, fiscalLimit_mu, fiscalLimit_std, ...
                epsA, epsR, epsD, epsG, epsBeta);
                
end