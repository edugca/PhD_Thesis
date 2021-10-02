cpi = edu_Database_WB_InflationEMEDEV_Fetch('data', 'hcpi_a', 'CPI_');
target_upper = edu_Database_WB_InflationEMEDEV_Fetch('cc', 'it_upper', 'TUP_');

dataset = synchronize(cpi, target_upper);

plot(dataset.Date, [dataset.CPI_BRA, dataset.TUP_BRA])