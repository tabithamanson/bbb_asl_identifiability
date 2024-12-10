# bbb_asl_identifiability

## Functions that solve the BBB ASL models:
- my_buxton_analytical.m                     
- my_solve_buxton_numerical.m                 
- my_solve_parallel_2CXM_numerical.m
- my_solve_series_2CXM_numerical.m
- my_solve_parallel_2CXM_T2_numerical.m      (multi-echo)
- my_solve_series_2CXM_T2_numerical.m        (multi-echo)

## Scripts for performing ID analysis:
- numerical_structural_id_by_sensitivity_SE_buxton.m
- numerical_structural_id_by_sensitivity_SE_parallel.m
- numerical_structural_id_by_sensitivity_SE_series.m
- numerical_structural_id_by_sensitivity_ME_parallel.m
- numerical_structural_id_by_sensitivity_ME_series.m

## Scripts for fitting simulations:
- sim_fit_ME_parallel_2CXM.m
- sim_fit_ME_series_2CXM.m

## kw_fitting_repository.m
Contains
- gm_whole_fit.m (fit for kw, ME parallel 2CXM)
- gm_whole_fit_series.m (fit for kw, ME series 2CXM)
- 4 x my_solve...m model functions (note these differ from those of the same name under "Functions that solve the BBB ASL models" - sorry about the naming)
- GM_tSNR_calculate_mcf_topup.m (used to find the relationship between noise and signal magnitude for one participant with 10 repeats)

