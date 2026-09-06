# PRA — Probabilistic Risk Assessment using L-moments

MATLAB code for estimating failure probabilities when the available samples are scarce and
may contain extremes. Distribution parameters are estimated by **L-moments**, which are far
less sensitive to extreme values than conventional moment estimators, and the fitted
distributions are propagated through surrogate models to count failures across operating
temperatures.

![Probabilistic risk assessment result](https://user-images.githubusercontent.com/64420389/226306620-dbd11f14-acc8-4030-bd68-46c6b359587f.png)

## Files

| File | What it does |
| --- | --- |
| `Identify_parent_distribution.m` | Selects the parent distribution for a sample. |
| `HIC_RBF_model.m`, `Seventeen_design_variable_RBF.m` | RBF surrogates — single-response, and the 17-design-variable case. |
| `RND_test.m`, `RND_test_17variable.m` | Random-sample propagation through the surrogates. |
| `onetemp_failures_sample.m`, `All_temperature_failures_sample.m` | Failure counts at a single temperature and across all temperatures. |
| `No_of_failures_1above.m` | Failure counts above a threshold. |
| `Plot_all_temp_comparison.m` | Comparison of results across temperatures. |
| `Box_16.m`, `Box_median.m`, `Area_shade_plot.m` | Result plots. |

## Requirements

MATLAB, with the SURROGATES Toolbox for the RBF models. The `.mat` sample and result files
the scripts load are not tracked here.
