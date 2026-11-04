# Shell validation against calculix_linear

Shell mesh density factor 0.5. Target: each case within its tolerance (3 %) of the shell model.

| case | compared | 1d-fem | shell | difference | passed |
|---|---|---|---|---|---|
| channel_ft | flexural buckling factor | 1.1207e+05 | 1.115e+05 | +0.51 % | yes |
| channel_ft | torsional buckling factor | 1.4758e+05 | 1.4632e+05 | +0.86 % | yes |
| ibeam_ltb_buckle | torsional buckling factor | 34204 | 33902 | +0.89 % | yes |
| ibeam_ltb_path | v at x/L = 0.5, final level | 0.0014432 | 0.0014453 | 0.14 % of max | yes |
| ibeam_ltb_path | phi at x/L = 0.5, final level | 0.0063865 | 0.0064388 | 0.81 % of max | yes |
| ibeam_large_twist | phi at x/L = 1, final level | 0.86124 | 0.86509 | 0.44 % of max | yes |
| composite_box | u at x/L = 1, final level | -8.2667e-19 | 4.2569e-08 | 0.01 % of max | yes |
| composite_box | w at x/L = 1, final level | 1.2 | 1.1954 | 0.39 % of max | yes |
| composite_box | phi at x/L = 1, final level | -0.087254 | -0.087048 | 0.24 % of max | yes |

Path differences are the largest difference over the load levels, relative to the largest shell value of that quantity (at least 1e-3 rad for a twist and 1e-4 L for a displacement).
