# Experiment Code

Supplement to: Kidney Exchange with Inhomogeneous Edge Existence Uncertainty (Hoda Bidkhori, John P Dickerson, Duncan C McElfresh, Ke Ren), UAI20.

These experiments test various robust and non-robust kidney exchange matching algorithms. Written for Pyhton 3, and requires Gurobi.

Example usage:

```
python -m robust_kex_experiment  --num-trials 1 --num-realizations 10 --alpha 0.5  --num-measurements 3  --use-omniscient --use-nonrobust --use-nonrobust-variable-p --use-nonrobust-fixed-p --num-measurements=5  --output-dir ./output_directory --input-dir ./intput_directory
```

Will produce an output file that looks like this. The first line contains experiment parameters, and each remaining line contains results for a single method on each graph:

```
Namespace(chain_cap=4, cycle_cap=3, input_dir='./intput_directory/', num_weight_measurements=5, num_weight_realizations=10, output_dir='./output_directory/',...)
   graph_name, trial_num, realization_num, cycle_cap, chain_cap, method, realized_score, runtime
unos_bimodal_apd_v64_i3_maxcard.input, 0, 0, 3, 4, nonrobust, 0E+00, 0.060791730880737300
unos_bimodal_apd_v64_i3_maxcard.input, 0, 0, 3, 4, nonrobust_fixed_p, 2E+00, 0.0398869514465332
unos_bimodal_apd_v64_i3_maxcard.input, 0, 0, 3, 4, nonrobust_variable_p, 0E+00, 0.04274106025695800
unos_bimodal_apd_v64_i3_maxcard.input, 0, 0, 3, 4, saa, 3E+00, 0.07482099533081060
unos_bimodal_apd_v64_i3_maxcard.input, 0, 0, 3, 4, omniscient, 6E+00, 0.0
unos_bimodal_apd_v64_i3_maxcard.input, 0, 1, 3, 4, nonrobust, 7E+00, 0.060791730880737300
unos_bimodal_apd_v64_i3_maxcard.input, 0, 1, 3, 4, nonrobust_fixed_p, 4E+00, 0.0398869514465332
unos_bimodal_apd_v64_i3_maxcard.input, 0, 1, 3, 4, nonrobust_variable_p, 7E+00, 0.04274106025695800
unos_bimodal_apd_v64_i3_maxcard.input, 0, 1, 3, 4, saa, 0E+00, 0.07482099533081060
...
```