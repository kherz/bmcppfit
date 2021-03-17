# Documentation

The *bmcppfit* application runs a Bloch-McConnell simulation on the .seq-file and tries to minimize the difference between the *fit_data* and the simulated data. Fit parameters can be choosen by the user.

## -p
The fit **p**arameter file.

The .yaml-file structure is similar to the on in the Pulseq-CEST project. See the example [file](../tests/matlab/yaml_fit.yaml) to see how to define them. The following parameters can be set:

- water_pool: f, t1 and t2
- mt_pool (optional): f, t1, t2, k, dw and lineshape
- cest_pool(s) (optional): f, t1, t2, k, dw
- b0: field strength [T]
- rel_b1 (optional): relative B1 (1 = 100 %)
- scale(optional): Initial magnetitazion at end of readout [0 ... 1]
- max_pulse_samples (optional): max number of pulseq pulse samples for simulation
- fit_data: The Data that should be fitted
- pulseq_file: pulseq file containng the sequence events
    - The number of ADC events must be equal to the number of entries in the fit_data
- weights (optional): weight vector if some offsets should be weighted more than others
    - The number of weights must be equal to the number of entries in the fit_data. The residuals at each data point will be multiplied by the weight at the same position. 
- threads: number of openMP threads
- fit_params: the parameters that should be fitted. Every parameter has a start value as well as an upper and lower bound. The following parameters can be fitted:
    - b0shift: B0 field shift (ppm)
    - scale
    - water_t1: t1 of water pool
    - water_t2: t2 of water pool
    - cest_n_t1: t1 of cest pool n
    - cest_n_t2: t2 of cest pool n
    - cest_n_k: exchange rate [Hz] of cest pool n
    - cest_n_dw: shift [ppm] of cest pool n
    - cest_n_f: fraction of cest pool n

If you have an experiment with two cest pools: 
```
cest_pool:
  amide: {f: 0.00064865, t1: 1.3, t2: 0.1, k: 30.0, dw: 3.5}
  NOE_1: {f: 0.0045, t1: 1.3, t2: 0.005, k: 16.0, dw: -3.5}
```
And you want to fit the exchange rate af the amide pool and the fraction of the NOE pool (and the T2 of the water pool), you set the following fit parameters in the parameter file:

```
fit_params:
  water_t2: {start: 0.0675, upper: 0.1, lower: 0.01}
  cest_1_k: {start: 60.0, upper: 1000.0, lower: 5.0}
  cest_2_f: {start: 0.0036, upper: 0.0225, lower: 0.00045}
 ```

 ## -f
 The **f**it algorithm options file.

 The non-linear least squares fit uses the [ceres-solver](http://ceres-solver.org/).
The -f option is optional, if no settings are provided, the standard ceres non-linear least square fit parameters are used.

From the [ceres options](http://ceres-solver.org/nnls_solving.html#solver-options), the following parameters can be changed:
- trust_region_strategy_type
- max_num_iterations
- max_solver_time_in_seconds
- initial_trust_region_radius
- max_trust_region_radius
- min_trust_region_radius
- function_tolerance 
- gradient_tolerance
- parameter_tolerance

For Bloch-McConnell fitting it can be very helpful to lower the function tolreance as residuals are usually quite small. 