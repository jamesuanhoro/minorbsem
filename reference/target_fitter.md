# target fitter function

A function that takes user input and fits the Stan model.

## Usage

``` r
target_fitter(
  data_list,
  seed,
  warmup,
  sampling,
  refresh,
  adapt_delta,
  max_treedepth,
  chains,
  ncores,
  show_messages,
  pa = FALSE
)
```

## Arguments

- data_list:

  Data list object passed to Stan

- seed:

  (positive integer) seed, set to obtain replicable results.

- warmup:

  (positive integer) The number of warmup iterations to run per chain.

- sampling:

  (positive integer) The number of post-warmup iterations to run per
  chain, retained for inference.

- refresh:

  (positive integer) How often to print the status of the sampler.

- adapt_delta:

  (real in (0, 1)) Increase to resolve divergent transitions.

- max_treedepth:

  (positive integer) Increase to resolve problems with maximum tree
  depth.

- chains:

  (positive integer) The number of Markov chains to run.

- ncores:

  (positive integer) The number of chains to run in parallel.

- show_messages:

  (Logical) If TRUE, show messages from Stan sampler, if FALSE, hide
  messages.

- pa:

  (LOGICAL) If TRUE: Path-analytic model; If FALSE (default): Generic
  model.

## Value

Fitted Stan model
