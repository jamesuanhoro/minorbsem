# Modify major parameters table helper function

A function that adds user friendly descriptions to the major parameters
table

## Usage

``` r
modify_major_params(
  major_parameters,
  idxs,
  group = "",
  op = "",
  from = "",
  to = ""
)
```

## Arguments

- major_parameters:

  Major paramters table

- idxs:

  Relevant rows indexes

- group:

  Parameter group

- op:

  lavaan style operator

- from:

  Variable from

- to:

  Variable to

## Value

Updated major parameters table
