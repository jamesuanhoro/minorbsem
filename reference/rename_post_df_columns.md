# Rename columns of posterior data.frame prior by parameter type

Rename columns of posterior data.frame prior by parameter type

## Usage

``` r
rename_post_df_columns(
  df,
  labels_1,
  labels_2,
  begin_name,
  operation,
  search_replace,
  search_term_1,
  search_term_2
)
```

## Arguments

- df:

  Posterior data.frame

- labels_1:

  Labels for first variable

- labels_2:

  Labels for second variable

- begin_name:

  String to begin columns names with

- operation:

  Operation, one of =~, ~, \~~

- search_replace:

  TRUE/FALSE for type of column naming

- search_term_1:

  Terms to use in replacement for first variable

- search_term_2:

  Terms to use in replacement for second variable

## Value

A renamed posterior data.frame
