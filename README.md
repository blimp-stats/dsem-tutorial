# Dynamic Structural Equation Modeling in Blimp

This repository is the supplemental material for the manuscript
"Dynamic Structural Equation Modeling in Blimp". The manuscript is in
preparation. To cite the manuscript, kindly use the following
citation:

Keller (2026). *Dynamic Structural Equation Modeling in Blimp*.
Manuscript in preparation.

This repository offers all the data and analysis code used in the
manuscript. Each model has a Blimp syntax file and the fitted output
it produced. Additionally, an R script is available that fits the
same models in `rblimp` and generates the conditional-effects
graphics, including the simple-slope and Johnson-Neyman plots.

`online-supplement.pdf` is the online supplement to the manuscript.
Section S1 explains the factored-regression representation, and
Section S2 presents a sensitivity analysis with person-specific
linear trends.

## Analysis-Example

The daily diary analysis, `paindiary.csv`, and every model in the
manuscript. Each model has a `.imp` syntax file and a matching
`.blimp-out` file.

| File | Description |
|---|---|
| `m1-random-ar1` | Random-slope univariate AR(1) for sleep quality |
| `m2-random-var1` | Random-slope bivariate VAR(1) for sleep quality and pain |
| `m3-moderated-carryover` | Moderation of sleep carry-over by prior-day relative pain |
| `m4-location-scale` | Model 3 plus person-specific log residual variances |
| `m4-time-varying-scale` | Prior-day relative pain predicting sleep residual variance |
| `m5-stress-prediction` | Stress predicting all eight person-specific quantities |
| `m6-trait-volatility` | Trait volatility as a predictor and as a moderator |
| `detrended-sensitivity` | Linear residual DSEM with person-specific trends |

`Analysis-Example-rblimp.R` fits all of the above in `rblimp` and
demonstrates `output()`, `summary()`, `posterior_plot()`,
`simple_plot()`, and `jn_plot()`.
