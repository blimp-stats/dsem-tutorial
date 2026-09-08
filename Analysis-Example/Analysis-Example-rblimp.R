#-------------------------------------------------------------------#
# Dynamic Structural Equation Modeling in Blimp
#
# Supplemental Material
#
# Analysis Example in `rblimp`, including plots
#
# Analysis-Example-rblimp.R
#
# Copyright Brian T. Keller 2026 all rights reserved
#-------------------------------------------------------------------#
#
# Person-specific quantities:
#   b_0j = typical sleep quality       g_0j = typical pain
#   b_1j = sleep carry-over            g_1j = pain carry-over
#   b_2j = pain-to-sleep spillover     g_2j = sleep-to-pain spillover
#   l_j  = sleep log residual variance
#   w_j  = pain log residual variance
#   d_Sj = sleep linear trend          d_Pj = pain linear trend
#
# Blimp truncates variable names to 9 characters, so latent names
# must be kept at or below that length.
#
#-------------------------------------------------------------------#

#### Setup ####

# Load packages
library(rblimp)
library(ggplot2)

## Set ggplot2's global theme
theme_set(
    theme_bw(base_size = 12, base_family = 'serif')
    + theme(
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(
            color = 'black', linewidth = 0.5,
            lineend = 'square'
        )
    )
)

## Point rblimp at the build used for the manuscript
set_blimp('/Applications/Blimp/blimp-nightly')

#### Read data ####
PainDiary <- read.csv('paindiary.csv')

## Sampler settings shared by every model
SEED   <- 12371
BURN   <- 40000
ITER   <- 80000
NCHAIN <- 4


#-------------------------------------------------------------------#
#### Model 1: Random-Slope Univariate AR(1) ####
#
# Sleep quality on day t is a deviation around the person's own
# typical level, and the strength of the day-to-day carry-over is
# itself allowed to differ across people.
#-------------------------------------------------------------------#

m1 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j;',
    model = '
        between:
            b_0j ~~ b_1j;
            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;

        within:
            sleep_lag = sleep.lag - b_0j;
            sleep ~ intercept@b_0j sleep_lag@b_1j;
    ',
    parameters = c(
        'b_1_sd   = sqrt(b_1j.totalvar)',
        'b_1_low  = b_1 - b_1_sd',
        'b_1_high = b_1 + b_1_sd'
    ),
    simple = c(
        'sleep_lag | b_1j @ b_1_low and b_0j @ 0',
        'sleep_lag | b_1j @ b_1 and b_0j @ 0',
        'sleep_lag | b_1j @ b_1_high and b_0j @ 0'
    ),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m1)

## Raw Blimp output
output(m1)

## Posterior plot for the sleep equation
posterior_plot(m1, 'sleep')

#### Model 1 conditional effects ####

## Sleep carry-over at the average person-specific slope and one
## between-person SD below and above it
(
    simple_plot(sleep ~ sleep_lag | b_1j + b_0j, m1)
    + labs(
        title = 'Conditional Sleep Carry-Over',
        subtitle = 'Person-specific carry-over slope as moderator',
        x = 'Prior-Day Sleep Quality (Deviation From Trait Level)',
        y = 'Sleep Quality (Deviation From Trait Level)'
    )
)


#-------------------------------------------------------------------#
#### Model 2: Random-Slope Bivariate VAR(1) ####
#
# Sleep quality and pain are modeled jointly. Each within-person
# equation includes its own lagged deviation (carry-over) and the
# other outcome's lagged deviation (spillover). All six
# person-specific quantities are random and freely covary.
#-------------------------------------------------------------------#

m2 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j;

            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;
            b_2j ~ intercept@b_2;
            g_0j ~ intercept@g_0;
            g_1j ~ intercept@g_1;
            g_2j ~ intercept@g_2;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;
            pain ~~ sleep;
    ',
    parameters = c(
        'b_1_sd   = sqrt(b_1j.totalvar)',
        'b_1_low  = b_1 - b_1_sd',
        'b_1_high = b_1 + b_1_sd',
        'b_2_sd   = sqrt(b_2j.totalvar)',
        'b_2_low  = b_2 - b_2_sd',
        'b_2_high = b_2 + b_2_sd',
        'g_1_sd   = sqrt(g_1j.totalvar)',
        'g_1_low  = g_1 - g_1_sd',
        'g_1_high = g_1 + g_1_sd',
        'g_2_sd   = sqrt(g_2j.totalvar)',
        'g_2_low  = g_2 - g_2_sd',
        'g_2_high = g_2 + g_2_sd'
    ),
    simple = c(
        'sleep_lag | b_1j @ b_1_low and b_0j @ 0',
        'sleep_lag | b_1j @ b_1 and b_0j @ 0',
        'sleep_lag | b_1j @ b_1_high and b_0j @ 0',
        'pain_lag | b_2j @ b_2_low and b_0j @ 0',
        'pain_lag | b_2j @ b_2 and b_0j @ 0',
        'pain_lag | b_2j @ b_2_high and b_0j @ 0',
        'pain_lag | g_1j @ g_1_low and g_0j @ 0',
        'pain_lag | g_1j @ g_1 and g_0j @ 0',
        'pain_lag | g_1j @ g_1_high and g_0j @ 0',
        'sleep_lag | g_2j @ g_2_low and g_0j @ 0',
        'sleep_lag | g_2j @ g_2 and g_0j @ 0',
        'sleep_lag | g_2j @ g_2_high and g_0j @ 0'
    ),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m2)

#### Model 2 conditional effects ####

## Each lagged coefficient probed at the average of its own
## person-specific slope and one between-person SD either side

## Sleep carry-over
(
    simple_plot(sleep ~ sleep_lag | b_1j + b_0j, m2)
    + labs(
        title = 'Sleep Carry-Over',
        subtitle = 'Person-specific carry-over slope as moderator',
        x = 'Prior-Day Sleep Quality (Deviation From Trait Level)',
        y = 'Sleep Quality (Deviation From Trait Level)'
    )
)

## Pain-to-sleep spillover
(
    simple_plot(sleep ~ pain_lag | b_2j + b_0j, m2)
    + labs(
        title = 'Pain-to-Sleep Spillover',
        subtitle = 'Person-specific spillover slope as moderator',
        x = 'Prior-Day Pain (Deviation From Trait Level)',
        y = 'Sleep Quality (Deviation From Trait Level)'
    )
)

## Pain carry-over
(
    simple_plot(pain ~ pain_lag | g_1j + g_0j, m2)
    + labs(
        title = 'Pain Carry-Over',
        subtitle = 'Person-specific carry-over slope as moderator',
        x = 'Prior-Day Pain (Deviation From Trait Level)',
        y = 'Pain (Deviation From Trait Level)'
    )
)

## Sleep-to-pain spillover
(
    simple_plot(pain ~ sleep_lag | g_2j + g_0j, m2)
    + labs(
        title = 'Sleep-to-Pain Spillover',
        subtitle = 'Person-specific spillover slope as moderator',
        x = 'Prior-Day Sleep Quality (Deviation From Trait Level)',
        y = 'Pain (Deviation From Trait Level)'
    )
)


#-------------------------------------------------------------------#
#### Model 3: Moderated Sleep Carry-Over ####
#
# Adds the product of the two lagged within-person deviations to the
# sleep equation. Because both deviations are centered at latent
# trait levels, uncertainty in those levels carries into the product.
# The interaction coefficient, b_3, is fixed across people.
#-------------------------------------------------------------------#

m3 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j;

            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;
            b_2j ~ intercept@b_2;
            g_0j ~ intercept@g_0;
            g_1j ~ intercept@g_1;
            g_2j ~ intercept@g_2;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j
                (sleep_lag*pain_lag)@b_3;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;
            pain ~~ sleep;
    ',
    simple = 'sleep_lag | pain_lag and b_1j @ b_1 and b_0j @ 0',
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m3)

#### Model 3 conditional effects ####

## Simple slopes of sleep_lag at low/mean/high pain_lag
(
    simple_plot(sleep ~ sleep_lag | pain_lag + b_0j, m3)
    + labs(
        title = 'Conditional Sleep Carry-Over',
        subtitle = 'Prior-day relative pain as moderator',
        x = 'Prior-Day Sleep Quality (Deviation From Trait Level)',
        y = 'Sleep Quality (Deviation From Trait Level)'
    )
)

## Johnson-Neyman plot for the sleep-lag by pain-lag interaction
(
    jn_plot(sleep ~ sleep_lag | pain_lag + b_0j, m3)
    + labs(
        title = 'Johnson-Neyman Plot for Sleep Carry-Over',
        subtitle = 'Red area represents 0 within 95% interval',
        x = 'Prior-Day Pain (Deviation From Trait Level)',
        y = 'Sleep Carry-Over'
    )
)


#-------------------------------------------------------------------#
#### Technical supplement: person-specific linear trends ####
#
# Replaces each stable trait level with a person-specific linear
# trajectory, so the dynamic terms operate on deviations from that
# fitted trajectory (a linear residual DSEM).
#-------------------------------------------------------------------#

m3d <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j
                          d_Sj d_Pj;',
    model = '
        sleep_resid_lag = sleep.lag - (b_0j + d_Sj*(day - 1));
        pain_resid_lag  = pain.lag  - (g_0j + d_Pj*(day - 1));

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j d_Sj d_Pj ~~
                b_0j b_1j b_2j g_0j g_1j g_2j d_Sj d_Pj;

            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;
            b_2j ~ intercept@b_2;
            g_0j ~ intercept@g_0;
            g_1j ~ intercept@g_1;
            g_2j ~ intercept@g_2;
            d_Sj ~ intercept@d_S;
            d_Pj ~ intercept@d_P;

        within:
            sleep ~ intercept@b_0j day@d_Sj
                sleep_resid_lag@b_1j pain_resid_lag@b_2j
                (sleep_resid_lag*pain_resid_lag)@b_3;
            pain ~ intercept@g_0j day@d_Pj
                pain_resid_lag@g_1j sleep_resid_lag@g_2j;
            pain ~~ sleep;
    ',
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m3d)


#-------------------------------------------------------------------#
#### Model 4: Random-Slope Location-Scale VAR(1) ####
#
# Adds person-specific log residual variances for sleep and pain.
# Modeling the variances on the log scale keeps them positive and
# lets them covary with the six location quantities.
#-------------------------------------------------------------------#

m4 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;

            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;
            b_2j ~ intercept@b_2;
            g_0j ~ intercept@g_0;
            g_1j ~ intercept@g_1;
            g_2j ~ intercept@g_2;
            l_j  ~ intercept@l_0;
            w_j  ~ intercept@w_0;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j
                (sleep_lag*pain_lag)@b_3;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;

            var(sleep) ~ intercept@l_j;
            var(pain)  ~ intercept@w_j;

            pain ~~ sleep;
    ',
    parameters = c(
        'sleep_var_median = exp(l_0)',
        'sleep_var_mean   = exp(l_0 + l_j.totalvar/2)',
        'pain_var_median  = exp(w_0)',
        'pain_var_mean    = exp(w_0 + w_j.totalvar/2)'
    ),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m4)

## Posterior plots for the two scale equations
posterior_plot(m4, 'l_j')
posterior_plot(m4, 'w_j')


#-------------------------------------------------------------------#
#### Model 4 extension: a time-varying scale predictor ####
#
# Prior-day relative pain is added to the sleep variance equation,
# asking whether sleep becomes more or less predictable after an
# unusually painful day.
#-------------------------------------------------------------------#

m4t <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;

            b_0j ~ intercept@b_0;
            b_1j ~ intercept@b_1;
            b_2j ~ intercept@b_2;
            g_0j ~ intercept@g_0;
            g_1j ~ intercept@g_1;
            g_2j ~ intercept@g_2;
            l_j  ~ intercept@l_0;
            w_j  ~ intercept@w_0;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j
                (sleep_lag*pain_lag)@b_3;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;

            var(sleep) ~ intercept@l_j pain_lag@a;
            var(pain)  ~ intercept@w_j;

            pain ~~ sleep;
    ',
    parameters = c('var_ratio = exp(a)'),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m4t)


#-------------------------------------------------------------------#
#### Model 5: Stress Predicting the Person-Specific Quantities ####
#
# Grand-mean-centered perceived stress predicts both trait levels,
# all four lagged slopes, and both log residual variances. Stress was
# measured once per participant, so it distinguishes people rather
# than occasions.
#-------------------------------------------------------------------#

m5 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;
        stress_c  = stress - stress_mean;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;

            b_0j   ~ intercept@b_0 stress_c@b_4;
            b_1j   ~ intercept@b_1 stress_c@b_5;
            b_2j   ~ intercept@b_2 stress_c@b_6;
            g_0j   ~ intercept@g_0 stress_c@g_4;
            g_1j   ~ intercept@g_1 stress_c@g_5;
            g_2j   ~ intercept@g_2 stress_c@g_6;
            l_j    ~ intercept@l_0 stress_c@l_1;
            w_j    ~ intercept@w_0 stress_c@w_1;
            stress ~ intercept@stress_mean;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j
                (sleep_lag*pain_lag)@b_3;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;

            var(sleep) ~ intercept@l_j;
            var(pain)  ~ intercept@w_j;

            pain ~~ sleep;
    ',
    parameters = c(
        'sleep_var_ratio = exp(l_1)',
        'pain_var_ratio  = exp(w_1)'
    ),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m5)


#-------------------------------------------------------------------#
#### Model 6: Trait Volatility as a Predictor and Moderator ####
#
# Moves the two trait--volatility associations out of the level-2
# residual covariance matrix and into the trait-level regressions,
# and asks whether pain volatility moderates the stress--typical-pain
# association. The volatility residuals therefore form a covariance
# block separate from the six location residuals.
#-------------------------------------------------------------------#

m6 <- rblimp(
    data      = PainDiary,
    clusterid = 'person',
    timeid    = 'day',
    latent    = 'person = b_0j b_1j b_2j g_0j g_1j g_2j l_j w_j;',
    model = '
        sleep_lag = sleep.lag - b_0j;
        pain_lag  = pain.lag  - g_0j;
        stress_c  = stress - stress_mean;
        l_c       = l_j - l_0;
        w_c       = w_j - w_0;

        between:
            b_0j b_1j b_2j g_0j g_1j g_2j ~~
                b_0j b_1j b_2j g_0j g_1j g_2j;
            l_j w_j ~~ l_j w_j;

            b_0j ~ intercept@b_0 stress_c@b_4 l_c@b_7;
            g_0j ~ intercept@g_0 stress_c@g_4 w_c@g_7
                (stress_c*w_c)@g_8;

            b_1j   ~ intercept@b_1 stress_c@b_5;
            b_2j   ~ intercept@b_2 stress_c@b_6;
            g_1j   ~ intercept@g_1 stress_c@g_5;
            g_2j   ~ intercept@g_2 stress_c@g_6;
            l_j    ~ intercept@l_0 stress_c@l_1;
            w_j    ~ intercept@w_0 stress_c@w_1;
            stress ~ intercept@stress_mean;

        within:
            sleep ~ intercept@b_0j sleep_lag@b_1j pain_lag@b_2j
                (sleep_lag*pain_lag)@b_3;
            pain  ~ intercept@g_0j pain_lag@g_1j sleep_lag@g_2j;

            var(sleep) ~ intercept@l_j;
            var(pain)  ~ intercept@w_j;

            pain ~~ sleep;
    ',
    parameters = c(
        # log(1.5) = 0.405465 rescales the log-variance coefficients
        # to a 50% increase in residual variance
        'b_7_50 = b_7*log(1.5)',
        'g_7_50 = g_7*log(1.5)',
        'g_8_50 = g_8*log(1.5)',
        'w_sd   = sqrt(w_j.totalvar)',
        'w_low  = -w_sd',
        'w_high = w_sd'
    ),
    simple = c(
        'stress_c | w_c @ w_low',
        'stress_c | w_c @ 0',
        'stress_c | w_c @ w_high'
    ),
    seed = SEED, burn = BURN, iter = ITER, chains = NCHAIN
)

## Posterior summary
summary(m6)

#### Model 6 conditional effects ####

## Conditional stress association at low/mean/high pain volatility
(
    simple_plot(g_0j ~ stress_c | w_c, m6)
    + labs(
        title = 'Conditional Stress--Typical-Pain Association',
        subtitle = 'Latent pain volatility as moderator',
        x = 'Perceived Stress (Centered)',
        y = 'Typical Pain'
    )
)

## Johnson-Neyman plot for the stress by pain-volatility interaction
(
    jn_plot(g_0j ~ stress_c | w_c, m6)
    + labs(
        title = 'Johnson-Neyman Plot for the Stress Association',
        subtitle = 'Red area represents 0 within 95% interval',
        x = 'Latent Pain Volatility (Centered Log Residual Variance)',
        y = 'Typical Pain | Perceived Stress'
    )
)


#-------------------------------------------------------------------#
# End of file
#-------------------------------------------------------------------#
