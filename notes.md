# Forecasting: Principles and Practice

Notes from reading through [Forecasting: Principles and Practice, 3rd Edition](https://otexts.com/fpp3/)

## __TODO__
Do all the exercises! Read https://r4ds.had.co.nz/index.html chapters 1-8 to get started.

Read these papers:
* Box-Cox transformation: https://www.jstor.org/stable/2984418?seq=1
* STL decomposition: https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/stl-a-seasonal-trend-decomposition-procedure-based-on-loess.pdf

## Chapter 1 - Getting Started

Some notation:
* $t$ specifies a timestep
* $y_t$ is the observation at timestep $t$
* $\mathcal{I}$ is the information we have observed
* $p(y_t|\mathcal{I})$ is the probability distribution that observation $y_t$ is drawn from, given what we know from $\mathcal{I}$
* Typically our estimate at time $t$ is the expected value of this distribution, denoted $\hat{y}_t=\mathbb{E}[y_t]_{p(y_t|\mathcal{I})}$; computationally we may prefer to take the mean or median of a number of predictions
* A useful shorthand to express what information was used in arriving at a prediction is $\hat{y}_{t+h|t-1}$, which tells us that we arrived at the prediction at time $t+h$ using the observations $(y_1, y_2, ..., y_{t-1})$.


## Chapter 2 - Time Series Graphics

Time series are typically made up of a number of patterns:
* __Trend__ - a long-term, continual increase/decrease in values
* __Seasonal__ - the impact due to factors such as time of the year, day of the week, etc; these have a fixed and known period
* __Cyclic__ - rises and falls caused by external factors that do not have a typical frequency; this may be influenced by something like the unemployment rate
* __Noise__ - random (or currently unexplained) perturbations that cause deviation from the underlying components

We can use scatterplots to visualise how correlated timeseries are (duh). We can also use these to compare time series to their own lags (previous values at a fixed point) to see if any "autocorrelation" exists.

Autocorrelation measures how correlated a timeseries is with its own previous points. For example, the measure of how correlated a timeseries is with the points $k$ timesteps back is:
$$r_k = \frac{\sum_{t=k+1}^{T}(y_t-\bar{y})(y_{t-k}-\bar{y})}{\sum_{t=1}^{T}(y_t-\bar{y})^2}$$

A significant trend will be evided in a plot of lags - values from a similar period will be similar in value - this may seen in combination with any seasonal effects.

If there is no statistically significant autocorrelation present in a timeseries then we are looking at white noise.


## Chapter 3 - Time Series Decomposition

Common adjustments include:
* __Calendar__ - artefacts of the calendar might induce variations, e.g. monthly totals will be affected by how many days there are in a month, we should consider daily average per month instead
* __Population__ - populations also change over time, it is worth keeping this in mind when exploring changes in e.g. hospital beds - it's more important to think about whether the number of beds _per capita_ is changing
* __Inflation__ - the value of a unit of money changes over time, if looking at long term prices it is worth adjusting by a price index so we can look at everything e.g. in its equivalent 2000 value; e.g. if the original price in year $t$ is $p_t$ and the price index has value $z_t$ then we will adjust everthing to be $p'_t = p_t/z_t \times z_{2000}$
* __Mathematical__ - we can use the usual mathematical transformations to alter the shape of a signal so that it's easier to model, e.g. if variations get larger as the magnitude of the signal increases then a Box-Cox transformation can help to standardise the variation across the full trace

If we assume time series are made up of various components, we can reason about how these components combine to give us the resulting signal. There are two main options:
* __Additive__ - this is the simplest option, at time $t$ we say our observation is given by $y_t = S_t + T_t + R_t$, where $S_t$ is the seasonal component, $T_t$ is a combination of the trend and cyclical components and $R_t$ is the remainder
* __Multiplicative__ - alternatively we can assume that the product of each component gives us the final value; this may be evident in cases where the size of variation depends on the height of the signal, the resulting observation is given by $y_t = S_t\times T_t\times R_t$. It is possible to transform the original signal to a point where it can be instead explained using additive modelling, for example by taking logs

Often we are interesting in subtracting/dividing out the seasonal variation to see a bigger-picture impact. For example, unemployment naturally varies with the season, but we may be more interested in seeing the general trend if we are in a recession.


The moving average forms the basis of extracting a trend-cycle component in the classical decomposition of a time series. It's given by:
$$\hat{T}_t = \frac{1}{m}\sum_{i=t-k}^{t+k}y_{i}$$
Where $m=2k+1$ is the size of the window.

We can combine together moving averages for different purposes. Ideally we want the window to be symmetric, so $m$ has to be an odd number, but we can get around this by doing a $2\times m$ moving average when $m$ is even, e.g. a $2\times4\mathrm{-MA}$:
$$
\begin{align}
\hat{T}_t &= \frac{1}{2}[\frac{1}{4}(y_{t-2} + y_{t-1} + y_{t} + y_{t+1}) + \frac{1}{4}(y_{t-1} + y_t + y_{t+1} + y_{t+2})]\\
&= \frac{1}{8}y_{t-2} + \frac{1}{4}y_{t-1} + \frac{1}{4}y_t + \frac{1}{4}y_{t+1} + \frac{1}{8}y_{t+2}
\end{align}
$$
This gives the first and last points of the window a slightly smaller weight. In the case of quarterly samples this would mean we could the first quarter twice (once in each subsequent year), but the total weighting is eqivalent to what all of the other quarters recieve.

We can achieve a better smoothing by weighting entries in the moving average:
$$\hat{T}_t = \frac{1}{m}\sum_{i=t-k}^{t+k}a_iy_{i}$$
Where $a_i = a_{-i}$ and all of the $a_i$ sum to 1. This allows entries to drop in and out of the moving average gradually, rather than at full weight.


### Classical Decomposition

For data with a seasonal period $m$ (e.g. 4 for quarterly data or 12 for monthly data with period of a year, 7 for daily data with period of a week), the classical additive decomposition is:

1. __Trend-cycle component__: $\hat{T}_t$ is obtained using a moving average; if $m$ is even then a $2\times m\textrm{-MA}$ is taken, otherwise a standard $m\textrm{-MA}$ suffices.
2. __Seasonal component__: We first find the trend-adjusted series $y'_t = y_t - \hat{T}_t$, then for each "season" in the time series (e.g. if we had monthly data, for each January, February, March, ...) we find the average value of the time-adjusted series, then roll back out to give us a value for each $t$ in our original time series; this is $\hat{S}_t$
3. __Remainder component__: $R_t = y_t - \hat{T}^t - \hat{S}_t$

For the multiplicative decomposition we just replace the subtractions with divisions.

This approach isn't recommended - there are many better ones. Some problems:
* The approach to determining the trend-cycle component means we lose the first and last set of values when modelling
* The trend-cycle component can over-smooth even regular sharp changes
* The approach to determining seasonality assumes this component is the same for every period
* Large outliers can impact seasonality and/or trend

### X11

This is largely based on the classical decomposition but includes algorthmic approaches which solve many of the above proble. For example, we don't lose any data in the trend-cycle component, it is better suited to sudden changes in signal, etc. _This isn't available in `statsmodels`_.

### SEATS

Seasonal Extraction in ARIMA Time Series (SEATS) is another alternative method, but it requires quarterly or monthly data. _This isn't available in `statsmodels`_.

### STL

Seasonal and Trend decomposition using Loess (STL) has lots of advantages over the above methods:
* Can handle any kind of seasonality
* The seasonal component can change over time
* There is control over the smoothness of the trend-cycle component
* It can be made robust to outliers (but these will appear in the remainder component)

Buuuut it can only be used for additive decomposition; we can get around this with Box-Cox transformations.

The key parameters to play with are:
* __`trend: int`__ - the number of consecutive observations used when estimating the trend-cycle component; the smaller this number is the more rapidly this component is allowed to change
* __`seasonal: int`__ - the number of consecutive periods (e.g. years) over which to determine the seasonal component


## Chapter 4 - Time Series Features

There are lots of common features to be derived from time series data. These include:
* __Simple statistics__: things such as min, first quartile, median, third quartile, max, ...
* __ACL__: the first autocorrelation coefficient, the sum-squared of the first 10 autocorrelation coefficients (useful if we want to know if there is _some_ autocorrelation, but when we don't know at which lag), the first autocorrelation coefficient of the diffs, ...


### STL

Recall that for additive decompositions we write our original series $y_t$ in the form:
$$
y_t = T_t + S_t + R_t
$$
Where $T_t$ is the trend-cycle component, $S_t$ the seasonal component and $R_t$ the remained. For data with a strong trend, we expect $T_t + R_t$ (the seasonally adjusted component) to have much more variation than $R_t$ alone (the latter will be distributed around zero, the former will follow an underlying trend all over the place). In contrast, the values should be roughly the same for data with little to no trend. Thus, we can attempt to quantify how "trendy" the data is by looking at the quantity:
$$
F_T = \max\Bigg(0, 1-\frac{\mathrm{Var}(R_t)}{\mathrm{Var}(T_t + R_t)}\Bigg)
$$

We can do the same to look at the strength of the seasonal component, this time instead using the de-trended series:
$$
F_S = \max\Bigg(0, 1-\frac{\mathrm{Var}(R_t)}{\mathrm{Var}(S_t + R_t)}\Bigg)
$$

A scatter plot (e.g.) can be used to identify which series exhibit the strongest seasonal component, strongest trend and so on. There are also a lot of other useful features to look at, e.g.:
* __Seasonal peak/trough__: for e.g. quarterly data, this tells us which quartier has the highest/lowest value, for instance a ski resort may peak in winter and have a trough in summer.
* ...


## Chapter 5 - The Forecaster's Toolbox

### Simple Forecasting Methods

#### Average Method

Here we just forecast using the mean of historical data:
$$
\hat{y}_{T+h|T} = \frac{1}{T}(y_1 + \ldots + y_T)
$$
Where $\hat{y}_{T+h|T}$ tells us we are predicting the value at timestamp $T+h$ (i.e. $h$ steps into the future) using the mean of all observations from $1$ to $T$.

#### Naive Method

This method is surprisingly effective for financial and economic forecasts - we just predict that any future value will be equal to the last observation, i.e.:
$$
\hat{y}_{T+h|T} = y_T
$$

This is best applied when the series we are forecasting is a random walk.

#### Seasonal Naive Method

This is similar to the above, but instead we predict that each future value will be the same as it was in the same season during the last period, e.g. for quarterly data we would say that for each Q1 we predict what the final Q1 in our training data was. For something like daily data we may have to specify whether the period is a week, month or year. This is a bit complicated to write out, but the concept is simple:
$$
\hat{y}_{T+h|T} = y_{T+h-m(k+1)}
$$
Where $m$ is the seasonal period (e.g. 12 for monthly data), $k$ is the integer part of $(h-1)/m$ - just how many complete periods have passed between the end of the training data and the timestamp we're forecasting.

#### Drift Method

This is very similar to the naive method, but allows us to roughly account for a long term trend. It boils down to predicting a straigh line from the very first observation to the very last, then projecting that forwards:
$$
\hat{y}_{T+h|T} = y_{T} + h\Bigg(\frac{y_T-y_1}{T-1}\Bigg)
$$

### Residuals

Residuals should really exhibit two patterns:
* __Uncorrelated__: if there are correlations then some information has not been picked up by the model
* __Zero-mean__: if the residuals aren't distributed about zero then there is some bias in the model

When it comes to prediction confidence intervals, it also helps if they satisfy these constraints:
* __Constant variance__: i.e. the variance of the distribution of residuals does not change over time
* __Normally distributed__

For example, consider predicting a closing stock price using the naive method - the residuals will just be the difference between successive values:
$$
e_t = y_t - \hat{y}_t = y_t - y_{t-1}
$$

We can use simple methods like studying the autocorrelation of each lag, or we can use more formal techniques such as the Box-Pierce test (simple) or the Ljung-Box test (preferred), both of which compare the distribution of autocorrelations of various lags. These tests are collectively known as portmanteau tests. For example, the Box-Pierce test uses the test statistic:
$$
Q = T\sum_{k=1}^lr_k^2
$$
Where $l$ is the maxmimum lag being tested (more on this in a second), $T$ is the number of observations, $r_k$ is the correlation of points with lag $k$. It is recommended to use $l=\min(10,T/5)$ for non-seasonal data or $l=\min(2m,T/5)$ for seasonal data. The $T/5$ term is just to account for the fact that the test performs poorly for large $l$. The Ljung-Box test is similarly defined, but a bit more complex.

In each case, for white noise the test statistics will have a $\chi^2$ distribution with $(l-K)$ degrees of freedom (where $K$ is the number of parameters in the model, e.g. the Naive method has no parameters, the drift method has one). We can thus use these statistics to perform a basic statistical test.

### Distributional forecasts

We can generally also provide a confidence interval about our prediction; for confidence level $c$ if we assume a normal distribution then we can write our prediction interval as:
$$
\hat{y}_{T+h|T}\pm c\hat{\sigma}_h
$$
Where $\hat{\sigma}_h$ is an estimate of the standard deviation at the $h$th forecast step.

For a one-step forecast, we can estimate the confidence intervals by looking at the standard deviation of the residuals across the training data. For a model which fits $K$ parameters, this is given by:
$$
\hat{\sigma} = \sqrt{\frac{1}{T-K}\sum_{t=1}^Te_t^2}
$$
There is no general solution for extrapolating this interval to multiple steps - the approach required depends on the form of modelling used. In general though, one must assume that the residuals are uncorrelated.

If we can't make the assumption of a normal distribution, one option is to look at estimating the confidence intervals by simulating lots of "futures" with _bootstrapped residuals_. This process assumes only constant variance over time and independence. For a one-step forecast the error is:
$$
e_t = y_t - \hat{y}_{t|t-1}
$$
Which we can rewrite as:
$$
y_t = \hat{y}_{t|t-1} + e_t
$$
We can then use this to simulate a future observation:
$$
y_{T+1} = \hat{y}_{T+1} + e_{T+1}
$$
We don't yet know what $e_{T+1}$ will be, so instead we can pick at random a residual value from our training dataset - this assumes that future errors will look like the ones we've seen before. We can repeat the process for $y_{T+2},\ldots,y_{T+h}$ to simulate a future time seires, and we can repeat the whole process to get an idea of what various possible futures would look like.

If we perform any transformations before modelling, when we transform back we may well find that our back-transformed point forecast is the median of the forecast distribution, rather than the mean. In order to recover the mean, we may have to bias-adjust our back-transformation as in many cases the mean is easier to work with than the median.

### Using Decompositions

Often it can be useful to separate season from the other components when modelling. Assuming an additive decomposition, we can write our time series as:
$$
y_t = T_t + S_t + R_t = S_t + A_t
$$
Where $A_t \equiv T_t + R_t$ is the _seasonal adjusted_ component. There's an equivalent definition for multiplicative decompositions. We can then forecast $\hat{S}_t$ and $\hat{A}_t$ separately. Often we use a very simple approach to forecasting the seasonal component, the naive seasonal method can suffice.

### Metrics

Some new observations:
* Minimising MAE leads to forecasts of the median, RMSE leads to forecasts of the mean.
* MAPE is only useful if zero is a meaningful value - this is not true e.g. in the case of Farenheit and Celsius where the zero-point is arbitrary
* sMAPE solves the problem of positive-negative asymmetry in MAPE, but can return negative values and remains unstable about zero - it is therefore discouraged
* We can define our error by comparing the value to a baseline model - if we divide our errors by the mean training error for the basline model then a value greater than one tells us our model is underperforming the baseline, a value less than one tells us we're outperforming it


#### For distributions

For a specified quantile (e.g. for $p=0.10$ we are looking at the 10th quantile, or the band defined by the 90% confidence interval), we can use the quantile score penalise forecasts more heavily which reach outside the band.

The Winkler score allows us to score the prediction interval directly, rather than looking at a series of quantiles. The score is a basic combination of two factors - if a given point lies wihtin the interval then the score is just the length of the interval at that point, but if then point lies outside the interval then we add a penalty. Consider a prediction interval with confidence $100(1-\alpha)%$ at time $t$, let it be defined by its limits $[l_{\alpha,t},u_{\alpha,t}]$; for a given observation $y_t$ the Winkler score is defined:
$$
W_{\alpha,t} = 
\begin{cases}
    (u_{\alpha,t} - l_{\alpha,t}) + \frac{2}{\alpha}(l_{\alpha,t} - y_t),& y_t < l_{\alpha,t} \\
    (u_{\alpha,t} - l_{\alpha,t}),& l_{\alpha,t}\leq y_t\leq u_{\alpha,t} \\
    (u_{\alpha,t} - l_{\alpha,t}) + \frac{2}{\alpha}(y_t - u_{\alpha,t}),& y_t > u_{\alpha,t}
\end{cases}
$$

If we care about the whole distribution instead of specific quantiles/intervals, we can use the Continuous Ranked Probability Score (CRPS), which averages the quantile score over all values of $p$.

We can use the skill score metric to generate a scale free comparison, e.g. for comparing how much better/worse method $M_a$ is to $M_b$ when using the CRPS metric, we look at the value of:
$$
\frac{\mathrm{CRPS}(M_a) - \mathrm{CRPS}(M_b)}{\mathrm{CRPS}(M_a)}
$$
If $M_b$ outperforms $M_a$ then this value will be positive.