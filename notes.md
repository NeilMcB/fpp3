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
