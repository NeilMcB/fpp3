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


## Chapter 7 - Time Series Regression Models

A linear regression model is defined:
$$
y_t = \sum_{i=0}^{K}\beta_ix_{i,t} + \epsilon_t
$$
Where we predict the observation $y_t$ at time $t$ using $K$ predictor observations $x_{i,t}$ (with the special case of $x_0=1$, whose influence is measured with coefficients $\beta_i$. $\epsilon_t$ captures any random devations from our perfectly linear model. Each coefficient measures the marginal effect of its corresponding predictor - that is the impact that predictor has once all other predictors have been accounted for.

There are a number of assumptions implicit in a linear model:
* The relationship between predictor variables and the target is truly linear
* The errors have zero mean (otherwise our model is biased)
* The errors are not (auto)correlated, otherwise we have not made use of all available information
* Each predictor is not random, otherwise there is no meaningful input to derive from them

If we are using the Least Squares method then the end goal is to optimise the coefficients $\beta_0, \ldots, \beta_k$ such that we minimise the "loss" defined by:
$$
\sum_t\epsilon_t^2 = \sum_t(y_t - \beta_0 - \beta_1x_{1,t} - \ldots - \beta_Kx_{K,t})^2
$$

Once we have our fitted coefficients we can generate predictions:
$$
\hat{y}_t = \sum_{i=0}^K\hat{\beta}_ix_{i,t}
$$
These are taken within the training set, rather than being genuine forecasts.

The coefficient of determination gives us a measure for how well a linear model fits the data - this is the square of the correlation between the observed and predicted values, or can be calculated as:
$$
R^2 = \frac{\sum(\hat{y}_t - \bar{y})^2}{\sum(y_t - \bar{y})^2}
$$
This can be interpreted as the amount of variance present in the observed data that is accounted for in the predicted data. For a good model we would expect this value to be close to 1. This doesn't really tell us much about the performance of a forecast, though.

We can alternatively look at the standard error of the regression on the residuals. Note that if we have $K$ predictors then we fit $K+1$ parameters (a coefficient for each predictor, plus an intercept); this gives us the formula:
$$
\hat{sigma}_e = \sqrt{\frac{1}{T-K-1}\sum_{t=1}^Te_t^2}
$$
This is roughly the average size of the error that the model produces. We can then use this to generate prediction intervals.


The residuals for linear models have some interesting properties:
$$
\sum_{t=1}^Te_t=0 \quad \mathrm{and} \quad \sum_{t=1}^Tx_{k,t}e_t=0 \quad \forall k \in K
$$

For time series data, linear models may exhibit some autocorrelation - observations are likely to take on values similar to what they were in the past. This doesn't mean the linear model is incorrect, but it does highlight that maybe we aren't taking advantage of all the information available to us.

We can also study the quality of a fit by looking at things like scatterplots of residuals versus predictors, residuals versus the target and so on. This may tell us whether the relationship is non-linear, if variance depends on target value, etc. Transformations may be required if any of these are present.


A time series is stationary if its values fluctuate around a common mean over time with constant variance. If we do not work with stationary data then we can make spurious regressions - for example if both house price in London and rice production in Ghana happen to be rising, a regression model may conclude that there is an inter-dependence between the two. High $R^2$ and strong residual autocorrelation may indicate such spurious regression has taken place.


We can add trend to a regression model by just adding $t$ as a variable, e.g.:
$$
y_t = \beta_0 + \beta_1t + \epsilon_t
$$
We can add indicator variables in the case of binary data. These may be particularly useful when accounting for external events which may lead to outliers - e.g. if we're looking at data on tourism volume in Brazil we may wish to add an Olympics or World Cup indicator.


Intervention variables are a generalisation of this concept - these help us to process cases where there is a sudden change in a series behaviour. If the effect only lasts for one period then we can use a "spike" variable to tells us the period where the event occurs, or if the effect is permanent we can use a "step" variable - before the timestep in question its value is zero then thereafter it is one.

For things like advertising spend which may take a while to take effect, we can introduce lagged variables for predictors too.

For seasonal data we can use a Fourier decomposition to generate features.

In order to select the predictors to use, we can just pick the combination which optimise our forecasting power. Common metrics to do this include (Corrected) AIC and Time Series Cross Validation (TSCV). For a small number of predictors we just try out all combinations and choose the best one. For a large number of predictors we can perform either backwards (forwards) stepwise regression, where we either start with all (none) of the predictors then iteratively remove (add) predictors and see whether model performance increases.

Ex-ante forecasts use only information that is available at the time the forecast is made. For example, if we are basing our regression of external predictors then we must also forecast these at the same time we forecast our target. Ex-post forecasts assume we're able to know the predictors ahead of time - the only thing we don't know about is the target. It can be useful to compare ex-ante to ex-post forecasts to see if any drop in performance is due to the model itself, or the method we use to forecast the predictors. In a small number of cases we can know predictors ahead of time, e.g. calendar features, trend features. One possible approach is to perform scenario testing; this is where we explore the impact on our model that differenct scenarios would have - in this case we assume we know perfectly what future values will be.

An easy way to get around this problem is to instead base our regression model on lagged predictors - if we care about a forcast 8 weeks into the future then regress on the predictor 8 weeks ago; this way when we need the forecast we can use today's value for our model.

log-log models give us an interpretation of the coefficients where $\beta_k$ provides an estimate of the elasticity of predictor $k$. If our dataset contains zeros we can instead transform using $\phi = \log(x+1)$; helpfully this means that we'll still have zeros in our transformed space where we have them in the original. We can use all the usual techniques for other forms of non-linear modelling: arbitrary feature transformations (of both predictors and timesteps for trends). We can also add piecewise non-linearities, e.g. by introducting predictors that only take on a value after a certain timestep/true value, e.g. if we let $x_{1,t} = x_t$ and define:
$$
x_{2,t} = (x-c)_+ = 
\begin{cases}
0,&x<c \\
x-c,&x\geq c
\end{cases}
$$
Basically this lets us introduce an additional coefficient $\beta_2$ if the value of $x$ crosses some threshold $c$.

This approach is often preferable for non-linear trends - extrapolating quadratic and above trends can get a bit out of hand for distanct predictions.


## Chapter 8 - Exponential Smoothing

Exponential smoothing applies successively smaller weights to prior observations (in an expoentially decreasing manner) - it gives us something in between the naive method and a simple mean. For a smoothing parameter $0\leq\alpha\leq1$, we can generate a forecast from:
$$
\hat{y}_{T+1|T} = \alpha y_T + \alpha(1-\alpha)y_{T-1} + \alpha(1-\alpha)^2y_{T-2} + \ldots
$$
The sum of all these weights will approximate one. If $\alpha$ is close to zero then a larger weight will be assigned to distant observations; the larger it gets the more emphasis we place on recent observations instead.

Simple Exponential Smoothing (SES) forecasts are said to be flat; for any horizon $h$, the value assigned will be the last level component:
$$
\hat{y}_{T+h|T} = \hat{y}_{T+1|T} = l_T, \quad h=1,2,3,\ldots
$$
Where:
$$
l_T = \alpha y_T + (1-\alpha)l_{T-1} = \alpha y_T + \alpha(1-\alpha) y_{T-1} + (1-\alpha)l_{T-2} = \ldots
$$

The choices of $\alpha$ and $l_0$ (the recursive base) can be chosen by minimising error (e.g. SSE) over the training data.


We can then incorporate trend into our exponential smoothing:
$$
\begin{align}
    &\mathrm{Forecast}&\quad \hat{y}_{t+h|t}&= l_t + hb_t\\
    &\mathrm{Level}&\quad l_t&= \alpha y_t + (1-\alpha)(l_{t-1} + b_{t-1})\\
    &\mathrm{Trend}&\quad b_t&= \beta^*(l_t-l_{t-1})+(1-\beta^*)b_{t-1}
\end{align}
$$
We introduce a new smoothing term $0\leq\beta^*\leq1$, so we use a smoothed estimate of both the level and the slope to generate our final estimate. This means that to go $h$ steps into the future we just take our level estimate and add $h$ times our estimate of the slope.

These basic trended forecasts tend to over-trend, so we can instead dampen the trend overtime to arrive at a flat slope:
$$
\begin{align}
    &\mathrm{Forecast}&\quad \hat{y}_{t+h|t}&= l_t + \Big(\sum_{i=1}^h\phi^i\Big)b_t\\
    &\mathrm{Level}&\quad l_t&= \alpha y_t + (1-\alpha)(l_{t-1} + \phi b_{t-1})\\
    &\mathrm{Trend}&\quad b_t&= \beta^*(l_t-l_{t-1})+(1-\beta^*)\phi b_{t-1}
\end{align}
$$
Where, again, $0\leq\phi\leq1$; the smaller the value the greater the damping effect.

We can also add a seasonal factor:
$$
\begin{align}
    &\mathrm{Forecast}&\quad \hat{y}_{t+h|t}&= l_t + \Big(\sum_{i=1}^h\phi^i\Big)b_t + s_{t+h-m(k+1)}\\
    &\mathrm{Level}&\quad l_t&= \alpha(y_t - s_{t-m}) + (1-\alpha)(l_{t-1} + \phi b_{t-1})\\
    &\mathrm{Trend}&\quad b_t&= \beta^*(l_t-l_{t-1})+(1-\beta^*)\phi b_{t-1}\\
    &\mathrm{Season}&\quad s_t&= \gamma(y_t-l_{t-1}-b_{t-1}) + (1-\gamma)s_{t-m}
\end{align}
$$
Where $m$ is the length of the season (e.g. 12 for monthly data). This time the coefficient is defined in the range $0\leq\gamma\leq1-\alpha$. There may be some incorrect parentheses etc going on in the above formulae.

The models we construct when generating exponentially smoothed forecasts are also knonw as state space models since we don't just output our forecast over time but also keep track of some underlying variables: level, trend and seasonality.


We can think about how model predictions are updated by looking at the "error correction" form of the forecast equation:
$$
\begin{align}   
    &\mathrm{Forecast}&\quad \hat{y}_{t+1|t}&=l_t\\
    &\mathrm{Smoothing}&\quad l_t&=\alpha y_t + (1-\alpha)l_{t-1}
\end{align}
$$
We can rearrange this smoothing equation:
$$
\begin{align}
    l_t&=l_{t-1}+\alpha(y_t - l_{t-1})\\
       &=l_{t-1}+\alpha(y_t - \hat{y}_{t|t-1})\\
       &=l_{t-1}+\alpha e_t
\end{align}
$$
In other words, as we step through our training data, we look at how good our previous forecast was then update our next forecast accordingly, with a weighting set by $\alpha$. If we over-predicted in the last timestep then we adjust by subtracting proportionately from the next step.

We can also use these concepts to build a generative model of the data. There are two stages to the process. Firstly, we assume there is some underlying latent model of the data ($l_t$) which evolves according to:
$$
l_t = l_{t-1} + \alpha\epsilon_t
$$
Where $\epsilon_t \sim \mathcal{N}(0,\sigma^2)$ is just white noise. We assume the same white noise is present when we actually observe this latent model too:
$$
y_t = l_t + \epsilon_t
$$
It's kind of like a Markov Chain model, but the influence of distant events is dampened out by the exponential smoothing.


Alternatively we can think about errors in a multiplicative sense; if we write out the error in its relative form:
$$
\epsilon_t = \frac{y_t - \hat{y}_{t|t-1}}{\hat{y}_{t|t-1}}
$$
Then we can make the substitution $l_{t-1} = \hat{y}_{t|t-1}$ and rearrange as before:
$$
y_t = l_{t-1}(1 + \epsilon_t)
$$
And similarly:
$$
l_t = l_{t-1}(1 + \alpha\epsilon_t)
$$

We can repeat this process for more complex exponential smoothing methods that include trend and season components.


## Chapter 9 - ARIMA Models

A timeseries is stationary if its statistical properties do not depend on the time at which we observe it - for instance, its mean and variance should not change with time. For example, data with seasonality and/or trend are not stationary. Data with a cyclic component may still be classed as stationary because this may not necessarily be regular. Strictly speaking, if ${y_t}$ is stationary then for all $s$ the distribution of $(y_t,\ldots,y_{t+s})$ does not depend on $t$.

We can try to stabilitise variance by taking logarithmss and mean by differencing.

A differenced model is denoted:
$$
y_t' = y_t - y_{t-1}
$$
If this results in a random walk dataset (i.e. stationary), we can denote this:
$$
y_t' = y_t - y_{t-1} = \epsilon_t
$$
From which we can recover the random walk equation:
$$
y_t = y_{t-1} + \epsilon_t
$$
Where $\epsilon\sim\mathcal{N}(0,\sigma^2)$ is white noise. We can introduce a non-zero mean (i.e. $\epsilon\sim\mathcal{N}(\mu,\sigma^2)$) which will give our model a propensity to move up or down with a given magnitude.


Sometimes we have to take second-order differences (but rarely do we go beyond):
$$
\begin{align}
    y_t'' &= y_t' - y_{t-1}'\\
    &= (y_t - y_{t-1}) - (y_{t-1} - y_{t-2})\\
    &= y_t - 2y_{t-1} + y_{t-2}
\end{align}
$$


We can also take seasonal differences - these compare the difference between the current value and the value at the same point last season:
$$
y_t' = y_t - y_{t-m}
$$
If this seasonally differenced data looks like white noise then we can model with:
$$
y_t' = y_{t-m} + \epsilon_t
$$

We can do both seasonal and ordinary differencing. The order doesn't really matter, but it can help to do seasonal differencing first then check if we need to difference any further.

There are a few statistical tests we can use to check for both stationarity and seasonal stationarity. We can also test to see how many differences are needed.


We can use the backshift operator to denote lagged observations, so:
$$
\begin{align}
    By_t &= y_{t-1}\\
    B(By_t) &= B^2y_t = By_{t-2}
\end{align}
$$
A useful example for monthly data is the backshift operator that gives us the same observation this time last year:
$$
B^{12}y_t = y_{t-12}
$$

We can then rewrite first order differences as:
$$
\begin{align}
    y_t' &= y_t - y_{t-1}
         &= y_t - By_t
         &= (1-B)y_t
\end{align}
$$
In general, the $d$th order difference can be written as:
$$
(1-B)^dy_t
$$

These operators are particularly useful when combining differences, e.g. for a first order difference combined with a seasonal difference of length $m$:
$$
y_t' = (1-B)(1-B^m)y_t
$$


### Autoregression

For time series, we can construct a linear model by regressing on past values of the variable we're interested in. An autoregressive model of order $p$ takes the form:
$$
y_t = c + \sum_{k=1}^p\phi_kB^ky_t + \epsilon_t
$$
Where $\epsilon_t$ is white noise. This can more concisely be denoted as an $\mathrm{AR}(p)$ model.

### Moving Average

Moving average models instead look at the prediction error over time, e.g. an $\mathrm{MA}(q)$ model looks like:
$$
y_t = c + \epsilon_t + \sum_{k=1}^q\theta_kB^k\epsilon_t
$$

### Non-Seasonal ARIMA

We can simply combine autoregressive and moving average models to give us the (non-seasonal) AutoRegressive Integrated Moving Average (ARIMA) model. This looks like:
$$
y_t'=c+\sum_{i=1}^p\phi_iB^iy_t'+\sum_{j=1}^q\theta_jB^j\epsilon_j
$$
We can summarise these models with three parameters; an $\mathrm{ARIMA}(p,d,q)$ model is defined by:
* $p$: order of the autoregressive part;
* $d$: degree of the first differencing performed;
* $q$: order of the moving average part.

These models require the data to be stationary (amongst a few other things).

We can generate periodic forecasts if $p\geq2$ and some other constraints within the model parameters are met.


We can estimate the appropriate value of $p$ by looking at how autocorrelated our time series is. In practice, if we have a significant first degree autocorrelation then we'll likely have a large second degree autocorrelation too - if 2 is strongly influenced by 1 and 3 by 2 then it'll apear as if 3 is also strongly influenced by 1. The equivalent holds for other multiples. Instead we can look at partial autocorrelations - each degree is adjusted to show the effects at that degree only. For example, to understand the correlation present between $y_t$ and $B^ky_t$ we have to remove the effects of the lags at $y_{t-1},\ldots,y_{t-k+1}$.

We can use the ACF and PACF plots to help us find a decent ARIMA model. For an $\mathrm{ARIMA}(p,d,0)$ model, the differenced data may exhibit:
* an exponentially decaying or sinusoidal ACF;
* a significant spike at lag $p$, but none beyond.
For an $\mathrm{ARIMA}(0,d,q)$ model, the differenced data may exhibit:
* an exponentially decaying or sinusoidal PACF;
* a significant spike at lag $q$, but none beyond.

Parameter tuning ($c,\phi_1,\ldots,\phi_p,\theta_1,\theta_q$) is done by performing Maximum Likelihood Estimation (MLE). Hyperparameter tuning is done by selecting the choice of $p$ and $q$ that minimise AIC, AICc or BIC. __Note__ that these scores at not comparible across different choices of $d$ and so this should be picked using other statistical measures in advance.

The general outline for fitting a non-seasonal ARIMA is:
1. Plot the data and identify any outliers
2. Check if the variance needs to be stabilised - use something like a Box-Cox transformation if so
3. If the data are non-stationary, take differences until the data are stationary (use a statistical test to check)
4. Examine the (P)ACF: is an $\mathrm{ARIMA}(p,d,0)$ or $\mathrm{ARIMA}(0,d,q)$ model appropriate?
5. Try the chosen model, then play around with the options to find the one which minimises AICc
6. Check the residuals of the fitted model - use a Portmanteau test to verify they look like white noise; if that fails then we need to modify the model again
7. __Done!__ Use the model to make some sickkkk forecasts

In order to generate a forecast upto time $y_{T+h}$ where we have observations upto time $T$, we rearrange the ARIMA function to have $y_t$ on the LHS. We then generate predictions for each of $\hat{y}_{T+1|T},\ldots,\hat{y}_{T+h|T}$ through the following process:
1. For $\hat{y}_{T+1|T}$, the only unknown is $\epsilon_{T+1}$ - we can just replace this with zero and we have our estimate
2. For $\hat{y}_{T+2|T}$, our unknowns are $\epsilon_{T+2}$, (maybe $\epsilon_{T+1}$ if $q$ is large enough) and $y_{T+1}$. The first two we replace with zero as before, and the last one we just sub in the estimate we generated first.
3. Repeat this cycle until we reach our desired forecast horizon.

We can also generate probability distributions over our forecast estimates.


### Seasonal ARIMA

For a seasonal ARIMA we can just add a second set of terms, so we have the non-seasonal parameters $(p,d,q)$ and the seasonal parameters $(P,D,Q)_m$, where as usual $m$ is the seasonal period.


## Chapter 10 - Dynamic Regression

In this model formulation we drop the assumption that regression errors are uncorrelated white noise; instead we allow the errors to be due to their own ARIMA model. So the original formulation:
$$
y_t = \sum_{i=0}^K\beta_kx_{k,t} + \epsilon_t
$$
Becomes, in the case of an $\mathrm{ARIMA}(1,1,1)$ error model, the following:
$$
\begin{align}
    y_t &= \sum_{i=0}^K\beta_kx_{k,t} + \eta_t
    (1 - \phi_1B)(1 - B)\eta_t = (1 + \theta_1B)\epsilon_t
\end{align}
$$
In other words, we still have our white noise term, but it's wrapped up in a more complex model.

All of our predictors still have to be stationary in this case for their coefficients to be meaningful. If any one of our variables is not stationary, it's common to difference all of them so as to keep each predictor at the same degree.

For timeseries with seasons of a large period, it can be more effective to perform a dynamic regression with Fourier components as predictors.