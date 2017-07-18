### Alpha-Beta filter for detrended fluctuation analysis (DFA)

This is a simple Python function (`ab_filter`) that calculates
continuous variation in scaling exponents by applying alpha-beta
filter to the DFA results. This procedure was reported in Echeverria
et al. "Interpretation of heart rate variability via detrended
fluctuation analysis and ab filter", *Chaos* 2003; 13(2):467-75.

The main output of a DFA routine is window sizes and fluctuation
levels. When plotted in log-log they show linear tendencies, which are
summarized by a scaling exponent. However, sometimes one might
encounter phenomena that deviate from truly linear relationship, for
example, cross-over is to be described by two linear fits and two
exponents.

Alpha-Beta filter takes as an input log-fluctuations and log-window
sizes to produce a *continous* plot of the scaling exponent variation
over window sizes.

Consult the help of the function `ab_filter` for further details on
usage and see information on alpha-beta filter elsewhere.

##### Package

A single principal function `ab_filter` constitutes the package (with
two ancillary functions). Dependency: python3 + NumPy.

