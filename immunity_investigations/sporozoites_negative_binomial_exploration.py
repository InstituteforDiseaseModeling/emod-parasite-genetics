import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom
import itertools
from math import sqrt
import statsmodels
## constraints mean of 5, 80% CDF < 25, 6% CDF > 100


plist = np.logspace(-9,-0.000001,1000)

likelihoods = []
zinb = statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialP()
for i,p in enumerate(plist):

    n = 5*(1-p)/p
    mean, var, skew, kurt = nbinom.stats(n, p,moments='mvsk')
    cdf_1 = nbinom.cdf(k=25,n=n,p=p)
    cdf_2 = 1-nbinom.cdf(k=100,n=n,p=p)

    if 3.0<mean<6.0:
        if 0.7 < cdf_1 < 0.95:
            if 0 < cdf_2 < 0.1:
                print(n,p)

                fig, ax = plt.subplots(1, 1)
                x = np.arange(nbinom.ppf(0.01, n, p),
                              nbinom.ppf(0.99, n, p))
                ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
                ax.vlines(x, 0, nbinom.pmf(x, n, p), colors='b', lw=5, alpha=0.5)
                rv = nbinom(n, p)
                ax.vlines(x, 0, rv.pmf(x), colors='k', linestyles='-', lw=1,
                        label='frozen pmf')
                ax.legend(loc='best', frameon=False)
                plt.show()



# ax.plot(x, nbinom.cdf(x, n, p), 'bo', ms=8, label='nbinom cdf')
