import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as pyplot
fig = pyplot.figure()
fig.set_size_inches((6, 3.5))
xdata = np.array("This data is available in source data for the figure", dtype=float)
ydata = np.array("This data is available in source data for the figure")
lowess = sm.nonparametric.lowess(ydata, xdata, frac=0.4)
pyplot.plot(xdata, ydata, 'o', label='data')
pyplot.plot(lowess[:, 0], lowess[:, 1], label='lowess fit')
pyplot.xlim((-1, 9))
pyplot.ylim((ydata.min() - ydata.max() * 0.1, ydata.max() + ydata.max() * 0.1))
pyplot.xlabel('Sampling time (in Day)')
pyplot.ylabel('Expresion abundance (in TPM)')
pyplot.xticks(np.arange(9), ['2', '4', '6', '8', '10', '12', '14', '16', '18'])
pyplot.legend(loc='best')
pyplot.title('Gene Expression')
fig.savefig("gene.pdf")
