
For better comparison with the MBN and ARG schemes, the accommodation coefficient was allowed to be a free parameter in the single-mode chaos expansions. How much of a problem arises when fixing that parameter to $\alpha \approx 1$?

```python
>>> import numpy as np
>>> import pandas as pd
...
>>> %matplotlib inline
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> sns.set(style='ticks', context='poster')
```

A special sampling of the activation parameter space with $n = 10000$ was created, with all the values for $\alpha$ set to $1$. The design and results are contained in the CSV's here.

```python
>>> design = pd.read_csv("SM_OLS_LHS_fixed_accom_design.csv",
...                      index_col=0)
>>> results = pd.read_csv("SM_OLS_LHS_fixed_accom_design_results.csv",
...                       index_col=0)
```

Some data munging - 

1) Replace "expansion_order" in column names with "OLS"

```python
>>> results.rename(columns=lambda x: x.replace("expansion_order", "OLS"),
...                inplace=True)
```

2) Map the columns "Smax\*" and "Nderiv" to their own DataFrames

```python
>>> Smax = pd.DataFrame({col.replace("Smax_", ""): results[col] \
...                      for col in results.columns if 'Smax' in col})
>>> Smax = np.power(10., Smax)
>>> Smax_ids = Smax.columns
...
>>> Nderiv = pd.DataFrame({col.replace("Nderiv_", ""): results[col] \
...                        for col in results.columns if 'Nderiv' in col})
>>> Nderiv['parcel'] = results['Neq_parcel']
>>> Nderiv = np.power(10., Nderiv)
>>> Nderiv_ids = Nderiv.columns
...
>>> params = [c for c in Nderiv.columns if c != 'parcel']
```

3) Convert `log` variables to normal

```python
>>> for col in design.columns:
...     if col.startswith("log"):
...         design[col] = np.power(10., design[col])
...         design.rename(columns={col: col.replace("log","")},
...                       inplace=True)
```

4) Merge `design` into the results DataFrames

```python
>>> Smax = pd.concat([Smax, design], axis=1)
>>> Nderiv = pd.concat([Nderiv, design], axis=1)
```

5) On rare cases, `Neq` or `Smax` from the parcel model simulations could be so small they were saved as `np.inf`.

```python
>>> Smax.drop(Smax[Smax.parcel == 0.].index, inplace=True)
>>> Nderiv.drop(Nderiv[Nderiv.parcel == 0.].index, inplace=True)
```

## Distribution of relative errors

```python
>>> @np.vectorize
... def rel_err(obs,exp):
...     return 100.*(obs - exp)/exp
...
>>> smax_rel_err = Smax[params].apply(rel_err, args=(Smax.parcel, ))
>>> nderiv_rel_err = Nderiv[params].apply(rel_err, args=(Nderiv.parcel, ))
...
>>> # There are some weird things that happen in the MBN code which produce
... # grossly wrong N, so we filter there
... f = lambda x: x > 1000.
>>> for col in ['MBN', 'OLS_2']:
...     inds = f(nderiv_rel_err[col])
...     nderiv_rel_err[col].ix[inds] = np.NaN
...
...
>>> # "melt" into variable: value format
... agg_smax_re= pd.melt(smax_rel_err,
...                      var_name="param", value_name="relative error")
>>> agg_nderiv_re = pd.melt(nderiv_rel_err,
...                         var_name="param", value_name="relative error")
...
>>> # Drop missing/na values
... agg_smax_re.dropna(inplace=True)
>>> agg_nderiv_re.dropna(inplace=True)
...
>>> print agg_nderiv_re.groupby('param').count()
       relative error
param                
ARG              9999
MBN              9983
OLS_2            9995
OLS_3            9999
OLS_4            9999
OLS_5            9999
```

```python
>>> fig, axes = plt.subplots(2, 1, sharex=True, squeeze=True,
...                         figsize=(10, 9))
>>> ax_smax, ax_nderiv = axes
...
>>> vp_kwargs = dict(color='grey', inner='quartile', lw=0.5)
...
>>> sns.violinplot(x="param", y="relative error", data=agg_smax_re,
...                ax=ax_smax, **vp_kwargs)
>>> sns.violinplot(x="param", y="relative error", data=agg_nderiv_re,
...                ax=ax_nderiv, **vp_kwargs)
...
>>> for ax in axes:
...     ax.hlines(0, -1, 1+len(params), 'grey')
...     ax.set_ylim(-150, 150)
...
>>> sns.despine(fig=fig)
```

```python

```
