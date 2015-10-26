```python
>>> import parcel_model as pm
...
>>> import numpy as np
>>> import pandas as pd
...
>>> # Use original paper archive
... # from plot_data_setup import get_data, PAPER_DIR
... # data = get_data('both', 'results')
...
... # Use local copy
... design = pd.read_csv("design.csv")
>>> results = pd.read_csv("results.csv")
>>> data = pd.concat([design, results], axis=1)
...
>>> %matplotlib inline
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> sns.set(style="ticks")
```

```python
>>> # Extract Smax and cloud droplet number concentrations
... Smax = np.power(10., data.Smax_parcel)
>>> Neq = np.power(10., data.Neq_parcel)
>>> Nkn = np.power(10., data.Nkn_parcel)
>>> results = pd.DataFrame({'Smax': Smax, 'Neq': Neq, 'Nkn': Nkn})
...
>>> # Re-name and project log design variables
... design_cols = ['logN', 'logmu', 'sigma', 'kappa',
...                'T', 'P', 'logV', 'accom']
>>> design = data[design_cols]
>>> for col in design_cols:
...     if col.startswith("log"):
...         design[col] = np.power(10., design[col])
...         design.rename(columns={col:col[3:]}, inplace=True)
...
>>> results['act_frac_eq'] = results.Neq / design.N
>>> results['act_frac_kn'] = results.Nkn / design.N
/Users/daniel/anaconda/lib/python2.7/site-packages/ipykernel/__main__.py:13: SettingWithCopyWarning: 
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead

See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
```

Compute the "derived" cloud droplet number - the number concentration predicted using Kohler theory and $S_\text{max}$.

```python
>>> activate = np.vectorize(pm.activation.lognormal_activation)
...
>>> Nderiv, act_frac_deriv = activate(results['Smax'], design['mu']*1e-6,
...                                   design.sigma, design.N, design.kappa,
...                                   T=design['T'], approx=True)
>>> results['Nderiv'] = Nderiv
>>> results['act_frac_deriv'] = act_frac_deriv
```

Correlations between equilibrium number concentration and inertially-limited estimates.

```python
>>> results[results['act_frac_eq'] < 1e-5]
               Neq           Nkn      Smax   act_frac_eq   act_frac_kn  \
1229  0.000000e+00  0.000000e+00  0.000003  0.000000e+00  0.000000e+00   
5670  1.841859e-07  1.841859e-07  0.007699  1.547593e-09  1.547593e-09   

            Nderiv  act_frac_deriv  
1229  2.051739e+00    4.153688e-04  
5670  2.230879e-07    1.874461e-09
```

```python
>>> fig, axes = plt.subplots(1, 2, sharex=True, sharey=True,
...                          figsize=(16,5))
...
>>> ax_deriv, ax_eq = axes
...
>>> # Nderiv
... sns.regplot('Nderiv', 'Nkn', results, ax=ax_deriv,
...             robust=True, ci=None)
>>> # Neq
... sns.regplot('Neq', 'Nkn', results, ax=ax_eq,
...             robust=True, ci=None)
...
>>> for ax in axes:
...     nn = np.linspace(0, 1e4, 100)
...     ax.plot(nn, nn, 'k-', zorder=-999)
...     ax.set_xlim(0, 1e4)
...     ax.set_ylim(0, 1e4)
```

```python
>>> too_big = results['Nkn'] > results['Neq']
...
>>> n_rows = 3
>>> n_cols = int(np.ceil(len(design.columns) / (n_rows*1.)))
>>> aspect, size = 16./10., 3.
>>> figsize = (n_cols * size * aspect, n_rows * size)
>>> fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
...
>>> for col, ax in zip(design.columns, list(axes.flat)):
...     all_data = design[~too_big][col]
...     bad_data = design[too_big][col]
...
...     summary_str = (
...         "{:6s}: {:2.2e} < bad < {:2.2e} | {:2.2e} < all < {:2.2e}"
...             .format(col, bad_data.min(), bad_data.max(),
...                     all_data.min(), all_data.max())
...     )
...     print summary_str
...
...     # Merge into'tidy' dataset
...     _df = pd.DataFrame({'all': all_data,
...                         'Nkn > Neq': bad_data})
...     _df = pd.melt(_df, var_name="set", value_name=col)
...     _df.dropna(inplace=True)
...
...     sns.distplot(bad_data, label='Nkn > Neq', color='r', ax=ax,
...                  norm_hist=False)
...     sns.distplot(all_data, label='all', color='k', ax=ax,
...                  norm_hist=False)
...
>>> ax = axes[0, 1]
>>> ax.legend(loc='upper center')
...
>>> sns.despine(fig)
>>> plt.tight_layout()
N     : 1.34e+02 < bad < 9.97e+03 | 1.06e+01 < all < 1.00e+04
mu    : 1.00e-03 < bad < 2.83e-03 | 1.34e-03 < all < 1.00e-01
sigma : 1.23e+00 < bad < 2.99e+00 | 1.20e+00 < all < 3.00e+00
kappa : 2.88e-03 < bad < 1.20e+00 | 8.66e-05 < all < 1.20e+00
T     : 2.40e+02 < bad < 3.10e+02 | 2.40e+02 < all < 3.10e+02
P     : 5.01e+04 < bad < 1.04e+05 | 5.00e+04 < all < 1.05e+05
V     : 1.59e-01 < bad < 9.89e+00 | 1.03e-02 < all < 1.00e+01
accom : 1.00e-01 < bad < 9.89e-01 | 1.00e-01 < all < 1.00e+00
```

```python
>>> len(_df.dropna())
10094
```

---
!!python/unicode 'scrolled': false
...

```python
>>> d = design.ix[421]
>>> print d
>>> aers = [pm.AerosolSpecies('test', pm.Lognorm(d.mu, d.sigma, d.N), d.kappa,
...                           bins=250), ]
>>> model = pm.ParcelModel(aers, V=d.V, T0=d['T'], S0=-0.0, P0=d.P,
...                        accom=d.accom,
...                        truncate_aerosols=True, console=False)
>>> par, aer = model.run(t_end=1800, output_dt=0.1, solver_dt=1.0,
...                      solver='cvode', output='dataframes',
...                      terminate=True, terminate_depth=10.)
...
>>> Smax = par['S'].max()
>>> T_fin = par['T'].iloc[0]
>>> pm.binned_activation(Smax, T_fin, aer['test'].iloc[-1].values, aers[0])
N         9857.694798
mu           0.001254
sigma        1.829266
kappa        0.343324
T          280.141008
P        65527.767602
V            6.015705
accom        0.319689
Name: 421, dtype: float64
(0.0024756855670004709,
 0.99998468905020799,
 403.92233261745946,
 0.79682097243948813)
```

```python
>>> kc = np.vectorize(pm.thermo.kohler_crit)
...
>>> irads = np.ones((len(aer['test'], )))*-9999
>>> N = 0.
>>> for i in xrange(len(aer['test'].columns)):
...     #i = -(i+1)
...
...     r = aer['test'].columns[i]
...     print i, r,
...
...     ts_radii = aer['test'][r]
...     r_dry = aers[0].r_drys[i]
...     kappa = aers[0].kappa
...     T = par['T'].iloc[0]
...
...     r_crit, s_crit = kc(T, r_dry, kappa, approx=False)
...
...     for irad, ts_rad in enumerate(ts_radii):
...         if ts_rad > r_crit:
...             print irad, ts_rad, r_crit
...             irads[i] = irad
...             break
...
...     if irads[i] >= 0:
...         N += aers[0].Nis[i]
...
...     print
0 r000
1 r001
2 r002
3 r003
4 r004
5 r005
6 r006
7 r007
8 r008
9 r009
10 r010
11 r011
12 r012
13 r013
14 r014
15 r015
16 r016
17 r017 0 1.46308181508e-10 1.37770173941e-10

18 r018 0 1.49525357373e-10 1.41795287243e-10

19 r019 0 1.52813462191e-10 1.45915089059e-10

20 r020 0 1.56174088204e-10 1.50131806885e-10

21 r021 0 1.59608866847e-10 1.54447720615e-10

22 r022 0 1.63119469997e-10 1.58865163775e-10

23 r023 0 1.66707611264e-10 1.63386524789e-10

24 r024 0 1.70375047325e-10 1.68014248266e-10

25 r025 0 1.7412357931e-10 1.72750836322e-10

26 r026 0 1.77955054227e-10 1.77598849935e-10

27 r027
28 r028
29 r029
30 r030
31 r031
32 r032
33 r033
34 r034
35 r035
36 r036
37 r037 0 2.26154323839e-10 1.97669068261e-10

38 r038 0 2.31144023157e-10 2.03103269197e-10

39 r039 0 2.36245581785e-10 2.08665306637e-10

40 r040 0 2.41461663739e-10 2.14358187866e-10

41 r041 0 2.46795005642e-10 2.20184990907e-10

42 r042 0 2.52248418874e-10 2.26148866197e-10

43 r043 0 2.57824791763e-10 2.32253038282e-10

44 r044 0 2.63527091838e-10 2.38500807565e-10

45 r045 0 2.69358368129e-10 2.44895552089e-10

46 r046 0 2.75321753527e-10 2.51440729364e-10

47 r047 0 2.81420467205e-10 2.58139878236e-10

48 r048 0 2.87657817098e-10 2.64996620798e-10

49 r049 0 2.94037202449e-10 2.71520331423e-10

50 r050 0 3.00562116428e-10 2.77907693757e-10

51 r051 0 3.07236148814e-10 2.8444531518e-10

52 r052 0 3.14062988764e-10 2.91136730451e-10

53 r053 0 3.2104642765e-10 2.97985557485e-10

54 r054 0 3.28190361993e-10 3.04995499304e-10

55 r055 0 3.35498796473e-10 3.12170346042e-10

56 r056 0 3.4297584704e-10 3.19513976995e-10

57 r057 0 3.50625744118e-10 3.27030362715e-10

58 r058 0 3.58452835913e-10 3.01390228863e-10

59 r059 0 3.66461591829e-10 3.09264411478e-10

60 r060 0 3.74656605996e-10 3.17323829789e-10

61 r061 0 3.83042600914e-10 3.25572841361e-10

62 r062 0 3.9162443122e-10 3.34015906268e-10

63 r063 0 4.0040708759e-10 3.42657589504e-10

64 r064 0 4.09395700768e-10 3.51502563454e-10

65 r065 0 4.18595545738e-10 3.60555610418e-10

66 r066 0 4.28012046047e-10 3.69821625194e-10

67 r067 0 4.37650778273e-10 3.79305617731e-10

68 r068 0 4.4751747666e-10 3.89012715832e-10

69 r069 0 4.57618037913e-10 3.98948167929e-10

70 r070 0 4.6795852617e-10 4.08382659607e-10

71 r071 0 4.78545178152e-10 4.17989630857e-10

72 r072 0 4.89384408496e-10 4.27822600676e-10

73 r073 0 5.00482815292e-10 4.37886885552e-10

74 r074 0 5.1184718581e-10 4.48187927042e-10

75 r075 0 5.23484502442e-10 4.5873129471e-10

76 r076 0 5.35401948866e-10 4.69522689144e-10

77 r077 0 5.47606916426e-10 4.80567945033e-10

78 r078 0 5.60107010751e-10 5.3290597985e-10

79 r079 0 5.72910058623e-10 5.45442291829e-10

80 r080 0 5.86024115084e-10 5.58273513461e-10

81 r081 0 5.99457470818e-10 5.7140658233e-10

82 r082 0 6.13218659793e-10 5.84848599221e-10

83 r083 0 6.27316467194e-10 5.98606831962e-10

84 r084 0 6.4175993764e-10 6.12688719353e-10

85 r085 0 6.56558383708e-10 6.27101875185e-10

86 r086 0 6.71721394769e-10 6.41854092362e-10

87 r087 0 6.8725884615e-10 6.5695334711e-10

88 r088 0 7.03180908629e-10 6.72407803292e-10

89 r089 0 7.19498058286e-10 6.8822581682e-10

90 r090 0 7.3622108671e-10 7.04415940177e-10

91 r091 0 7.53361111589e-10 7.20986927035e-10

92 r092 0 7.70929587677e-10 7.37947736993e-10

93 r093 0 7.88938318175e-10 7.55307540419e-10

94 r094 0 8.07399466528e-10 7.73075723409e-10

95 r095 0 8.26325568651e-10 7.91261892861e-10

96 r096 0 8.45729545607e-10 8.0987588167e-10

97 r097 0 8.65624716752e-10 8.62261099675e-10

98 r098 0 8.86024813355e-10 8.83289676967e-10

99 r099 3 9.06983383391e-10 9.06973494449e-10

100 r100
101 r101 0 9.50398447629e-10 9.43057114508e-10

102 r102 0 9.72964302467e-10 9.64457852286e-10

103 r103 0 9.96110430741e-10 9.92562883761e-10

104 r104
105 r105
106 r106
107 r107
108 r108
109 r109
110 r110
111 r111
112 r112
113 r113
114 r114
115 r115
116 r116
117 r117
118 r118
119 r119
120 r120
121 r121
122 r122
123 r123
124 r124
125 r125
126 r126
127 r127
128 r128
129 r129
130 r130
131 r131
132 r132
133 r133
134 r134
135 r135
136 r136
137 r137
138 r138
139 r139
140 r140
141 r141
142 r142
143 r143
144 r144
145 r145
146 r146
147 r147
148 r148
149 r149
150 r150
151 r151
152 r152
153 r153
154 r154
155 r155
156 r156
157 r157
158 r158
159 r159
160 r160
161 r161
162 r162
163 r163
164 r164
165 r165
166 r166
167 r167
168 r168
169 r169
170 r170
171 r171
172 r172
173 r173
174 r174
175 r175
176 r176
177 r177
178 r178
179 r179
180 r180
181 r181
182 r182
183 r183
184 r184
185 r185
186 r186
187 r187
188 r188
189 r189
190 r190
191 r191
192 r192
193 r193
194 r194 173 1.58296971106e-08 1.57781977777e-08

195 r195 161 1.63441999082e-08 1.63162753673e-08

196 r196 153 1.69633114218e-08 1.68693553826e-08

197 r197 146 1.75497934376e-08 1.74460337912e-08

198 r198 140 1.81757046346e-08 1.80371021384e-08

199 r199 134 1.86752858913e-08 1.86404988403e-08

200 r200 129 1.93089883937e-08 1.92753047395e-08
```

Something weird is going on with the kinetic limitation calculation for this particular example, so let's dig into it. First, plot the original size distribution and the final wet size distribution.

```python
>>> a = aers[0]
>>> r_dry = a.r_drys*1e6
>>> r_wet = aer['test'].iloc[-1]*1e6
...
>>> plt.plot(r_dry)
>>> plt.plot(r_wet)
[<matplotlib.lines.Line2D at 0x10d3d43d0>]
```

Now, diagnose the equilibrium droplet number concentration and the kinetic droplet number concentration. Droplets counted under the "equilibrium" description are simply those whose supersaturation is less than or equal the maximum supersaturation achieved in the parcel. Counted under the "kinetic" description are all the droplets *larger* than the smallest particle that has strictly activated (grown larger than its critical value).

```python
>>> kohler_crit = pm.activation.kohler_crit
>>> T = par['T'].iloc[-1]
>>> Smax = par['S'].max()
>>> rs = aer['test'].iloc[-1]
...
>>> a = aers[0]
>>> kappa = a.kappa
>>> r_drys = a.r_drys
>>> Nis = a.Nis
>>> N_tot = np.sum(Nis)
...
>>> r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa, False) \
...                        for r_dry in r_drys])
>>> s_crits = np.array(s_crits)
>>> r_crits = np.array(r_crits)
...
>>> ## Equilibrium calculation
... activated_eq = (Smax >= s_crits)
>>> N_eq = np.sum(Nis[activated_eq])
>>> eq_frac = N_eq/N_tot
...
>>> ## Kinetic calculation
... is_r_large = rs >= r_crits
>>> if not np.any(is_r_large):
...     N_kn, kn_frac = 0., 0.
...     phi = 1.0
>>> else:
...     smallest_ind = np.where(is_r_large)[0][0]
...     N_kn = np.sum(Nis[smallest_ind:])
...     kn_frac = N_kn/N_tot
...
...     ## Unactivated - all droplets smaller than their critical size
...     droplets = range(smallest_ind, len(Nis))
...     Nis_drops = Nis[droplets]
...     r_crits_drops = r_crits[droplets]
...     rs_drops = rs[droplets]
...     too_small = (rs_drops < r_crits_drops).values
...
...     N_unact = np.sum(Nis_drops[too_small])
...
...     phi = N_unact/N_kn
...
>>> alpha = N_kn/N_eq
```

```python
>>> hasattr(rs.values, 'values')
```

```python

```
