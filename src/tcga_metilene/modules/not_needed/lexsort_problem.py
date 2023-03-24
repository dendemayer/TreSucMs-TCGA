import pandas as pd
import numpy as np
"""
in case of identical index entries in different MI levels, a sorting is
necessary to avoid the:
<stdin>:1: PerformanceWarning: indexing past lexsort depth may impact
performance.
"""
mux = pd.MultiIndex.from_arrays([
    list('aaaabbbbbccddddd'),
    list('tuvwtuvwtuvwtuvw')
], names=['one', 'two'])

df = pd.DataFrame({'col': np.arange(len(mux))}, mux)
df.loc[('c', 'u'), :]
# this gives <stdin>:1: PerformanceWarning: indexing past lexsort depth may
# impact performance.

# sorting the index before slicing:
df2 = df.sort_index()
df2.loc[('c', 'u'), :]
