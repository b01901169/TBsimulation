import numpy as np

ageBrac = [
[0   ,   0  ],
[1   ,   4  ], 
[5   ,   9  ], 
[10  ,   14 ], 
[15  ,   19 ], 
[20  ,   24 ], 
[25  ,   29 ], 
[30  ,   34 ], 
[35  ,   39 ], 
[40  ,   44 ], 
[45  ,   49 ], 
[50  ,   54 ], 
[55  ,   59 ], 
[60  ,   64 ], 
[65  ,   69 ], 
[70  ,   74 ], 
[75  ,   79 ], 
[80  ,   84 ], 
[85  ,   89 ], 
[90  ,   94 ], 
[95  ,   98 ], 
[99  ,   109],
]

matlab_death_rate = [
    0.0305,
    0.0253,
    0.0249,
    0.0248,
    0.0249,
    0.0249,
    0.0249,
    0.0250,
    0.0250,
    0.0251,
    0.0253,
    0.0256,
    0.0261,
    0.0269,
    0.0282,
    0.0300,
    0.0324,
    0.0353,
    0.0394,
    0.0452,
    0.0537,
    1.0000
]

matlab_mu = [
    0.0059,
    0.0006,
    0.0002,
    0.0001,
    0.0002,
    0.0002,
    0.0003,
    0.0003,
    0.0003,
    0.0004,
    0.0006,
    0.0009,
    0.0015,
    0.0023,
    0.0036,
    0.0055,
    0.0079,
    0.0108,
    0.0150,
    0.0211,
    0.0297,
    1.0000
]

def getRateParameters(yrsOfAnalysis=25):
    death_rate = np.zeros(110)
    tmp_mu = np.zeros(110)
    mu = np.zeros((110,yrsOfAnalysis))
    for i in range(22):
        for age in range(ageBrac[i][0], ageBrac[i][1]+1):
            death_rate[age] = matlab_death_rate[i]
            tmp_mu[age] = matlab_mu[i]
    for y in range(yrsOfAnalysis):
        mu[:, y] = tmp_mu

    return mu, death_rate
