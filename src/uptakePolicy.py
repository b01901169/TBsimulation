import numpy as np

uptakeAgeBracs = [
    0,	4,
    5,	9,
    10,	14,
    15,	19,
    20,	24,
    25,	29,
    30,	34,
    35,	39,
    40,	44,
    45,	49,
    50,	54,
    55,	59,
    60,	64,
    65,	69,
    70,	74,
    75,	79,
    80,	84,
    85,	109,
    ];
uptakeAgeBracs = np.resize(uptakeAgeBracs, (18, 2)).tolist()

uptakeUrbanKnowledge = [
    0.0282572,       0.0265703,
    0.03421145,      0.034611175,
    0.0391262,       0.0411383,
    0.04314695,      0.046267925,
    0.0464192,       0.0501163,
    0.04908845,      0.052799675,
    0.13294343,      0.12393383,
    0.1225361,       0.12221708,
    0.08431062,      0.08530207,
    0.08464646,      0.08522047,
    0.07465641,      0.08945654,
    0.09882135,      0.11919841,
    0.0631112,       0.0487283,
    0.06619445,      0.046276675,
    0.0699842,       0.0437063,
    0.07462595,      0.041133425,
    0.0802652,       0.0386743,
    0.08704745,      0.036445175,
    ];
uptakeUrbanKnowledge = np.resize(uptakeUrbanKnowledge, (18, 2)).tolist()

finerUptakeAgeBracs = []
for i in range(86): # TODO
    finerUptakeAgeBracs.append([i, i])
finerUptakeAgeBracs[85][1] = 109

finerUptakeUrbanKnowledge = np.zeros((110, 2))
for i in range(110):
    index = int(i / 5)
    if i > 85:
        finerUptakeUrbanKnowledge[i] = uptakeUrbanKnowledge[-1]
    else:
        finerUptakeUrbanKnowledge[i] = uptakeUrbanKnowledge[index]


finerUptakeUrbanKnowledge = finerUptakeUrbanKnowledge.tolist()

FOI_cal = 0.0045
modContactMat = FOI_cal * 1 * np.array([[0.69, 0.24, 0.36, 0.76, 0.79, 0.53, 0.97, 0.83, 0.38, 0.22, 0.36, 0.22, 0.22, 0.05, 0.09],
[0.59, 1.82, 0.43, 0.22, 0.79, 0.77, 0.42, 1.15, 1.09, 0.43, 0.21, 0.32, 0.15, 0.08, 0.09],
[0.25, 0.53, 1.69, 0.32, 0.09, 1.11, 0.69, 0.51, 0.58, 0.22, 0.14, 0.1, 0.2, 0.07, 0.1],
[0.18, 0.44, 0.79, 1.4, 0.38, 0.2, 0.67, 0.75, 0.47, 0.74, 0.47, 0.1, 0.1, 0.08, 0.08],
[0.42, 0.57, 0.17, 0.85, 1.02, 0.58, 0.31, 0.2, 0.27, 0.46, 0.22, 0.14, 0.1, 0.07, 0.19],
[0.61, 0.42, 0.44, 0.12, 0.49, 0.49, 0.42, 0.24, 0.29, 0.29, 0.17, 0.22, 0.14, 0.12, 0.1],
[0.57, 0.68, 0.32, 0.37, 0.28, 0.35, 0.8, 0.47, 0.25, 0.17, 0.15, 0.13, 0.15, 0.08, 0.05],
[0.74, 0.99, 0.51, 0.29, 0.46, 0.21, 0.21, 0.76, 0.49, 0.17, 0.14, 0.21, 0.2, 0.16, 0.13],
[0.18, 0.66, 0.69, 0.48, 0.23, 0.38, 0.27, 0.47, 0.45, 0.35, 0.19, 0.05, 0.19, 0.11, 0.24],
[0.2, 0.15, 0.27, 0.51, 0.44, 0.25, 0.38, 0.29, 0.4, 0.58, 0.27, 0.25, 0.13, 0.04, 0.29],
[0.36, 0.17, 0.2, 0.27, 0.2, 0.3, 0.23, 0.18, 0.12, 0.15, 0.33, 0.17, 0.09, 0.05, 0.23],
[0.15, 0.19, 0.19, 0.13, 0.28, 0.3, 0.2, 0.11, 0.19, 0.19, 0.35, 0.52, 0.24, 0.07, 0.22],
[0.18, 0.41, 0.17, 0.09, 0.11, 0.2, 0.23, 0.32, 0.23, 0.14, 0.2, 0.32, 0.27, 0.17, 0.11],
[0.26, 0.3, 0.26, 0.26, 0.15, 0.04, 0.15, 0.04, 0.15, 0.11, 0.04, 0.44, 0.3, 0.37, 0.22],
[0.07, 0.07, 0.17, 0.37, 0, 0.03, 0.07, 0.03, 0.3, 0.2, 0.07, 0.1, 0.13, 0.43, 0.7]])

modTransMatAgeBrac = [(i*5, i*5+5) for i in range(14)] + [(70, 110)]

