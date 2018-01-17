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
for i in range(86):
    finerUptakeAgeBracs.append([i, i])
finerUptakeAgeBracs[85][1] = 109

finerUptakeUrbanKnowledge = np.zeros((86, 2))
for i in range(86):
    index = int(i / 5)
    finerUptakeUrbanKnowledge[i] = uptakeUrbanKnowledge[index]

finerUptakeUrbanKnowledge = finerUptakeUrbanKnowledge.tolist()

