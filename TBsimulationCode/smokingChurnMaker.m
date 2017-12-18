
function [oneMonQuittingUrbanProb, oneMonQuittingRuralProb, oneMonStartingUrbanProb, oneMonStartingRuralProb] = smokingChurnMaker
%Smoking Churn maker
%function call: [oneMonQuittingUrbanProb, oneMonQuittingRuralProb,
%oneMonStartingUrbanProb, oneMonStartingRuralProb] = smokingChurnMaker

%this does not actually go up by fives, but its named that way
fiveYrageBrac = [...
    5
    10
    15
    18
    20
    23
    25
    30
    35
    45
    55
    65
    75
    85
    100
    ];

%whole matrix

%note that urban and rural are the same now.
fiveYrQuittingUrban=[...
    0.835653	0.839122
    0.896673	0.897863
    0.772399	0.855893
    0.423173	0.776012
    0.352349	0.788724
    0.285179	0.728726
    0.144885	0.678519
    0.061432	0.569138
    0.022686	0.341696
    0.026297	0.114469
    0.067979	0.052994
    0.055044	0.054102
    0.023821	0.087116
    0.085997	0.345978
    0           0
    ];

fiveYrQuittingRural=[...
    0.835653	0.839122
    0.896673	0.897863
    0.772399	0.855893
    0.423173	0.776012
    0.352349	0.788724
    0.285179	0.728726
    0.144885	0.678519
    0.061432	0.569138
    0.022686	0.341696
    0.026297	0.114469
    0.067979	0.052994
    0.055044	0.054102
    0.023821	0.087116
    0.085997	0.345978
    0           0
    ];

fiveYrStartingUrban=[...
    0.007581	0.007097667
    0.0209565	0.0150575
    0.3646475	0.0897115
    0.587424833	0.236930167
    0.500550667	0.249215
    0.445547333	0.243918167
    0.394225833	0.2039635
    0.335455	0.135234833
    0.296201833	0.077028167
    0.255669167	0.0674715
    0.131254	0.081969
    0.140175667	0.069094333
    0.218427667	0.061037333
    0.076402667	0.022063167
    0           0
    ];

fiveYrStartingRural=[...
    0.007581	0.007097667
    0.0209565	0.0150575
    0.3646475	0.0897115
    0.587424833	0.236930167
    0.500550667	0.249215
    0.445547333	0.243918167
    0.394225833	0.2039635
    0.335455	0.135234833
    0.296201833	0.077028167
    0.255669167	0.0674715
    0.131254	0.081969
    0.140175667	0.069094333
    0.218427667	0.061037333
    0.076402667	0.022063167
    0           0
    ];


%make one year linear interpolations
oneYrAgeBrac = [5:1:100]';
oneMonQuittingUrbanProb = interp1(fiveYrageBrac,fiveYrQuittingUrban,oneYrAgeBrac);
oneMonQuittingRuralProb = interp1(fiveYrageBrac,fiveYrQuittingRural,oneYrAgeBrac);
oneMonStartingUrbanProb = interp1(fiveYrageBrac,fiveYrStartingUrban,oneYrAgeBrac);
oneMonStartingRuralProb = interp1(fiveYrageBrac,fiveYrStartingRural,oneYrAgeBrac);

%plot(fiveYrageBrac,FiveYrQuitting,'o',oneYrAgeBrac,ans,'-');



