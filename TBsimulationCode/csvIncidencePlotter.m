function csvIncidencePlotter(fullPlottingMatrix, folderName, plotResolution)

%function csvDeathsPlotter(plottingMatrix, plotResolution)

%

%plottingMatrix is a (number of time periods) x 17 matrix for plotting,

%made from TBsimulation with:

%  col 1: nonsmoking natural deaths

%  col 2: smoking natural deaths

%  col 3: nonsmoking untreatedTB deaths

%  col 4: smoking untreatedTB deaths

%  col 5: catI nonsmoking deaths 

%  col 6: catI smoking deaths

%  col 7: catII nonsmoking deaths

%  col 8: catII smoking deaths

%  col 9: catIV nonsmoking deaths

%  col 10: catIV smoking deaths

%  col 11: nonsmoking age-related deaths (reach age 99)

%  col 12: smoking age-related deaths (reach age 99)

%  col13 - 17: all active mdrTB related deaths (see code in TBsimulation_july6 if you

%  need details)

%  col 18: all active sensTB related deaths

%

%folderName is where you want the graphs to print

%plotResolution should be '-r70'





%--------------%%GRAB CHARACTERISTIC COLUMNS%%-------------------%



%not need since indidence matrix isn't subset



%--------------%%MAKE GRAPHS%%-------------------%



   numPeriods = size(fullPlottingMatrix,1);

    

    yrlyPlottingMat = zeros(floor(numPeriods/12),size(fullPlottingMatrix,2));

    for t = 1:floor(numPeriods/12)

        yrlyPlottingMat(t,:) = sum(fullPlottingMatrix(( ( (t-1)*12 )+1):(12*t),:));

        yrlyPlottingMat(t,end) = fullPlottingMatrix( ( (t-1)*12 )+7 ,end) ;  %the last column is the total num of alive, which is not summed

    end

    

    %%%%%%%%%%%%%incidence graph%%%%%%%%%%

    %PLOT

    toPlot = yrlyPlottingMat(120:end,1:2) ./ [yrlyPlottingMat(120:end,3) ,yrlyPlottingMat(120:end,3)  ] ;
    x = [1986:1:2046];
    x = x(1:size(toPlot,1));

    bar(x,toPlot);

    title('TB Incidence');

    legend('sensTBincidence', 'mdrTBincidence');

    xlabel('Year');

    ylabel('Incidence');

    print('-dpng',plotResolution,[folderName '/' 'TBincidence']);

    tablePrinter('sensTBincidence, mdrTBincidence, numAlive', yrlyPlottingMat, 'TBincidence', folderName);

    close all



	%do incidence by rate (ie, total active out of 100,000) and put

	%the WHO numbers on too for comparison

%    dataTime = [124	; 125	; 126	; 127	; 128	; 129	; 130	; 131	; 132	; 133	; 134	; 135	; 136	; 137	; 138	; 139	; 140	; 141	; 142	; 143	; 144];

    dataTime = [1990	; 1991	; 1992	; 1993	; 1994	;  1995	; 1996	; 1997	; 1998	; 1999	; 2000	; 2001	; 2002	; 2003	; 2004	; 2005	; 2006	; 2007	; 2008	; 2009	; 2010	];

    actTBinci = [216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 215	; 214	; 212	; 209	; 205	; 201	; 196	; 190	; 185 ];

	distToUpperLimit = [39	; 37	; 35	; 33	; 32	; 30	; 29	; 27	; 26	; 25	; 24	; 23	; 23	; 23	; 22	; 22	; 22	; 21	; 21	; 21	; 20	];

	distToLowerLimit = [35	;  33	; 32	; 30	; 29	; 27	; 26	; 25	; 24	; 23	; 22	; 22	; 22	; 22	; 22	;21	    ; 21	; 20	; 20	; 19	; 18	];   



	%make simulation incidences into rates out of 100000 and plot

    smallx = [1990:1:2011];

	toPlot = toPlot * 100000;

	%p = plot(x,toPlot);

    p = plot(smallx,toPlot(5:26));

	set(p,'Color','red','LineWidth',1);

	title('Total TB incidence rate, out of 100000, with WHO estimates');

    xlabel('Year');

    ylabel('Incidence');

	hold on;

	errorbar(dataTime, actTBinci,distToLowerLimit,distToUpperLimit);

    ylim([0 600]);

	print('-dpng',plotResolution,[folderName '/' 'TBincidence_WHOcomparison.png']);

	close all