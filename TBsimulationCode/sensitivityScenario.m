disp('################ SENSITIVITY PARAMETERS #############')
format longG
if sensitivityScenarioNum == 1
    %do nothing
elseif (sensitivityScenarioNum >= 2 && sensitivityScenarioNum <= 2+25)
    %geneX and DST cost
    %     TBparams.geneXpertCost = 18.67629239;  %base 18.67629239, range 14.9410 to 23 (volume discount to about the sa
    %     TBparams.DSTcost = 26.57162211;  %base 26.57162211, range go from 20 to 33 (about 75% to 125%).  This is a similar percentage range given by Menzies for DST costs in other countries.
    
    geneXpertCostArray = linspace(14.94103391, 23, 5);
    DSTcostArray = linspace(20, 33, 5);
    
    counter = 2;
    for i = 1:5
        for j = 1:5
            costs(1:2,counter) = [geneXpertCostArray(1,i);DSTcostArray(1,j)];
            counter = counter + 1;
        end
    end
    TBparams.geneXpertCost = costs(1, sensitivityScenarioNum);
    TBparams.DSTcost = costs(2, sensitivityScenarioNum);
    fprintf('GeneXcost is: %d \n', costs(1, sensitivityScenarioNum))
    fprintf('DSTcost is: %d \n', costs(2, sensitivityScenarioNum))
    
elseif (sensitivityScenarioNum == 28 )
    
    % SSsens+ 0.90 (0.89 to 0.91)
    % SSspec+ 0.98 (0.98 to 0.99)
    % DSTsens+ 0.94 (0.92 to 0.96)
    % DSTspec+ 0.97 (0.96 to 0.98)
    
    %worst GeneX
    TBparams.geneXsensDS=0.89;
    TBparams.geneXspecDS=0.98;
    TBparams.geneXsensMDR=0.92;
    TBparams.geneXspecMDR=0.96;
    disp('geneXsensDS, geneXspecDS, geneXsensMDR, geneXspecMDR')
    [0.89,0.98,0.92,0.96]
elseif (sensitivityScenarioNum == 29)
    %best GeneX
    TBparams.geneXsensDS=0.91;
    TBparams.geneXspecDS=0.99;
    TBparams.geneXsensMDR=0.96;
    TBparams.geneXspecMDR= 0.98;
    disp('geneXsensDS, geneXspecDS, geneXsensMDR, geneXspecMDR')
    [0.91, 0.99, 0.96, 0.98]
elseif (sensitivityScenarioNum >= 30 && sensitivityScenarioNum <= 53)
    
    geneXpertCostArray = linspace(0, 1, 5);
    DSTcostArray = linspace(0, 1, 5);
    
    counter = 1;
    for i = 1:5
        for j = 1:5
            valMat(1:2,counter) = [geneXpertCostArray(1,i);DSTcostArray(1,j)];
            counter = counter + 1;
        end
    end
    %pertinent counter values are 2 to 25 inclusive
    TBparams.privateReferToRNTCP = valMat(1, sensitivityScenarioNum-28);
    TBparams.PPMeffectiveness = valMat(2, sensitivityScenarioNum-28);
    fprintf('privateReferToRNTCP is: %d \n', valMat(1, sensitivityScenarioNum-28))
    fprintf('PPMeffectiveness is: %d \n', valMat(2, sensitivityScenarioNum-28))
    
elseif (sensitivityScenarioNum == 54)
    %increase diagnosis to 150% of what it was before
    TBMAC1_increaseTrt = 1
    TBmac1ScalingFac = repmat(1.5,1,120);
    fprintf('diagnosis factor is increased by 1.5 starting in 2015')
    
elseif (sensitivityScenarioNum == 55)  %rationale: additional burden on RNTCP reduces quality
    TBparams.CatIdeath = TBparams.CatIdeath*2
    TBparams.CatIdefault = TBparams.CatIdefault*2
    TBparams.CatIIdeath = TBparams.CatIIdeath*2
    TBparams.CatIIdefault = TBparams.CatIIdefault*2
elseif (sensitivityScenarioNum == 56)
    TBparams.DSTperiodLength = 2
    
elseif (sensitivityScenarioNum >= 57 && sensitivityScenarioNum <= 59)
    multiplier = [1.2,1.5,2];
    
    TBparams.perPatientPPMcost = 34.72322628*multiplier(sensitivityScenarioNum-56);
    fprintf('perPatientPPMcost x %d \n',multiplier(sensitivityScenarioNum-56))
    
elseif (sensitivityScenarioNum == 60)  %just seeing if this has any effect at all...
    TBparams.CatIdeath = TBparams.CatIdeath*20
    TBparams.CatIdefault = TBparams.CatIdefault*20
    TBparams.CatIIdeath = TBparams.CatIIdeath*20
    TBparams.CatIIdefault = TBparams.CatIIdefault*20
    
    
end





disp('################ END SENSITIVITY PARAMETERS #############')

