function[forceOfInfection] = forceOfInfectionCalculator(susecptibleWithLat, riskSensTB, riskMdrTB, totInfectious)
%finding force of infection OLD NOT USING ANYMORE
% %count how many people in each age grp at last time period
% 
%     alive = sum(totPplInAgeBrac);    
%     fracPplInAge = totPplInAgeBrac ./repmat(alive, 1, 15);
%     contactsPerDay = (fracPplInAge)*[...
%         5.45
%         7.64
%         6.66
%         6.44
%         5.71
%         5.74
%         6.02
%         6.35
%         5.66
%         4.42
%         3.31
%         3.29
%         2.61
%         1.95
%         2.84
%         ];
% 
% forceOfInfection = contactsPerDay*365*0.0029;


SensTransmitted_calc = rand(length(susecptibleWithLat),1) <= riskSensTB(susecptibleWithLat);
MDRtransmitted_calc  = rand(length(susecptibleWithLat),1) <= riskMdrTB(susecptibleWithLat);
forceOfInfection = sum(max(SensTransmitted_calc, MDRtransmitted_calc)) / totInfectious;
