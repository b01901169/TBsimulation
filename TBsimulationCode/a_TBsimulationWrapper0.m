%TBsimulationWrapper

%open opening folder where 'thingsToRun' lives
%==

%set up stuff to send me an email when it is done
% Define these variables appropriately:
% mail = 'TBsimulation@gmail.com'; %Your GMail email address
% password = 'rntcpinindia'; %Your GMail password
% 
% % Then this code will set up the preferences properly:
% setpref('Internet','E_mail',mail);
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username',mail);
% setpref('Internet','SMTP_Password',password);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');

% Send the email. Note that the first input is the address you are sending the email to
%sendmail('suensze@gmail.com','Test from MATLAB','Hello! This is a test from MATLAB!')

%defaults for graphs.  The usual run: run this with calibrationRun off in the overwriter
%TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 1, 60, 0, 'inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake')
%defaults for stata output or latent-to-active analysis.  Run this with calibrationRun on in the overwriter
%TBsimulation_july6('.', 'This is the calibration simulation run.',230, 4800,  '-r0', 1, 1, 1561, 'inf0p0020_lat2p55_noTreat_fullCatIV_empUptake_230')
%to do life expectancy/mortality stuff
%TBsimulation_july6('.', 'This is the calibration simulation run.',230, 8000,  '-r0', 1, 60, 0, 'inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake')


%=======

%%%%this reinitializes the random number generator used by rand with a seed based on the current time.
rng('shuffle');
    %%%base case
%  %TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 'NA', 2013,0, 2.13, 0.0024, 0.0042, 2.4, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')
%remember if you change function call after burn in you need to put in the special section of overwriter or redo the burn in

%%% %base case with that writes preloaded burnIn seed
%  %TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 'writeBurnIn', 2013,0, 2.13, 0.0024, 0.0042, 2.4, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')


    %%%%%% %base case with that uses preloaded burnIn seed
%  %TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 'loadBurnIn', 2013,0, 2.13, 0.0024, 0.0042, 2.4, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')

    %Counterfactual: no dots
%  %TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 'loadBurnIn', 2013,0, 2.13, 0.0024, 0.0042, 2.4, 'noTreat_fullCatIV_empUpta')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCENARIOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % base case
% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','writeBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','loadBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','loadBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'GeneXallRNTCP_2014')

% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','loadBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'GeneXDST_2014')
% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','loadBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'PPM0p7_0p6_2014')
% TBsimulation_jan23('.', 's01',160, 20000,  '-r70','loadBurnIn', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'PPM0p8_2014')


TBsimulation_jan23('.', 'b01',160, 200,  '-r70','NA', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
% 
% TBsimulation_jan23('.', 'b01',160, 20000,  '-r70','LEbuilder', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
% TBsimulation_jan23('.', 'b01',160, 20000,  '-r70','LEbuilder', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'GeneXallRNTCP_2014')
% TBsimulation_jan23('.', 'b01',160, 20000,  '-r70','LEbuilder', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'GeneXDST_2014')
% TBsimulation_jan23('.', 'b01',160, 20000,  '-r70','LEbuilder', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'PPM0p7_0p6_2014')
% TBsimulation_jan23('.', 'b01',160, 20000,  '-r70','LEbuilder', 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'PPM0p8_2014')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extras %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% sendmail('suensze@gmail.com','TBsimulation Run','simulation run is done!!')
