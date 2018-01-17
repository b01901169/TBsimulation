import matlab.engine
import os
import uptakePolicy

if __name__ == "__main__":
    eng_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/TBsimulationCode"
    eng = matlab.engine.start_matlab()
    eng.cd(eng_path)

    NumPpl = 2000
    yrsOfAnalysis = 25

    Group1 = eng.TBsimulation_jan23('.','p01', 130+yrsOfAnalysis, NumPpl, '-r70','loadBurnIn', 2018,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,'base', 'customized',
                                    matlab.double(uptakePolicy.finerUptakeUrbanKnowledge), matlab.int8(uptakePolicy.finerUptakeAgeBracs))

