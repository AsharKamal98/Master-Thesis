import numpy as np
import subprocess
import time
import sys

start_time = time.time()


model = "TSM"
num_in_param = 13

LesHouches_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/LesHouches.in.{}".format(model)
SPheno_spc_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/SPheno.spc.{}".format(model)

HB_script_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsBounds-5.3.2beta/RunHiggsBounds.sh"
HB_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsBounds-5.3.2beta"
HB_output_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/HiggsBounds_results.dat"

HS_script_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsSignals-2.2.3beta/RunHiggsSignals.sh"
HS_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsSignals-2.2.3beta"
HS_output_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/HiggsSignals_results.dat"


def SearchGrid(n):
    DataFile_InParam = open("DataFile_InParam", "w")
    DataFile_InParam.writelines(f'{"INPUT PARAMETERS / PATTERNS"} \n {"Lambda1":<{10}} {"Lambda2":<{10}} {"Lambda3":<{10}} {"Lambda4":<{10}} {"Lambda5":<{10}} {"Lambda6":<{10}} {"TanBeta":<{10}} \n')
    DataFile_InParam.close()

    DataFile_Labels = open("DataFile_Labels", "w")
    DataFile_Labels.writelines(f'{"OUTPUT / LABELS"} \n {"T parameter":<{20}} {"S parameter":<{20}} {"U parameter":<{20}} {"HB Result":<{20}} {"HS chi^2()":<{20}} {"HS chi^2()":<{20}} \n')
    DataFile_Labels.close()

    DataFile_Masses = open("DataFile_Masses", "w")
    DataFile_Masses.writelines(f'{"PARTICLE MASSES"} \n {"mH(1)":<{20}} {"mH(2)":<{20}} {"mH(3)":<{20}} \n')
    DataFile_Masses.close()


    # Set input parameters
    lam4 = 0.10
    lam5 = 0.10
    lam6 = 0.10
    lam7 = 0.10
    lam8 = 0.10
    lam9 = 0.10
    lam10 = 0.10
    lam11 = 0.10
    mT = 0.10
    mS = 0.10


    lam1_range = np.linspace(-0.01, 1.15 ,n)
    lam2_range = np.linspace(-0.05, 0.29 ,n)
    lam3_range = np.linspace(-0.3, 1.3 ,n)
    lam4_range = np.linspace(-0.7, -0.3 ,n)
    lam5_range = np.linspace(-0.1, 0.6 ,n)
    lam6_range = np.linspace(-0.01, 1.15 ,n)
    lam7_range = np.linspace(-0.05, 0.29 ,n)
    lam8_range = np.linspace(-0.3, 1.3 ,n)
    lam9_range = np.linspace(-0.7, -0.3 ,n)
    lam10_range = np.linspace(-0.1, 0.6 ,n)
    lam11_range = np.linspace(-0.1, 0.6 ,n)
    mT_range = np.linspace(-0.1, 0.6 ,n)
    mS_range = np.linspace(-0.1, 0.6 ,n)

    #lam1, lam2, lam3, lam4, lam5, m12, tanbeta = np.meshgrid(lam1_range, lam2_range, lam3_range, lam4_range, lam5_range, m12_range, tanbeta_range)
    lam1, lam2, lam3   = np.meshgrid(lam1_range, lam2_range, lam3_range)
    #print(lam1.flatten())

    forloop_len = len(lam1.flatten())
    for i in range(forloop_len):
        subprocess.run(["rm", "-f", SPheno_spc_path])
        Analysis([lam1.flatten()[i],lam2.flatten()[i],lam3.flatten()[i],lam4,lam5,lam6,lam7,lam8,lam9,lam10,lam11,mT,mS])
        print("Finished {}/{} iterations".format(i+1, forloop_len))


    return "Done!"




def Analysis(in_param_list):

    if len(in_param_list) != num_in_param:
        print("Number of free parameters does not match given input parameters. \nExiting code!")
        sys.exit()

    DataFile_InParam = open("DataFile_InParam", "a")
    for i in range(num_in_param):
        DataFile_InParam.writelines(f'{in_param_list[i]:<{13}}')
    DataFile_InParam.writelines('\n')
    DataFile_InParam.close()


                                # === INPUT PARMETERS === #
    # Find MINPAR block in LesHouches file
    InputFile = open(LesHouches_path, "r")
    l = InputFile.readlines()
    MINPARindex = [idx for idx, s in enumerate(l) if 'Block MINPAR' and '# Input' in s][0]

    # Define input parameters in LesHouches file
    for i in range(num_in_param):
        l[MINPARindex+1+i] = " {}   {}     # {}\n".format(i, in_param_list[i], l[MINPARindex+1+i].split()[-1])


    InputFile = open(LesHouches_path, "w")
    InputFile.writelines(l)
    InputFile.close()

    
#----------------------------------------------------------------------------------------
                                    # === RUN/READ PACKAGES ===#
    RunSPheno(model)
    DataFile_Labels = open("DataFile_Labels", "a")

    try:
        spheno_output1, spheno_output2 = ReadSPheno()
    except:
        print("SPheno does not like these paramater choices, let's go to the next iteration!")
        DataFile_Labels.writelines(f'{"Not valid point!":{30}} \n')
        DataFile_Labels.close()
        return None


    RunHiggsBounds()
    higgsbounds_output = ReadHiggsBounds()

    RunHiggsSignals()
    higgssignals_output = ReadHiggsSignals() 

#---------------------------------------------------------------------------------------

                                    # === WRITE DATA FILE === #
    for i in range(3):
        DataFile_Labels.writelines(f'{spheno_output2[i]:<{20}}')
    DataFile_Labels.writelines(f'{higgsbounds_output:<{20}} {higgssignals_output[0]:<{20}} {higgssignals_output[1]:<{20}} \n')
    DataFile_Labels.close()

    DataFile_Masses = open("DataFile_Masses", "a")
    for i in range(4):
        DataFile_Masses.writelines(f'{spheno_output1[i]:<{20}}')
    DataFile_Masses.writelines('\n')
    DataFile_Masses.close()


    return None



def RunSPheno(model):
    subprocess.run(["./bin/SPheno{}".format(model)], stdout=subprocess.DEVNULL)

    return None

def ReadSPheno():
    OutputFile = open(SPheno_spc_path, "r")
    l = OutputFile.readlines()

    index1 = [idx for idx, s in enumerate(l) if 'Block MASS' in s][0]
    spheno_output1 = [l[index1+i+2].split()[1] for i in range(4)] # May use list comprehenseion over for loop!

    index2 = [idx for idx, s in enumerate(l) if 'Block SPhenoLowEnergy' in s][0]
    spheno_output2 = [l[index2+1+i].split()[1] for i in range(3)] # May use list comprehenseion over for loop!

    OutputFile.close()
    return spheno_output1, spheno_output2


def RunHiggsBounds():
    subprocess.run(["bash", HB_script_path, HB_path])
    return None

def RunHiggsSignals():
    subprocess.run(["bash", HS_script_path, HS_path])
    return None

def ReadHiggsBounds():
    OutputFile = open(HB_output_path)
    l = OutputFile.readlines()
    index1 = [idx for idx, s in enumerate(l) if '#cols' in s][0]
    index2 = l[index1].split().index("HBresult")
    higgsbounds_output = l[index1+2].split()[index2-1]
    OutputFile.close()
    return higgsbounds_output

def ReadHiggsSignals():
    OutputFile = open(HS_output_path)
    l = OutputFile.readlines()
    index1 = [idx for idx, s in enumerate(l) if '#cols:' in s][0]
    index2 = l[index1].split().index("csq(mu)")
    higgssignals_output = [l[index1+2].split()[index2-1+i] for i in range(2)]
    OutputFile.close()
    return higgssignals_output


print(SearchGrid(2))


print("The script took {} seconds to run".format(time.time()-start_time))
