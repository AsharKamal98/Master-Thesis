import numpy as np
import math

from scipy.stats import qmc
import matplotlib.pyplot as plt
import pandas as pd

import subprocess
import time
import sys
from tqdm import tqdm

start_time = time.time()

#=============================================================== USER INPUT ========================================================
model = "TSM"
v = 246 # Fix! Read from LH
num_in_param = 13  # Fix! Can be found from df
num_free_param = 10  # Fix!
exp_num_data_points = 12
num_data_points = 2**(exp_num_data_points)

d = { \
    'Input parameter name': ['lam1','lam2','lam3','lam4','lam5','lam6','lam7','lam8','lam9','lam10','lam11','mT','mS','mH','mC','mN1','mN2'], \
    'LesHouches number': [1,2,3,4,5,6,7,8,9,10,11,12,13,None,None,None,None], \
    'Range start' : [-1.15, -1.0, 0, 0, -1.6, -1.15, -0.5, None, None, -0.29, -0.29, None, None, 125.25, 3300, 2200, 4400], \
    'Range end' : [1.15, 1.0, 0, 0, 1.6, 1.15, 0.5, None, None, 0.29, 0.29, None, None, 125.25, 3399, 2299, 4500], \
    'Dependence' :[None, None, None, None, None, None, None, '(1/v**2) * 2**(3/2) * (mC-mN1)**(1/2) * (mN2-mC)**(1/2)', 'mH/(2*v**2)', None, None, '(1/4) * (2*mC - lam10*v**2)', '(1/2) * (mN1 + mN2 - mC - lam7*v**2)', None, None, None, None] \
    }
df = pd.DataFrame(data=d)
#print(df)


# PATHS
LesHouches_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/LesHouches.in.{}".format(model)
SPheno_spc_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/SPheno.spc.{}".format(model)

HB_script_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsBounds-5.3.2beta/RunHiggsBounds.sh"
HB_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsBounds-5.3.2beta"
HB_output_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/HiggsBounds_results.dat"

HS_script_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsSignals-2.2.3beta/RunHiggsSignals.sh"
HS_path = "/home/etlar/m22_ashar/.Mathematica/Applications/HiggsSignals-2.2.3beta"
HS_output_path = "/home/etlar/m22_ashar/.Mathematica/Applications/SPheno-4.0.5/HiggsSignals_results.dat"

#===================================================================================================================================

def SearchGrid(n):
    DataFile_InParam = open("DataFile_InParam", "w")
    DataFile_InParam.writelines(f'{"INPUT PARAMETERS / PATTERNS"} \n {"Lambda1":<{10}} {"Lambda2":<{10}} {"Lambda3":<{10}} {"Lambda4":<{10}} {"Lambda5":<{10}} {"Lambda6":<{10}} {"TanBeta":<{10}} \n')
    DataFile_InParam.close()

    DataFile_Labels = open("DataFile_Labels", "w")
    DataFile_Labels.writelines(f'{"OUTPUT / LABELS"} \n {"T parameter":<{20}} {"S parameter":<{20}} {"U parameter":<{20}} {"HB Result":<{20}} {"HS chi^2()":<{20}} {"HS chi^2()":<{20}} {"Real masses":<{20}} \n')
    DataFile_Labels.close()

    DataFile_Masses = open("DataFile_Masses", "w")
    DataFile_Masses.writelines(f'{"PARTICLE MASSES"} \n {"mH(1)":<{20}} {"mH(2)":<{20}} {"mH(3)":<{20}} \n')
    DataFile_Masses.close()


    input_samples = Sampling()
    
    df_free1 = df[~df['Range start'].isna()]
    df_free2 = df_free1[df_free1['Range start'] != df_free1['Range end']]

    #num_data_points = 1
    for i in tqdm(range(num_data_points)):
        subprocess.run(["rm", "-f", SPheno_spc_path])

        d = dict((df_free2['Input parameter name'].to_numpy()[j], input_samples[i,j]) for j in range(num_free_param))
        d['v'] = v
        d['lam3'] = 0
        d['lam4'] = 0
        d['mH'] = float(df['Range start'][13])

        lam8 = eval(df['Dependence'][7], d)
        lam9 = eval(df['Dependence'][8], d)
        mT = eval(df['Dependence'][11], d)
        mS = eval(df['Dependence'][12], d)
        Analysis([input_samples[i,0], input_samples[i,1], 0, 0, input_samples[i,2], input_samples[i,3], input_samples[i,4], lam8, lam9, input_samples[i,5], input_samples[i,6], mT, mS])
        #if i%100 == 0:
        #print("Finished {}/{} iterations".format(i+1, num_data_points))
        #print(input_samples[0,-3:])


    return "Done!"

def Sampling():
    df_free1 = df[~df['Range start'].isna()]
    df_free2 = df_free1[df_free1['Range start'] != df_free1['Range end']]

    free_param_range = df_free2[['Range start', 'Range end']].to_numpy()
    sampler = qmc.Sobol(d=num_free_param)
    sample = sampler.random_base2(m=exp_num_data_points)
    input_samples = qmc.scale(sample, free_param_range[:,0], free_param_range[:,1])

    return input_samples


def Analysis(in_param_list):

    if len(in_param_list) != num_in_param:
        print("Number of free parameters does not match given input parameters. \nExiting code!")
        sys.exit()

#---------------------------------------------------------------------------------------
                                # === INPUT PARMETERS === #
    # Find MINPAR block in LesHouches file
    InputFile = open(LesHouches_path, "r")
    l = InputFile.readlines()
    MINPARindex = [idx for idx, s in enumerate(l) if 'Block MINPAR' and '# Input' in s][0]

    # Define input parameters in LesHouches file
    for i in range(num_in_param):
        l[MINPARindex+1+i] = " {}   {}     # {}\n".format(i+1, in_param_list[i], l[MINPARindex+1+i].split()[-1])

    InputFile = open(LesHouches_path, "w")
    InputFile.writelines(l)
    InputFile.close()
    
#----------------------------------------------------------------------------------------
                                    # === RUN/READ HEP PACKAGES ===#
    RunSPheno(model)
    try:
        spheno_output1, spheno_output2 = ReadSPheno()
        real_masses = 1

        RunHiggsBounds()
        higgsbounds_output = ReadHiggsBounds()
        RunHiggsSignals()
        higgssignals_output = ReadHiggsSignals()

    except:
        print("Not valid parameters, next iteration!")        
        spheno_output1 = ['Negative mass square(s)']
        spheno_output2 = ['X', 'X', 'X']
        real_masses = 0
        higgsbounds_output = 'X'
        higgssignals_output = ['X','X']

#---------------------------------------------------------------------------------------
                                    # === WRITE DATA FILES === #
    DataFile_InParam = open("DataFile_InParam", "a")
    for i in range(num_in_param):
        DataFile_InParam.writelines(f'{round(in_param_list[i],5):<{13}}')
    DataFile_InParam.writelines('\n')
    DataFile_InParam.close()

    DataFile_Labels = open("DataFile_Labels", "a")
    for i in range(3):
        DataFile_Labels.writelines(f'{spheno_output2[i]:<{20}}')
    DataFile_Labels.writelines(f'{higgsbounds_output:<{20}} {higgssignals_output[0]:<{20}} {higgssignals_output[1]:<{20}} {real_masses:<{20}} \n')
    DataFile_Labels.close()

    DataFile_Masses = open("DataFile_Masses", "a")
    for i in range(len(spheno_output1)):
        DataFile_Masses.writelines(f'{spheno_output1[i]:<{20}}')
    DataFile_Masses.writelines('\n')
    DataFile_Masses.close()


    return None



def RunSPheno(model):
    subprocess.run(["./bin/SPheno{}".format(model)], stdout=subprocess.DEVNULL)
    #subprocess.run(["./bin/SPheno{}".format(model)])
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


print(SearchGrid(exp_num_data_points))


#DataFile = open("DataFile_Labels", "r")
#l = DataFile.readlines()
#pattern = (l[2].split())
#print(pattern)
#if '0' in pattern:
#    print("inside!")
#DataFile.close()

print("The script took {} seconds to run".format(time.time()-start_time))
