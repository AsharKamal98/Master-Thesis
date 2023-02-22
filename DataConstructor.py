import numpy as np
import subprocess
import time
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
    DataFile = open("DataFile", "w")
    DataFile.writelines(f'{"Lambda1 / Mh(1)":<{30}} {"Lambda2 / Mh(2)":<{30}} {"Lambda3 / Mh(3)":<{30}} {"Lambda4 / Mhplus(1)":<{30}} {"Lambda5 / HBresult":<{30}} {"M12 / chi2(mu)":<{30}} {"TanBeta / chi2(H)":<{30}} \n')
    DataFile.close()

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



#def Analysis(Lambda1Input,Lambda2Input,Lambda3Input,Lambda4Input,Lambda5Input,M12Input,TanBeta):
def Analysis(in_param_list):

    DataFile = open("DataFile", "a")
    #DataFile.writelines(f'{Lambda1Input:<{20}} {Lambda2Input:<{20}} {Lambda3Input:<{20}} {Lambda4Input:<{20}} {Lambda5Input:<{20}} {M12Input:<{20}} {TanBeta:<{20}} \n')
    for i in range(num_in_param):
        DataFile.writelines(f'{in_param_list[i]:<{30}}')
    DataFile.writelines('\n')


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
    try:
        spheno_output = ReadSPheno()
    except:
        print("SPheno does not like these paramater choices, let's go to the next iteration!")
        DataFile.writelines(f'{"Not valid point!":{30}} \n \n')
        DataFile.close()
        return None


    RunHiggsBounds()
    higgsbounds_output = ReadHiggsBounds()

    RunHiggsSignals()
    higgssignals_output = ReadHiggsSignals() 

#---------------------------------------------------------------------------------------

                                    # === WRITE DATA FILE === #
    DataFile.writelines(f'{spheno_output[0]:<{20}} {spheno_output[1]:<{20}} {spheno_output[2]:<{20}} {spheno_output[3]:<{20}} {higgsbounds_output:<{20}} {higgssignals_output[0]:<{20}} {higgssignals_output[1]:<{20}} \n \n')
    DataFile.close()

    return None



def RunSPheno(model):
    subprocess.run(["./bin/SPheno{}".format(model)], stdout=subprocess.DEVNULL)

    return None

def ReadSPheno():
    OutputFile = open(SPheno_spc_path, "r")
    l = OutputFile.readlines()
    index = [idx for idx, s in enumerate(l) if 'Block MASS' in s][0]
    spheno_output = [l[index+i+2].split()[1] for i in range(4)] # May use list comprehenseion over for loop!
    OutputFile.close()
    return spheno_output


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


print(SearchGrid(3))
print("The script took {} seconds to run".format(time.time()-start_time))
