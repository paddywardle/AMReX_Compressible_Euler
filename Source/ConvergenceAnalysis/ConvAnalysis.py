import numpy as np

def ReadFile(filename):
     
    data = []
    with open(filename) as f:
        for line in f.readlines():
            line = line.strip().split()
            line = [eval(i) for i in line]
            data.append(line)
    
    return np.asarray(data)

def ConvAnalysis(FilenameExact, FilenameAMR):

    ExactData = ReadFile(FilenameExact)
    AMRData = ReadFile(FilenameAMR)
    DataDiffs = AMRData[:,1:] - ExactData[:,1:]

    L1_norm = np.linalg.norm(DataDiffs, ord=1, axis=0)/len(AMRData)
    L2_norm = np.linalg.norm(DataDiffs, ord=2, axis=0)/len(AMRData)
    Linf_norm = np.linalg.norm(DataDiffs, ord=np.inf, axis=0)/len(AMRData)

    return L1_norm, L2_norm, Linf_norm

if __name__ == "__main__":

    L1_norm, L2_norm, Linf_norm = ConvAnalysis("Sim100_Test1.dat", "Test1_100cells.dat")

    print(L1_norm, L2_norm, Linf_norm)