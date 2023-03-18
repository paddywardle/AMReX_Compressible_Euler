import numpy as np

def ReadFile(filename):
     
    data = []
    with open(filename) as f:
        for line in f.readlines():
            line = line.strip().split()
            line = [eval(i) for i in line]
            data.append(line)
    
    return np.asarray(data)

def WriteData(TestNum, results):

    OutputFile = f"../RiemannExactResults/{TestNum}/ConvData/{TestNum}_results.txt"

    with open(OutputFile, "w") as f:
        for exp in results:
            for val in exp:
                f.write(str(val) + " ")
            f.write("\n")
    

def ConvAnalysis(TestNum):

    results = [["Version", "Norm", "Density", "Velocity", "Pressure", "Specific Internal Energy"]]

    ExactFilename = f"../RiemannExactResults/{TestNum}/{TestNum}"
    NoAMRFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_NoAMR"
    AMRFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_AMR"
    NumCells = ["100cells", "200cells", "400cells"]

    for cells in NumCells:

        ExactData = ReadFile(ExactFilename+"_"+cells+".dat")
        AMRData = ReadFile(NoAMRFilename+"_"+cells+".dat")
        DataDiffs = AMRData[:,1:] - ExactData[:,1:]

        L1_norm = np.linalg.norm(DataDiffs, ord=1, axis=0)/len(AMRData)
        L2_norm = np.linalg.norm(DataDiffs, ord=2, axis=0)/len(AMRData)
        Linf_norm = np.linalg.norm(DataDiffs, ord=np.inf, axis=0)/len(AMRData)

        results.append([cells, "L1", L1_norm[0], L1_norm[1], L1_norm[2], L1_norm[3]])
        results.append([cells, "L2", L2_norm[0], L2_norm[1], L2_norm[2], L2_norm[3]])
        results.append([cells, "Linf", Linf_norm[0], Linf_norm[1], Linf_norm[2], Linf_norm[3]])

    ExactData = ReadFile(ExactFilename+"_"+NumCells[0]+".dat")
    AMRData = ReadFile(AMRFilename+"_"+NumCells[0]+".dat")
    DataDiffs = AMRData[:,1:] - ExactData[:,1:]

    L1_norm = np.linalg.norm(DataDiffs, ord=1, axis=0)/len(AMRData)
    L2_norm = np.linalg.norm(DataDiffs, ord=2, axis=0)/len(AMRData)
    Linf_norm = np.linalg.norm(DataDiffs, ord=np.inf, axis=0)/len(AMRData)

    results.append(["100_AMR", "L1", L1_norm[0], L1_norm[1], L1_norm[2], L1_norm[3]])
    results.append(["100_AMR", "L2", L2_norm[0], L2_norm[1], L2_norm[2], L2_norm[3]])
    results.append(["100_AMR", "Linf", Linf_norm[0], Linf_norm[1], Linf_norm[2], Linf_norm[3]])    

    return results

if __name__ == "__main__":

    tests = ["Test1", "Test2", "Test3", "Test4", "Test5"]

    for test in tests:
        results = ConvAnalysis(test)
        WriteData(test, results)
        
