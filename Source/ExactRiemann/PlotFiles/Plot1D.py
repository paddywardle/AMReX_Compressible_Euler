import matplotlib.pyplot as plt
import numpy as np

def ReadFile(filename):
     
    data = []
    with open(filename) as f:
        for line in f.readlines():
            line = line.strip().split()
            line = [eval(i) for i in line]
            data.append(line)
    
    return np.asarray(data)

def Plot1D(TestNum, col):

    ExactFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_400cells.dat"
    NoAMRFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_NoAMR"

    ExactData = ReadFile(ExactFilename)
    NoAMR100 = ReadFile(NoAMRFilename + "_100cells.dat")
    NoAMR200 = ReadFile(NoAMRFilename + "_200cells.dat")
    NoAMR400 = ReadFile(NoAMRFilename + "_400cells.dat")

    fig, ax = plt.subplots(figsize=(13,8))
    
    ax.scatter(NoAMR100[:,0], NoAMR100[:,col], c="b", marker="x", label="100 cells")
    ax.scatter(NoAMR200[:,0], NoAMR200[:,col], c="c", marker="+", label="200 cells")
    ax.scatter(NoAMR400[:,0], NoAMR400[:,col], facecolors='none', edgecolors='r', label="400 cells")
    ax.plot(ExactData[:,0], ExactData[:,col], "k", label="Exact")
    ax.legend()
    ax.set_xlabel("X-Domain")

    if TestNum == "Test1":
        ax.set_title("Test 1")
    elif TestNum == "Test2":
        ax.set_title("Test 2")
    elif TestNum == "Test3":
        ax.set_title("Test 3")
    elif TestNum == "Test4":
        ax.set_title("Test 4")
    elif TestNum == "Test5":
        ax.set_title("Test 5")

    if col == 1:
        ax.set_ylabel("Density")
        plt.savefig(f"../RiemannExactResults/{TestNum}/plots/{TestNum}_density.png")
    elif col == 2:
        ax.set_ylabel("Velocity")
        plt.savefig(f"../RiemannExactResults/{TestNum}/plots/{TestNum}_velocity.png")
    elif col == 3:
        ax.set_ylabel("Pressure")
        plt.savefig(f"../RiemannExactResults/{TestNum}/plots/{TestNum}_pressure.png")
    elif col == 4:
        ax.set_ylabel("Specific Internal Energy")
        plt.savefig(f"../RiemannExactResults/{TestNum}/plots/{TestNum}_internal_energy.png")

    plt.close()

def Subplot1D(TestNum):

    ExactFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_400cells.dat"
    NoAMRFilename = f"../RiemannExactResults/{TestNum}/{TestNum}_NoAMR"

    ExactData = ReadFile(ExactFilename)
    NoAMR100 = ReadFile(NoAMRFilename + "_100cells.dat")
    NoAMR200 = ReadFile(NoAMRFilename + "_200cells.dat")
    NoAMR400 = ReadFile(NoAMRFilename + "_400cells.dat")

    fig, ax = plt.subplots(2, 2, figsize=(13,11))
    
    ax[0,0].scatter(NoAMR100[:,0], NoAMR100[:,1], c="b", marker="x", label="100 cells")
    ax[0,0].scatter(NoAMR200[:,0], NoAMR200[:,1], c="c", marker="+", label="200 cells")
    ax[0,0].scatter(NoAMR400[:,0], NoAMR400[:,1], facecolors='none', edgecolors='r', label="400 cells")
    ax[0,0].plot(ExactData[:,0], ExactData[:,1], "k", label="Exact")
    ax[0,0].legend(borderpad=2, fontsize="large")
    ax[0,0].set_xlabel("X-Domain")
    ax[0,0].set_ylabel("Density")

    ax[0,1].scatter(NoAMR100[:,0], NoAMR100[:,2], c="b", marker="x", label="100 cells")
    ax[0,1].scatter(NoAMR200[:,0], NoAMR200[:,2], c="c", marker="+", label="200 cells")
    ax[0,1].scatter(NoAMR400[:,0], NoAMR400[:,2], facecolors='none', edgecolors='r', label="400 cells")
    ax[0,1].plot(ExactData[:,0], ExactData[:,2], "k", label="Exact")
    ax[0,1].legend(borderpad=2, fontsize="large")
    ax[0,1].set_xlabel("X-Domain")
    ax[0,1].set_ylabel("Velocity")

    ax[1,0].scatter(NoAMR100[:,0], NoAMR100[:,3], c="b", marker="x", label="100 cells")
    ax[1,0].scatter(NoAMR200[:,0], NoAMR200[:,3], c="c", marker="+", label="200 cells")
    ax[1,0].scatter(NoAMR400[:,0], NoAMR400[:,3], facecolors='none', edgecolors='r', label="400 cells")
    ax[1,0].plot(ExactData[:,0], ExactData[:,3], "k", label="Exact")
    ax[1,0].legend(borderpad=2, fontsize="large")
    ax[1,0].set_xlabel("X-Domain")
    ax[1,0].set_ylabel("Pressure")

    ax[1,1].scatter(NoAMR100[:,0], NoAMR100[:,4], c="b", marker="x", label="100 cells")
    ax[1,1].scatter(NoAMR200[:,0], NoAMR200[:,4], c="c", marker="+", label="200 cells")
    ax[1,1].scatter(NoAMR400[:,0], NoAMR400[:,4], facecolors='none', edgecolors='r', label="400 cells")
    ax[1,1].plot(ExactData[:,0], ExactData[:,4], "k", label="Exact")
    ax[1,1].legend(borderpad=2, fontsize="large")
    ax[1,1].set_xlabel("X-Domain")
    ax[1,1].set_ylabel("Pressure")
    ax[1,1].set_ylabel("Specific Internal Energy")

    plt.tight_layout()
    
    plt.savefig(f"../RiemannExactResults/{TestNum}/plots/{TestNum}_subplots.png")

    plt.close()


if __name__ == "__main__":

    tests = ["Test1", "Test2", "Test3", "Test4", "Test5"]

    for test in tests:
        Plot1D(test, 1)
        Plot1D(test, 2)
        Plot1D(test, 3)
        Plot1D(test, 4)
        Subplot1D(test)
    
