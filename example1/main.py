import numpy as np
import matplotlib.pyplot as plt

def plot_results(name):
    results = np.loadtxt(f"{name}.txt")

    eigvals = np.unique(results[:,2])
    assert(len(np.unique(eigvals)) == 3)

    R = np.array([[0,-1],[1,0]])
    
    res1 = results[results[:,2] == eigvals[0]][:,:2]
    res2 = results[results[:,2] == eigvals[1]][:,:2]
    res3 = results[results[:,2] == eigvals[2]][:,:2]

    res1 = res1 @ R.T
    res2 = res2 @ R.T
    res3 = res3 @ R.T

    colors = [ '#5F9ED1', '#C85200', '#CFCFCF'];

    fig = plt.figure()
    plt.plot(res1[:,0], res1[:,1], 's', color=colors[0], label=f"$\lambda_0 = {eigvals[0]}$", markersize=0.1)
    plt.plot(res2[:,0], res2[:,1], 's', color=colors[1], label=f"$\lambda_1 = {eigvals[1]}$", markersize=0.1)
    plt.plot(res3[:,0], res3[:,1], 's', color=colors[2], label=f"$\lambda_2 = {eigvals[2]}$", markersize=0.1)
    legend = plt.legend(markerscale=100)
    plt.xticks([], [])
    plt.yticks([], [])
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    plt.savefig(f"./figures/{name}.png", dpi=450, format="png");

plot_results("classic_0.100000")
plot_results("classic_0.980000")

plot_results("complex_0.100000")
plot_results("complex_0.980000")

