import matplotlib.pyplot as plt
import numpy as np


def plot_map(sol_array, X, Y, title=None):
    plt.figure()
    plt.imshow(sol_array, extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)], origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlim(0, np.max(X))
    plt.ylim(0, np.max(Y))
    plt.xlabel("x")
    plt.ylabel("y")
    if title:
        plt.title(title)

def plot_3D(sol_array, X, Y, title=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, sol_array, cmap='viridis')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if title:
        ax.set_zlabel(title)

def plot_field(U, V, X, Y, title=None):
    plt.figure()
    plt.quiver(X, Y, U, V, linewidth=None)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Y), np.max(Y))
    plt.xlabel("x")
    plt.ylabel("y")
    if title:
        plt.title(title)
    # Show plot with grid
    plt.grid()

def plot_xprofile(sol_arrays, x, ylabel=None, title=None, label_lst=None):
    plt.figure(figsize=(10, 5))

    for i, arr in enumerate(sol_arrays):
        if label_lst:
            plt.plot(x, arr, label=label_lst[i], marker="o")
        else:
            plt.plot(x, arr, marker="o") 
    plt.xlabel("x")
    plt.ylabel(ylabel) 
    # plt.ylim(0, 6)
    plt.title(title)
    if label_lst:
        plt.legend()
    plt.grid(True)
   
def plot_yprofile(sol_array, y, xlabel=None, title=None):
    plt.figure(figsize=(5, 6))
    plt.plot(sol_array, y, marker="o")
    if xlabel:
        plt.xlabel(xlabel)
    plt.ylabel("y")
    plt.title(title)
    plt.grid(True)
