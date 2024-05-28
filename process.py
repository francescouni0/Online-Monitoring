import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import matplotlib as mpl
import scienceplots
import seaborn as sns
import os
import scipy.stats as stats


root_files = []
csv_files = []

# Get the list of root and csv files in the build/train directory
for file in os.listdir("build/train"):
    if file.endswith(".root"):
        root_files.append(file)
    elif file.endswith(".csv"):
        csv_files.append(file)

# Sort the lists alphabetically
root_files.sort()
csv_files.sort()

# Process the files in pairs
for root_file, csv_file in zip(root_files, csv_files):
    root_file_path = os.path.join("build/train", root_file)
    csv_file_path = os.path.join("build/train", csv_file)

    file = uproot.open(root_file_path)
    tree = file["Hits"]
    data = tree.arrays(library="pd")

    df_gamma = pd.read_csv(csv_file_path, delimiter=',', header=1, names=['Energy', 'X', 'Y', 'Z','pX', 'pY', 'pZ','evt'])

    theta = np.arccos(df_gamma['pX'] / np.sqrt(df_gamma['pX']**2 + df_gamma['pY']**2 + df_gamma['pZ']**2))
    df_gamma['theta'] = theta

    df_gamma = df_gamma[np.isclose(df_gamma['theta'], 1, atol=0.5*1)]

    hist2d, x_edges, y_edges = np.histogram2d(df_gamma['X'], df_gamma['Y'], bins=(25, 1), range=[[-149.9, 149.9], [-100, 100]])

    slice_index = 0

    # Extract the slice from the 2D histogram
    hist_slice = hist2d[:, slice_index]
    hist_slice = hist_slice / np.sum(hist_slice)    

    bin_centers_x = (x_edges[:-1] + x_edges[1:]) / 2

    #plt.hist(bin_centers_x, weights=hist_slice, color='blue', bins=25, histtype='step', density=True)
#
    ## Add labels and title
    #plt.xlabel('X[mm]')
    #plt.ylabel('Frequency')
    #plt.title('Densità gamma')
    #plt.show()

    ######################################################################################################

    edep = data['EDep']
    xloc = data['xloc'] 

    hist, x_edges = np.histogram(xloc, weights=edep, bins=25)
    bin_centers = (x_edges[:-1] + x_edges[1:]) / 2

    #plt.hist(bin_centers, weights=hist, color='goldenrod', bins=25, histtype='step', range=(-150, 150))
    #plt.xlabel('X[mm]')
    #plt.ylabel('Absolute Dose [a.u]')
    #plt.xlim(-150, 150)
#
    #plt.show()


####################################################
a=np.hstack((hist_slice, hist))


plt.hist(bin_centers,a[0:24], bins=25, histtype='step')
plt.hist(bin_centers,a[25:49], bins=25, histtype='step')
plt.show()

np.savetxt('build/train/merged.csv', a, delimiter=',')