import numpy as np
from scipy.spatial.distance import pdist, squareform
import scipy
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial.distance import pdist, squareform
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as mpl
new_rc_params = {'text.usetex': False,
                 "svg.fonttype": 'none'
                 }
mpl.rcParams.update(new_rc_params)
import os
import collections
import seaborn as sns
from operator import itemgetter

def branching_func(list_lines, coord_list ):
    """function defining connected paths on a graph"""
    # TODO: improve readability / efficiency of the implementation
    coord_set = set(coord_list)
    def path_finder(set_koords, sub_path, to_find):
        list_to_find = []
        def append_sublist(sub_list, append_list):
            if sub_list[::-1] not in append_list and sub_list not in append_list:
                append_list.append(sub_list)
        if (to_find[0], to_find[1]) in set_koords:
            set_koords.remove((to_find[0], to_find[1]))
            sub_path.append(to_find)
            [append_sublist(sub_list, list_to_find) for sub_list in list_lines if to_find in sub_list]
            for i in list_to_find:
                if i[0] == to_find:
                    path_finder(set_koords, sub_path, i[1])
                else:
                    path_finder(set_koords, sub_path, i[0])
    paths = []
    for koords in coord_list:
        if tuple(koords) in coord_set:
            to_find = [float(koords[0]), float(koords[1])]
            sub_path = []
            path_finder(coord_set, sub_path, to_find)
            paths.append(sub_path)
    keys = []
    values = []
    values2 = []
    for i in paths:
        for j in i:
            keys.append(tuple(j))
            values.append([len(i)])
            values2.append(i)
    dicts = dict(zip(keys, values))
    dicts_con = dict(zip(keys, values2))
    return keys, dicts, dicts_con

def plot_graph(KD, Radiation):
    files = os.listdir()
    path = 'Data' + '\\' + rf'{Radiation}Gy\{KD}'
    files = os.listdir(path)
    files_xls = [f for f in files if f[-4:] == 'xlsx']
    files_xls = np.sort(files_xls)
    data = np.random.choice(files_xls)

    df = pd.read_excel(path + "\\" + data)
    df_1 = df.loc[df['Slice'] == 1]
    df_2 = df.loc[df['Slice'] == 2]
    df_3 = df.loc[df['Slice'] == 3]
    dfs = [df_1, df_2, df_3]

    distance_um = 30

    for counter_slice, df_i in enumerate(dfs):
        X = df_i['X']
        Y = df_i['Y']

        xmin, xmax = np.min(X), np.max(X)
        ymin, ymax = np.min(Y), np.max(Y)

        X, Y = X - xmin, Y - ymin
        coord_list = list(df_i.iloc[:, [1, 2]].itertuples(index=False, name=None))
        set_koords = set(coord_list)
        coord_array = np.asarray(coord_list)

        # Create a matrix with the corresponding euclidean distances for each value pair
        distmat = squareform(pdist(coord_array, 'euclidean'))

        n = np.sort(distmat, axis=1)
        n1 = np.argsort(distmat, axis=1)
        number_of_neigh = np.sum(np.where(n > distance_um, 0, 1), axis=1)
        number_of_neigh = number_of_neigh - np.ones(np.shape(number_of_neigh))
        n = np.where(n > distance_um, np.NaN, n1)
        # return a dictionary which has the coordinate as the key and the number of neighbors as the value

        # slicing with first index corresponding to row index
        X_n = []
        for i in np.arange(len(X)):
            # get all non
            X_n.append(n[i][np.logical_not(np.isnan(n[i]))])
            if np.isnan(n[:, i]).all():
                meigh_max = i
        # print(X_n)
        X_n = np.asarray(X_n, dtype=object)
        coordinates_connections = np.zeros((len(X), meigh_max, 2, 2))

        for i in np.arange(len(X)):
            for j in np.arange(len(X_n[i])):
                coordinates_connections[i, j, :, 0] = np.array(
                    [coord_array[i, :][0], coord_array[int(X_n[i][j]), :][0]])
                coordinates_connections[i, j, :, 1] = np.array(
                    [coord_array[i, :][1], coord_array[int(X_n[i][j]), :][1]])

        # create line artists
        lines_ = coordinates_connections.reshape((len(X) * meigh_max, 2, 2))
        connections_list = []

        # check and exclud empty entries and self connecting lines
        for i in lines_:
            if i.any() and [[i[0][0], i[0][1]] != [i[1][0], i[1][1]]]:
                connections_list.append(i.tolist())
    ### PLOT OF SINGLE GRAPHS ###
    lines = LineCollection(connections_list, color='crimson')
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.scatter(coord_array[:, 0], coord_array[:, 1], c='black', s=10)
    ax.add_artist(lines)
    plt.show()

def read_celldata():
    KDs = ['AOX1', 'BIRC5', 'C12ORF56', 'CENPE', 'CLSPN', 'E2F2', 'FGF1', 'FOXG1', 'ID1', 'ID3', 'PSMC3IP', 'PTK2', 'RAD51AP1', 'RTKN2', 'shCtrlOlga', 'Untrans', 'XRCC3']
    radiation = [0, 4]
    # KDs = ['AOX1', 'BIRC5', 'C12ORF56', 'CENPE']
    df_all = pd.DataFrame()

    for KD in KDs:
        for rad in radiation:
            path = 'Data' + '\\' + rf'{rad}Gy\{KD}'
            files = os.listdir(path)
            files_xls = [f for f in files if f[-4:] == 'xlsx']
            files_xls = np.sort(files_xls)

            for data in files_xls:
                df = pd.read_excel(path + "\\" + data)
                df_1 = df.loc[df['Slice'] == 1]
                df_2 = df.loc[df['Slice'] == 2]
                df_3 = df.loc[df['Slice'] == 3]
                dfs = [df_1, df_2, df_3]

                distance_um = 30

                for counter_slice, df_i in enumerate(dfs):
                    X = df_i['X']
                    Y = df_i['Y']

                    xmin, xmax = np.min(X), np.max(X)
                    ymin, ymax = np.min(Y), np.max(Y)

                    X, Y = X-xmin, Y-ymin
                    coord_list= list(df_i.iloc[:, [1,2]].itertuples(index=False, name=None))
                    set_koords = set(coord_list)
                    coord_array = np.asarray(coord_list)

                    # Create a matrix with the corresponding euclidean distances for each value pair
                    distmat = squareform(pdist(coord_array, 'euclidean'))


                    n = np.sort(distmat, axis=1)
                    n1 = np.argsort(distmat, axis=1)
                    number_of_neigh = np.sum(np.where(n > distance_um, 0, 1), axis=1)
                    number_of_neigh = number_of_neigh - np.ones(np.shape(number_of_neigh))
                    n = np.where(n>distance_um, np.NaN, n1)
                    # return a dictionary which has the coordinate as the key and the number of neighbors as the value

                    # slicing with first index corresponding to row index
                    X_n = []
                    for i in np.arange(len(X)):
                        # get all non
                        X_n.append(n[i][np.logical_not(np.isnan(n[i]))])
                        if np.isnan(n[:, i]).all():
                            meigh_max = i
                    # print(X_n)
                    X_n = np.asarray(X_n, dtype=object)
                    coordinates_connections = np.zeros((len(X), meigh_max, 2, 2))

                    for i in np.arange(len(X)):
                        for j in np.arange(len(X_n[i])):
                            coordinates_connections[i, j, :, 0] = np.array([coord_array[i,:][0], coord_array[int(X_n[i][j]), :][0]])
                            coordinates_connections[i, j, :, 1] = np.array([coord_array[i,:][1], coord_array[int(X_n[i][j]), :][1]])

                    # create line artists
                    lines_ = coordinates_connections.reshape((len(X)*meigh_max, 2, 2))
                    connections_list = []

                    # check and exclud empty entries and self connecting lines
                    for i in lines_:
                        if i.any() and [[i[0][0], i[0][1]] != [i[1][0], i[1][1]]]:
                            connections_list.append(i.tolist())

                    #### Dictionaries #####
                    # Dictionary: neighbors, keys:coordinate. values: coordinates of neighbor nodes
                    dict_neigh = collections.defaultdict(list)
                    for i in connections_list:
                        if (i[0][0], i[0][1]) != (i[1][0], i[1][1]):
                            dict_neigh[(i[0][0], i[0][1])].append((i[1][0], i[1][1]))

                    # Dictionary: Isolated node, keys:coordinates. values:boolean value whether node is isolated or not
                    dict_iso = collections.defaultdict(list)
                    # Dictionary: Distance to core, keys:coordinates. values:euclidean distance to core
                    dict_dist = collections.defaultdict(list)
                    # Dictionary: Clustercoeff, keys:coordinates. values:cluster coeff
                    dict_cluster = collections.defaultdict(list)

                    for i in np.arange(len(coord_list)):
                        counter = 0
                        coord = coord_list[i]
                        dict_dist[coord].append(np.sqrt((coord[0]**2 + abs(ymax - coord[1])**2)))
                        if len(dict_neigh[coord]) == 0:
                            dict_iso[coord].append(1)
                        else:
                            dict_iso[coord].append(0)

                        if len(dict_neigh[coord]) in [0, 1]:
                            dict_cluster[coord].append(0)

                        else:
                            for j in (dict_neigh[coord]):
                                list_neigh = dict_neigh[coord].copy()
                                list_neigh.remove(j)
                                for k in dict_neigh[j]:
                                    if k in dict_neigh[coord]:
                                        counter += 1

                            # C = 2n / (k * (k-1)) with n=number of edges; k=number of possible edges
                            # scaled with a factor of 1/2 to avoid counting edges twice
                            dict_cluster[coord].append(counter / (len(dict_neigh[coord]) * (len(dict_neigh[coord]) - 1)))
                    # Dictionary: connected nodes, keys:coord and values:number of continously connected node
                    _, dict_con, _ = branching_func(connections_list, coord_list)

                    ### PLOT OF SINGLE GRAPHS ###
                    # lines = LineCollection(connections_list, color='crimson')
                    # fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                    # ax.scatter(coord_array[:, 0], coord_array[:, 1], c='black', s=10)
                    # ax.add_artist(lines)
                    # plt.show()
                    # print(df_all)

                    ### Calculate global values ####
                    connected = list(dict_con.values())
                    cluster = list(dict_cluster.values())
                    iso = list(dict_iso.values())
                    neigh = list(dict_neigh.values())
                    neigh_number = [len(item) for item in neigh]

                    df_all = df_all.append({'Label': data, 'KD': KD, 'Radiation':rad, 'Slice': counter_slice, 'Cellnumber': len(X)
                                            , 'Connectedness':np.mean(connected)/len(X), 'Clustercoefficient':np.mean(cluster)
                                            , 'Perc_of_Isolated_Nodes':np.mean(iso), 'Giantcluster':np.max(connected)/len(X)
                                            , 'Average_Number_of_Neighbors': np.mean(neigh_number)
                                            }, ignore_index=True)
    df_all.to_csv('Analysis_all_Graphs.csv',index=False, header=True)

def plot_graph(Radiation, KD):
    files = os.listdir()
    path = 'Data' + '\\' + rf'{Radiation}Gy\{KD}'
    files = os.listdir(path)
    files_xls = [f for f in files if f[-4:] == 'xlsx']
    files_xls = np.sort(files_xls)
    data = np.random.choice(files_xls)

    df = pd.read_excel(path + "\\" + data)
    df_1 = df.loc[df['Slice'] == 1]
    df_2 = df.loc[df['Slice'] == 2]
    df_3 = df.loc[df['Slice'] == 3]
    dfs = [df_1, df_2, df_3]

    distance_um = 30

    X = df_1['X']
    Y = df_1['Y']

    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)

    X, Y = X - xmin, Y - ymin
    coord_list = list(df_1.iloc[:, [1, 2]].itertuples(index=False, name=None))
    set_koords = set(coord_list)
    coord_array = np.asarray(coord_list)

    # Create a matrix with the corresponding euclidean distances for each value pair
    distmat = squareform(pdist(coord_array, 'euclidean'))

    n = np.sort(distmat, axis=1)
    n1 = np.argsort(distmat, axis=1)
    number_of_neigh = np.sum(np.where(n > distance_um, 0, 1), axis=1)
    number_of_neigh = number_of_neigh - np.ones(np.shape(number_of_neigh))
    n = np.where(n > distance_um, np.NaN, n1)
    # return a dictionary which has the coordinate as the key and the number of neighbors as the value

    # slicing with first index corresponding to row index
    X_n = []
    for i in np.arange(len(X)):
        # get all non
        X_n.append(n[i][np.logical_not(np.isnan(n[i]))])
        if np.isnan(n[:, i]).all():
            meigh_max = i
    # print(X_n)
    X_n = np.asarray(X_n, dtype=object)
    coordinates_connections = np.zeros((len(X), meigh_max, 2, 2))

    for i in np.arange(len(X)):
        for j in np.arange(len(X_n[i])):
            coordinates_connections[i, j, :, 0] = np.array(
                [coord_array[i, :][0], coord_array[int(X_n[i][j]), :][0]])
            coordinates_connections[i, j, :, 1] = np.array(
                [coord_array[i, :][1], coord_array[int(X_n[i][j]), :][1]])

    # create line artists
    lines_ = coordinates_connections.reshape((len(X) * meigh_max, 2, 2))
    connections_list = []

    # check and exclud empty entries and self connecting lines
    for i in lines_:
        if i.any() and [[i[0][0], i[0][1]] != [i[1][0], i[1][1]]]:
            connections_list.append(i.tolist())

    # Calculate Delaunay graph using Scipy
    tri = scipy.spatial.Delaunay(coord_array)

    ### PLOT OF SINGLE GRAPHS ###
    lines = LineCollection(connections_list, color='crimson')
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    axs[0].scatter(coord_array[:, 0], coord_array[:, 1], c='black', s=10)
    axs[0].add_artist(lines)
    axs[1].triplot(coord_array[:, 0], coord_array[:, 1], tri.simplices, color='crimson')
    axs[1].plot(coord_array[:, 0], coord_array[:, 1], 'o', color='black', markersize=2)
    axs[1].set_title("Delaunay Graph", fontsize=15)
    axs[0].set_title("Geometric Graph", fontsize=15)
    # plt.show()


if __name__ == "__main__":
    plot_graph(0, 'Untrans')
    # plot_graph(KD='Untrans', Radiation=0)