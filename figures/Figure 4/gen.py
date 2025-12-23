import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import ScalarFormatter, FuncFormatter, SymmetricalLogLocator
import numpy as np

def matplotlib_config(df: pd.DataFrame) -> None:
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size']   = 14
    global fig, unique_name, color_map
    fig = plt.figure(figsize = (18, 18))
    fig.set_dpi(300)
    unique_name = df['Method_name_simple'].unique()
    unique_name.sort()
    color_map = dict(zip(unique_name, mcolors.TABLEAU_COLORS))
    plt.rc('text', usetex = False)
    # plt.rc('font', family = 'Times New Roman')


def graph_MPoX(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'mpox']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    # display(df1)
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]


    ax = fig.add_subplot(321)

    scatter = ax.scatter(df1['SPview'], df1['Time'], s = df1['Maxmemory'] / 16384, alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 12, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names)


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_title('MPoX')
    ax.text(-0.03, 1.05, 'a', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)

def graph_Variola(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'Variola_']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    # display(df1)
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]

    ax = fig.add_subplot(322)

    scatter = ax.scatter(df1['SPview'], df1['Time'], s = df1['Maxmemory'] ** 0.5 * 5, alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 12, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names)


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_xlim(-3200, -2100)
    ax.set_ylim(1.5, 100)
    ax.set_title('Variola virus', fontdict = {'fontstyle': 'italic'})
    ax.text(-0.03, 1.05, 'b', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)

def graph_Mycoplasma(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'Mycoplasma_']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    # display(df1)
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]


    ax = fig.add_subplot(323)

    scatter = ax.scatter(df1['SPview'], df1['Time'], s = df1['Maxmemory'] ** 0.5 , alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 12, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names)


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_title('Mycoplasma', fontdict = {'fontstyle': 'italic'})
    ax.text(-0.03, 1.05, 'c', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)

def graph_Streptococcus(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'Streptococcus_']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    # display(df1)
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]


    ax = fig.add_subplot(324)

    scatter = ax.scatter(df1['SPview'], df1['Time'], s = np.sqrt(df1['Maxmemory']) * 2, alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 12, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "{:.2f}M".format(x / 10 ** 6) if abs(x) >= 10 ** 6 else "{:.1f}".format(x)))
    # ax.ticklabel_format(axis = 'x', style = 'plain')


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_title('Streptococcus', fontdict = {'fontstyle': 'italic'})
    ax.text(-0.03, 1.05, 'd', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)

def graph_ecoli(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'ecoli_']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]


    ax = fig.add_subplot(326)

    df1 = df1.sort_values(by = 'SPview', ascending = False)
    scatter = ax.scatter(df1['SPview'], df1['Time'], s = df1['Maxmemory'] / 16384, alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 12, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names) 


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_xscale('symlog', linthresh = 1e6)
    ax.set_yscale('log')
    xlocator = SymmetricalLogLocator(transform = ax.xaxis.get_transform(), subs = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 15.0])
    ax.xaxis.minorticks_on()
    ax.grid(axis = 'x', which = 'minor', linestyle = '-', linewidth = 0.4, color = 'gray')
    # xlocator.view_limits(vmin = min(list(df1['SPview'])), vmax = max(list(df1['SPview'])))
    ax.xaxis.set_minor_locator(xlocator)
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda x, _: "{:.1f}M".format(x / 10 ** 6) if abs(x) >= 10 ** 6 else "{:.0f}".format(x)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "{:.1f}M".format(x / 10 ** 6) if abs(x) >= 10 ** 6 else "{:.0f}".format(x)))
    ax.set_title('Escherichia coli', fontdict = {'fontstyle': 'italic'})
    ax.text(-0.03, 1.05, 'f', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)

def graph_Neisseria_meningitidis(df: pd.DataFrame) -> None:
    df1 = df[df['Test_data'] == 'Neisseria_meningitidis_']
    del df1['Q_raw'], df1['TC_raw'], df1['Time_view'], df1['Time_raw'], df1['Memory_raw'], df1['Unnamed: 0']
    # display(df1)
    colors = df1['Method_name_simple'].map(color_map)
    graph_map = {}
    for item in df1['Method_name_simple'].to_list():
        graph_map[item] = color_map[item]


    ax = fig.add_subplot(325)

    scatter = ax.scatter(df1['SPview'], df1['Time'], s = df1['Maxmemory'] / 16384, alpha = 0.6, c = colors, edgecolors = 'w', linewidth = 1.5, marker = 'o')

    for i in range(len(df1)):
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], "{:.2f} GB".format(df1['Maxmemory'].iloc[i] / 1048576), fontsize = 12, ha = 'center', va = 'bottom')
        plt.text(df1['SPview'].iloc[i], df1['Time'].iloc[i], '+', fontsize = 14, ha = 'center', va = 'center')


    legend_handles_names = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 12, label = name.replace('_', '-'), alpha = 0.6)
                            for name, color in graph_map.items()]
    legend_names = ax.legend(handles = legend_handles_names, title = "Method Name", loc = "upper left")
    ax.add_artist(legend_names)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "{:.0f}M".format(x / 10 ** 6) if abs(x) >= 10 ** 6 else "{:.1f}".format(x)))


    ax.set_xlabel('SP')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_title('Neisseria meningitidis', fontdict = {'fontstyle': 'italic'})
    ax.text(-0.03, 1.05, 'e', transform = ax.transAxes, fontsize = 14, fontweight = 'bold', ha = 'center', va = 'center')
    ax.grid(True)


def graph_final() -> None:
    plt.subplots_adjust(wspace = 0.1, hspace = 0.2)
    plt.tight_layout()
    plt.subplots_adjust()
    plt.plot()
    plt.savefig('Figure4.svg')


if __name__ == '__main__':
    df = pd.read_csv('data.csv')
    matplotlib_config(df)
    data_list = list(set(list(df["Test_data"])))
    print(data_list)
    if 'mpox' in data_list: graph_MPoX(df)
    else: raise TypeError
    if 'Variola_' in data_list: graph_Variola(df)
    else: raise TypeError
    if 'Mycoplasma_' in data_list: graph_Mycoplasma(df)
    else: raise TypeError
    if 'Streptococcus_' in data_list: graph_Streptococcus(df)
    else: raise TypeError
    if 'ecoli_' in data_list: graph_ecoli(df)
    else: raise TypeError
    if 'Neisseria_meningitidis_' in data_list: graph_Neisseria_meningitidis(df)
    else: raise TypeError
    graph_final()
