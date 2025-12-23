import pandas as pd
import copy
from collections import ChainMap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import ScalarFormatter, FuncFormatter, SymmetricalLogLocator
import numpy as np
from typing import Callable


def matplotlib_config(df: pd.DataFrame) -> None:
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size']   = 14
    global fig, unique_name, colors
    fig = plt.figure(figsize = (20, 16))
    fig.set_dpi(300)
    unique_name = df['Method_name_simple'].unique()
    unique_name.sort()
    colors = plt.get_cmap('tab10')(range(10))
    colors2 = plt.get_cmap('Dark2_r')(range(10))
    colors = np.concatenate((colors, colors2))
    plt.rc('text', usetex = False)

def graph_final(fig_name: str) -> None:
    plt.subplots_adjust(wspace = 0.1, hspace = 0.2)
    plt.tight_layout()
    plt.subplots_adjust()
    plt.plot()
    plt.savefig(fig_name)

def read_sim_csv() -> None:
    global df
    df_raw = pd.read_csv("data.csv", encoding = 'gb2312')
    del df_raw['SAMEMAlign_Version'], df_raw['Sub_BlockSize']
    df = df_raw

def graph_sim(this_data_name: str, same_graph: bool = False) -> None:
    # graph matrix is:
    #          mafft wmsa abpoa FMAlign
    # tree-Q    []    []   []    []
    # tree-TC   []    []   []    []
    # sim-Q     []    []   []    []
    # sim-TC    []    []   []    []

    Q  = 'Q_view'
    TC = 'TC_view'
    data_name      = this_data_name + '_like_diff_'
    data_name_sim  = data_name + 'similarity'

    block_keys_sim = [80, 85, 90]
    block_list_sim = [i for i in range(len(block_keys_sim))]

    def valid_keys_in_data(df: pd.DataFrame, block_keys: list) -> dict:
        keys = df['BlockSize'].unique().tolist()
        ans  = {}
        for key in keys:
            ok = True
            for id_ in block_keys:
                # print(id_)
                # display(df[(df['BlockSize'] == key)]['Similarity'])
                if df[(df['BlockSize'] == key) & (df['Similarity'] == id_)]['Vaild_tests'].iloc[0] != 9:
                    ok = False
            if ok:
                ans.update({str(key): {}})
                for id_ in block_keys:
                    # print(id_)
                    # display(df[(df['BlockSize'] == key)]['Similarity'])
                    ans[str(key)].update({str(id_): {Q: df[(df['BlockSize'] == key) & (df['Similarity'] == id_)]['Q'].iloc[0],
                                                     TC: df[(df['BlockSize'] == key) & (df['Similarity'] == id_)]['TC'].iloc[0]}})
        return ans
    def dict_to_list_sim(d: dict, insert_keys: list, block_size: list, data: str) -> list:
        L = []
        for jtem in block_size:
            this_cycle = []
            for item in insert_keys:
                this_cycle.append(d[str(item)][str(jtem)][data])
            L.append(this_cycle)
        return L
    # Sim Data
    mafft_low  = df[(df['Method_name_simple'] == 'fftns1') & (df['Test_data'] == data_name_sim)]
    mafft_high = df[(df['Method_name_simple'] == 'fftns1_highsim') & (df['Test_data'] == data_name_sim)]
    wmsa_low   = df[(df['Method_name_simple'] == 'wmsa') & (df['Test_data'] == data_name_sim)]
    wmsa_high  = df[(df['Method_name_simple'] == 'wmsa_highsim') & (df['Test_data'] == data_name_sim)]
    abpoa_low  = df[(df['Method_name_simple'] == 'abpoa') & (df['Test_data'] == data_name_sim)]
    abpoa_high = df[(df['Method_name_simple'] == 'abpoa_highsim') & (df['Test_data'] == data_name_sim)]

    # check the data is valid in all thresholds
    d_mafft_low  = valid_keys_in_data(mafft_low, block_keys_sim)
    d_mafft_high = valid_keys_in_data(mafft_high, block_keys_sim)
    d_wmsa_low   = valid_keys_in_data(wmsa_low, block_keys_sim)
    d_wmsa_high  = valid_keys_in_data(wmsa_high, block_keys_sim)
    d_abpoa_low  = valid_keys_in_data(abpoa_low, block_keys_sim)
    d_abpoa_high = valid_keys_in_data(abpoa_high, block_keys_sim)
    same_keys = set(d_mafft_low.keys()) & set(d_mafft_high.keys()) & set(d_wmsa_low.keys()) & set(d_wmsa_high.keys()) & set(d_abpoa_low.keys()) & set(d_abpoa_high.keys())
    same_keys2 = []
    for item in same_keys:
        if int(item) <= 5000:
            same_keys2.append(int(item))
    same_keys = sorted(list(same_keys2))
    block_list_sim = [i for i in range(len(same_keys))]
    # Q
    Q_mafft_low  = dict_to_list_sim(d_mafft_low, same_keys, block_keys_sim, Q)
    Q_mafft_high = dict_to_list_sim(d_mafft_high, same_keys, block_keys_sim, Q)
    Q_wmsa_low   = dict_to_list_sim(d_wmsa_low, same_keys, block_keys_sim, Q)
    Q_wmsa_high  = dict_to_list_sim(d_wmsa_high, same_keys, block_keys_sim, Q)
    Q_abpoa_low  = dict_to_list_sim(d_abpoa_low, same_keys, block_keys_sim, Q)
    Q_abpoa_high = dict_to_list_sim(d_abpoa_high, same_keys, block_keys_sim, Q)
    legend = []
    # ax = plt.subplot(4, 3, 7)
    if not same_graph:
        ax = plt.subplot(221)
    else:
        ax = plt.subplot(223)
    for i, item in enumerate(Q_mafft_low):
        ax.errorbar(block_list_sim, item, color = colors[i], linestyle = '-', fmt = '-o')
        legend.append('deMEM-FFTNS1-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(Q_mafft_high):
        ax.errorbar(block_list_sim, item, color = colors[i], linestyle = '--', fmt = '-o')
        legend.append('deMEM-FFTNS1-H {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(Q_wmsa_low):
        ax.errorbar(block_list_sim, item, color = colors[i + 3], linestyle = '-', fmt = '-o')
        legend.append('deMEM-WMSA-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(Q_wmsa_high):
        ax.errorbar(block_list_sim, item, color = colors[i + 3], linestyle = '--', fmt = '-o')
        legend.append('deMEM-WMSA-H {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(Q_abpoa_low):
        ax.errorbar(block_list_sim, item, color = colors[i + 6], linestyle = '-', fmt = '-o')
        legend.append('deMEM-abPOA-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(Q_abpoa_high):
        ax.errorbar(block_list_sim, item, color = colors[i + 6], linestyle = '--', fmt = '-o')
        legend.append('deMEM-abPOA-H {}{} similarity'.format(block_keys_sim[i], '%'))
    # ax.legend(legend, loc = 'lower right')
    ax.set_yscale('symlog')
    ax.set_ylim(0.9, 0.99)
    ax.set_ylabel('Q')
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.set_xticks(block_list_sim)
    ax.set_xticklabels(same_keys)
    # TC
    TC_mafft_low  = dict_to_list_sim(d_mafft_low, same_keys, block_keys_sim, TC)
    TC_mafft_high = dict_to_list_sim(d_mafft_high, same_keys, block_keys_sim, TC)
    TC_wmsa_low   = dict_to_list_sim(d_wmsa_low, same_keys, block_keys_sim, TC)
    TC_wmsa_high  = dict_to_list_sim(d_wmsa_high, same_keys, block_keys_sim, TC)
    TC_abpoa_low  = dict_to_list_sim(d_abpoa_low, same_keys, block_keys_sim, TC)
    TC_abpoa_high = dict_to_list_sim(d_abpoa_high, same_keys, block_keys_sim, TC)
    legend = []
    # ax = plt.subplot(4, 3, 7)
    if not same_graph:
        ax = plt.subplot(222)
    else:
        ax = plt.subplot(224)
    for i, item in enumerate(TC_mafft_low):
        ax.errorbar(block_list_sim, item, color = colors[i], linestyle = '-', fmt = '-o')
        legend.append('deMEM-FFTNS1-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(TC_mafft_high):
        ax.errorbar(block_list_sim, item, color = colors[i], linestyle = '--', fmt = '-o')
        legend.append('deMEM-FFTNS1-H {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(TC_wmsa_low):
        ax.errorbar(block_list_sim, item, color = colors[i + 3], linestyle = '-', fmt = '-o')
        legend.append('deMEM-WMSA-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(TC_wmsa_high):
        ax.errorbar(block_list_sim, item, color = colors[i + 3], linestyle = '--', fmt = '-o')
        legend.append('deMEM-WMSA-H {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(TC_abpoa_low):
        ax.errorbar(block_list_sim, item, color = colors[i + 6], linestyle = '-', fmt = '-o')
        legend.append('deMEM-abPOA-L {}{} similarity'.format(block_keys_sim[i], '%'))
    for i, item in enumerate(TC_abpoa_high):
        ax.errorbar(block_list_sim, item, color = colors[i + 6], linestyle = '--', fmt = '-o')
        legend.append('deMEM-abPOA-H {}{} similarity'.format(block_keys_sim[i], '%'))
    ax.legend(legend, bbox_to_anchor = (1, 0.75))
    ax.set_yscale('log')
    ax.set_ylabel('TC')
    ax.set_ylim(0.5, 0.9)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.set_xticks(block_list_sim)
    ax.set_xticklabels(same_keys)



if __name__ == '__main__':
    read_sim_csv()
    matplotlib_config(df)
    graph_sim('mt')
    graph_sim('sars_cov_2', True)
    graph_final('Figure7.svg')
