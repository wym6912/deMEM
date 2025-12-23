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
    global fig, unique_name, colors
    fig = plt.figure(figsize = (12, 6))
    fig.set_dpi(300)
    unique_name = df['Method_name_simple'].unique()
    unique_name.sort()
    colors = plt.get_cmap('tab10')(range(10))
    colors2 = plt.get_cmap('Dark2_r')(range(2))
    colors = np.concatenate((colors, colors2))
    plt.rc('text', usetex = False)

def graph_final() -> None:
    plt.subplots_adjust(wspace = 0.1, hspace = 0.2)
    plt.tight_layout()
    plt.subplots_adjust()
    plt.plot()
    plt.savefig('Figure6.svg')

def read_BestQTC_csv() -> None:
    '''
    Q_dict  = {'abpoa': {'RNA-8191': 300}, 'fftns1': {'RNA-2047': 500, 'RNA-4095': 100, 'RNA-8191': 20}, 
                'wmsa': {'RNA-1023': 700, 'RNA-2047': 1500, 'RNA-255': 1000, 'RNA-511': 1500, 'RNA-8191': 400}, 
                'abpoa_highsim': {}, 'fftns1_highsim': {}, 'wmsa_highsim': {}}
    TC_dict = {'abpoa': {}, 'fftns1': {'RNA-1023': 1500}, 
                'wmsa': {'RNA-1023': 700, 'RNA-2047': 1500, 'RNA-255': 300, 'RNA-4095': 1500, 'RNA-511': 1500, 'RNA-8191': 1100}, 
                'abpoa_highsim': {}, 'fftns1_highsim': {}, 'wmsa_highsim': {}}
    '''
    global ext_dict, ext_dict_view
    ext_dict = {'abpoa': {}, 'fftns1': {'RNA-2047': [500], 'RNA-4095': [100]}, 'wmsa': {}, 
                'abpoa_highsim': {}, 'fftns1_highsim': {'RNA-1023': [1500]}, 
                'wmsa_highsim': {'RNA-1023': [700], 'RNA-2047': [1500], 'RNA-255': [1000, 300], 'RNA-511': [1500]}}
    ext_dict_view = {'abpoa': [2000], 'fftns1': [100, 500, 2000], 'wmsa': [2000], 'abpoa_highsim': [2000], 'fftns1_highsim': [1500, 2000], 
                        'wmsa_highsim': [300, 700, 1000, 1500, 2000]}

def read_RNA_csv() -> None:
    df_raw = pd.read_csv("data.csv", encoding = 'gb2312')
    global df
    df_raw['Test_data_base'] = df_raw['Test_data'].str.split('-').str[0]
    df = df_raw[df_raw['Test_data_base'] == 'RNA']
    df = df[df['Test_data'] != 'RNA-8191'] # not drawn 8191
    del df['SAMEMAlign_Version'], df['SP'], df['SPview'], df['Sub_BlockSize'], df['Test_data_base']
    df_final = pd.DataFrame()
    for val in set(df['Test_data']):
        for method in ext_dict:
            try:
                for ext_block_val in ext_dict[method][val]:
                    selected_rows = df[(df['Test_data'] == val) & (df['BlockSize'] == ext_block_val) & (df['Vaild_tests'] == 10) & (df['Method_name_simple'] == method)]
                    df_final = pd.concat([df_final, selected_rows], ignore_index = True)
                    print(f'Method {method} has {ext_block_val} better than raw method in test {val}')
            except:
                pass
    # display(df_final)
    selected_rows = df[(df['Test_data'] != 'RNA-8191') & (df['BlockSize'] == 2000)]
    df_final = pd.concat([df_final, selected_rows], ignore_index = True)
    # display(df_final)
    df = df_final
    # display(df)

def check_valid() -> bool:
    ans = True
    select_rows = df[df['BlockSize'] != 2000]
    for i in range(len(select_rows)):
        now_value = select_rows.iloc[i]
        raw_value = df[(df['BlockSize'] == 2000) & (df['Test_data'] == now_value['Test_data']) & (df['Method_name_simple'] == now_value['Method_name_simple'])].iloc[0]
        if not (now_value['Q'] > raw_value['Q'] or now_value['TC'] > raw_value['TC']): 
            print("Error!! Test case {} is lower than {} - Q: {} vs {}, TC: {} vs {}".format(now_value['Test_data'], raw_value['Test_data'], now_value['Q'], raw_value['Q'], now_value['TC'], raw_value['TC']))
            ans = False
    return ans

def graph_RNA() -> None:
    # graph matrix is:
    #       mafft wmsa abpoa
    # Q     []    []   []
    # TC    []    []   []
    df_raw = pd.read_csv("data.csv", encoding = 'gb2312')
    # Re-add datas
    df = pd.DataFrame()
    for item in ext_dict_view:
        for blocksize in ext_dict_view[item]:
            selected_rows = df_raw[(df_raw['Test_data'] != 'RNA-8191') & (df_raw['BlockSize'] == blocksize) & (df_raw['Vaild_tests'] == 10) & (df_raw['Method_name_simple'] == item)]
            df = pd.concat([df, selected_rows], ignore_index = True)
    del df['SAMEMAlign_Version'], df['SP'], df['SPview'], df['Sub_BlockSize']
    mafft_low  = df[df['Method_name_simple'] == 'fftns1']
    mafft_high = df[df['Method_name_simple'] == 'fftns1_highsim']
    wmsa_low   = df[df['Method_name_simple'] == 'wmsa']
    wmsa_high  = df[df['Method_name_simple'] == 'wmsa_highsim']
    abpoa_low  = df[df['Method_name_simple'] == 'abpoa']
    abpoa_high = df[df['Method_name_simple'] == 'abpoa_highsim']
    Q       = 'Q_view'
    TC      = 'TC_view'
    RNA255  = 'RNA-255'
    RNA511  = 'RNA-511'
    RNA1023 = 'RNA-1023'
    RNA2047 = 'RNA-2047'
    RNA4095 = 'RNA-4095'

    block_keys = [RNA255, RNA511, RNA1023, RNA2047, RNA4095]
    block_list = [i for i in range(len(block_keys))]

    getQTC = lambda df, key, data_name, blocksize: df[(df['BlockSize'] == blocksize) & (df['Test_data'] == data_name)][key].iloc[0] if not df[(df['BlockSize'] == blocksize) & (df['Test_data'] == data_name)][key].empty else None

    def RNA_table_to_list(df: pd.DataFrame, data: list, key: str, blocksize: int, data_rule: Callable[[pd.Series, str, str, int], (float | None)]) -> dict:
        d = {}
        for i in data:
            value = data_rule(df, key, i, blocksize)
            if value:
                d.update({i: value.split("±")})
        # d = df[(df['BlockSize'].isin(data)) & (df['Test_data'] == data_name)][key].to_list()
        return d
    def draw_RNA_QTC(ax: plt.Axes, group_ids: list, data: list, title: str, legend: list, drawlegend: bool, log: bool = False) -> None:
        real_id_gap = 0.2
        group_gap = 1.6

        positions = []
        for i in range(len(data)):
            real_id_positions = [pos + (i * real_id_gap) for pos in range(len(group_ids))]
            adjusted_positions = [real_id_pos + group_gap * i for i, real_id_pos in enumerate(real_id_positions)]
            positions.append(adjusted_positions)

        for i, dataset in enumerate(data):
            box = ax.boxplot(dataset, patch_artist=True, positions=positions[i], widths=0.15, labels=group_ids, flierprops = dict(marker = 'o', color = 'none', markersize = 0))
            
            for j, patch in enumerate(box['boxes']):
                patch.set_facecolor(colors[i])

            for whisker, cap, median in zip(box['whiskers'], box['caps'], box['medians']):
                whisker.set_color('black')
                cap.set_color('black')
                median.set_color('black')


        xticks_positions = [i + (len(data) - 1) * real_id_gap / 2 + i * group_gap for i in range(len(group_ids))]

        ax.set_xticks(xticks_positions)
        ax.set_xticklabels(group_ids, rotation = 45, ha = 'right')

        plt.subplots_adjust(left = 0.05, right = 0.95)
        if drawlegend:
            handles = [plt.Line2D([0], [0], marker = 's', color = 'w', markerfacecolor = colors[i], markersize = 10, label = legend[i]) for i in range(len(data))]
            ax.legend(handles = handles, loc = 'upper right', fontsize = 8)
        if log:
            ax.set_yscale('logit')

        ax.set_title(title)
        ax.set_ylabel(title)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
        if log:
            ax.set_yscale('logit')
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
            ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.3f}"))
            ax.minorticks_off()

    # Q
    legend = []
    graph_all = []
    for blocksize in ext_dict_view['fftns1']:
        if blocksize == 2000: continue
        graph_raw = []
        for item in block_keys:
            graph_raw.append(mafft_low[(mafft_low['BlockSize'] == blocksize) & (mafft_low['Test_data'] == item)]['Q_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        legend.append('deMEM-FFTNS1-L Block size = {}'.format(blocksize))
    for blocksize in ext_dict_view['fftns1_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(mafft_high[(mafft_high['BlockSize'] == blocksize) & (mafft_high['Test_data'] == item)]['Q_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: legend.append('deMEM-FFTNS1-H Block size = {}'.format(blocksize))
        else: legend.append('FFTNS1')
    for blocksize in ext_dict_view['wmsa']:
        if blocksize == 2000: continue
        graph_raw = []
        for item in block_keys:
            graph_raw.append(wmsa_low[(wmsa_low['BlockSize'] == blocksize) & (wmsa_low['Test_data'] == item)]['Q_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        legend.append('deMEM-WMSA-L Block size = {}'.format(blocksize))
    for blocksize in ext_dict_view['wmsa_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(wmsa_high[(wmsa_high['BlockSize'] == blocksize) & (wmsa_high['Test_data'] == item)]['Q_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: legend.append('deMEM-WMSA-H Block size = {}'.format(blocksize))
        else: legend.append('WMSA')
    for blocksize in ext_dict_view['abpoa_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(abpoa_high[(abpoa_high['BlockSize'] == blocksize) & (abpoa_high['Test_data'] == item)]['Q_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: raise TypeError
        else: legend.append('abPOA')
    
    ax = plt.subplot(121)

    draw_RNA_QTC(ax, block_keys, graph_all, 'Q', legend, False, True)

    # TC
    legend = []
    graph_all = []
    for blocksize in ext_dict_view['fftns1']:
        if blocksize == 2000: continue
        graph_raw = []
        for item in block_keys:
            graph_raw.append(mafft_low[(mafft_low['BlockSize'] == blocksize) & (mafft_low['Test_data'] == item)]['TC_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        legend.append('deMEM-FFTNS1-L Block size = {}'.format(blocksize))
    for blocksize in ext_dict_view['fftns1_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(mafft_high[(mafft_high['BlockSize'] == blocksize) & (mafft_high['Test_data'] == item)]['TC_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: legend.append('deMEM-FFTNS1-H Block size = {}'.format(blocksize))
        else: legend.append('FFTNS1')
    for blocksize in ext_dict_view['wmsa']:
        if blocksize == 2000: continue
        graph_raw = []
        for item in block_keys:
            graph_raw.append(wmsa_low[(wmsa_low['BlockSize'] == blocksize) & (wmsa_low['Test_data'] == item)]['TC_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        legend.append('deMEM-WMSA-L Block size = {}'.format(blocksize))
    for blocksize in ext_dict_view['wmsa_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(wmsa_high[(wmsa_high['BlockSize'] == blocksize) & (wmsa_high['Test_data'] == item)]['TC_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: legend.append('deMEM-WMSA-H Block size = {}'.format(blocksize))
        else: legend.append('WMSA')
    for blocksize in ext_dict_view['abpoa_highsim']:
        graph_raw = []
        for item in block_keys:
            graph_raw.append(abpoa_high[(abpoa_high['BlockSize'] == blocksize) & (abpoa_high['Test_data'] == item)]['TC_raw'].iloc[0])
        graph_all.append([eval(_) for _ in graph_raw])
        if blocksize != 2000: raise TypeError
        else: legend.append('abPOA')
    
    ax = plt.subplot(122)
    # display(graph_all)
    draw_RNA_QTC(ax, block_keys, graph_all, 'TC', legend, True)



if __name__ == '__main__':
    read_BestQTC_csv()
    read_RNA_csv()
    if not check_valid():
        print("Error! Data not OK")
    else:
        matplotlib_config(df)
        graph_RNA()
        graph_final()
