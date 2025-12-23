import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import ScalarFormatter, FuncFormatter, SymmetricalLogLocator
import numpy as np
from typing import Callable

def matplotlib_config(df: pd.DataFrame) -> None:
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 14
    global fig, unique_name, colors
    fig = plt.figure(figsize = (20, 18))
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
    plt.savefig('Figure5.svg')


def mt_graph(df: pd.DataFrame) -> None:
    # graph matrix is:
    #       mafft wmsa abpoa
    # SP     []    []   []
    # Time   []    []   []
    # Space  []    []   []
    mafft_low  = df[df['Method_name_simple'] == 'fftns1']
    mafft_high = df[df['Method_name_simple'] == 'fftns1_highsim']
    wmsa_low   = df[df['Method_name_simple'] == 'wmsa']
    wmsa_high  = df[df['Method_name_simple'] == 'wmsa_highsim']
    abpoa_low  = df[df['Method_name_simple'] == 'abpoa']
    abpoa_high = df[df['Method_name_simple'] == 'abpoa_highsim']
    block_list = list(df['BlockSize'].unique()); block_list.sort()
    SPview = 'SPview'
    Time   = 'Time'
    Memory = 'Maxmemory'
    mt1x   = 'mt1x'
    mt20x  = 'mt20x'

    getSP = lambda df, i, key, data_name: df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].iloc[0] if not df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].empty else None
    getMemory = lambda df, i, key, data_name: df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].iloc[0] if not df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].empty else None
    getTime = lambda df, i, key, data_name: df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].iloc[0] if not df[(df['BlockSize'] == i) & (df['Test_data'] == data_name)][key].empty else None

    def mt_table_to_list(df: pd.DataFrame, data: list, key: str, data_name: str, data_rule: Callable[[pd.Series, int, str, str], (float | None)]) -> dict:
        d = {}
        for i in data:
            value = data_rule(df, i, key, data_name)
            if value:
                d.update({i: float(value)})
        # d = df[(df['BlockSize'].isin(data)) & (df['Test_data'] == data_name)][key].to_list()
        return d
    '''
    SP
    '''
    # MAFFT SP
    d1 = mt_table_to_list(mafft_low, block_list, SPview, mt1x, getSP)
    d2 = mt_table_to_list(mafft_low, block_list, SPview, mt20x, getSP)
    d3 = mt_table_to_list(mafft_high, block_list, SPview, mt1x, getSP)
    d4 = mt_table_to_list(mafft_high, block_list, SPview, mt20x, getSP)
    fftns_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    fftns_list = [i for i in range(len(fftns_keys))]
    p1 = [d1[i] for i in fftns_keys]
    p2 = [d2[i] for i in fftns_keys]
    p3 = [d3[i] for i in fftns_keys]
    p4 = [d4[i] for i in fftns_keys]
    ax = plt.subplot(331)
    ax.errorbar(fftns_list, p2, color = colors[0], fmt = '-o')
    ax.errorbar(fftns_list, p4, color = colors[1], fmt = '-o')
    ax.errorbar(fftns_list, p1, color = colors[2], linestyle = '--', fmt = '-o')
    ax.errorbar(fftns_list, p3, color = colors[3], linestyle = '--', fmt = '-o')
    ax.set_xticks(fftns_list)
    ax.set_xticklabels(fftns_keys, rotation = 45, ha = 'right')
    ax.set_yscale('symlog', linthresh = 250)
    ax.set_ylim(top = -140, bottom = -250)
    ax.set_xlabel('Block Size')
    ax.set_ylabel('SP value')
    ax.set_title('deMEM-FFTNS1')
    ax.legend(['mt20x deMEM-FFTNS1','mt20x deMEM-FFTNS1-H', 'mt1x deMEM-FFTNS1', 'mt1x deMEM-FFTNS1-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'a', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # WMSA SP
    d1 = mt_table_to_list(wmsa_low, block_list, SPview, mt1x, getSP)
    d2 = mt_table_to_list(wmsa_low, block_list, SPview, mt20x, getSP)
    d3 = mt_table_to_list(wmsa_high, block_list, SPview, mt1x, getSP)
    d4 = mt_table_to_list(wmsa_high, block_list, SPview, mt20x, getSP)
    wmsa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    wmsa_list = [i for i in range(len(wmsa_keys))]
    p1 = [d1[i] for i in wmsa_keys]
    p2 = [d2[i] for i in wmsa_keys]
    p3 = [d3[i] for i in wmsa_keys]
    p4 = [d4[i] for i in wmsa_keys]
    ax = plt.subplot(332)
    ax.errorbar(wmsa_list, p2, color = colors[4], fmt = '-o')
    ax.errorbar(wmsa_list, p4, color = colors[5], fmt = '-o')
    ax.errorbar(wmsa_list, p1, color = colors[6], fmt = '-o', linestyle = '--')
    ax.errorbar(wmsa_list, p3, color = colors[7], fmt = '-o', linestyle = '--')
    ax.set_xticks(wmsa_list)
    ax.set_xticklabels(wmsa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('symlog', linthresh = 250)
    ax.set_ylim(top = -140, bottom = -250)
    ax.set_xlabel('Block Size')
    ax.set_title('deMEM-WMSA')
    ax.legend(['mt20x deMEM-WMSA', 'mt20x deMEM-WMSA-H', 'mt1x deMEM-WMSA', 'mt1x deMEM-WMSA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'b', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # abPOA SP
    d1 = mt_table_to_list(abpoa_low, block_list, SPview, mt1x, getSP)
    d2 = mt_table_to_list(abpoa_low, block_list, SPview, mt20x, getSP)
    d3 = mt_table_to_list(abpoa_high, block_list, SPview, mt1x, getSP)
    d4 = mt_table_to_list(abpoa_high, block_list, SPview, mt20x, getSP)
    abpoa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    abpoa_list = [i for i in range(len(abpoa_keys))]
    p1 = [d1[i] for i in abpoa_keys]
    p2 = [d2[i] for i in abpoa_keys]
    p3 = [d3[i] for i in abpoa_keys]
    p4 = [d4[i] for i in abpoa_keys]
    ax = plt.subplot(333)
    ax.errorbar(abpoa_list, p2, color = colors[8], fmt = '-o')
    ax.errorbar(abpoa_list, p4, color = colors[9], fmt = '-o')
    ax.errorbar(abpoa_list, p1, color = colors[10], fmt = '-o', linestyle = '--')
    ax.errorbar(abpoa_list, p3, color = colors[11], fmt = '-o', linestyle = '--')
    ax.set_xticks(abpoa_list)
    ax.set_xticklabels(abpoa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('symlog', linthresh = 250)
    ax.set_ylim(top = -140, bottom = -250)
    ax.set_xlabel('Block Size')
    ax.set_title('deMEM-abPOA')
    ax.legend(['mt20x deMEM-abPOA', 'mt20x deMEM-abPOA-H', 'mt1x deMEM-abPOA', 'mt1x deMEM-abPOA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'c', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    '''
    Time
    '''
    # MAFFT Time
    d1 = mt_table_to_list(mafft_low, block_list, Time, mt1x, getTime)
    d2 = mt_table_to_list(mafft_low, block_list, Time, mt20x, getTime)
    d3 = mt_table_to_list(mafft_high, block_list, Time, mt1x, getTime)
    d4 = mt_table_to_list(mafft_high, block_list, Time, mt20x, getTime)
    p1 = [d1[i] for i in fftns_keys]
    p2 = [d2[i] for i in fftns_keys]
    p3 = [d3[i] for i in fftns_keys]
    p4 = [d4[i] for i in fftns_keys]
    ax = plt.subplot(334)
    ax.errorbar(fftns_list, p2, color = colors[0], fmt = '-o')
    ax.errorbar(fftns_list, p4, color = colors[1], fmt = '-o')
    ax.errorbar(fftns_list, p1, color = colors[2], linestyle = '--', fmt = '-o')
    ax.errorbar(fftns_list, p3, color = colors[3], linestyle = '--', fmt = '-o')
    ax.set_xticks(fftns_list)
    ax.set_xticklabels(fftns_keys, rotation = 45, ha = 'right')
    ax.set_xlabel('Block Size')
    ax.set_ylabel('Time/s')
    ax.set_yscale('log')
    ax.set_ylim(5, 600000)
    # ax.legend(['mt20x deMEM-FFTNS1', 'mt20x deMEM-FFTNS1-H', 'mt1x deMEM-FFTNS1', 'mt1x deMEM-FFTNS1-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    # ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'd', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # WMSA Time
    d1 = mt_table_to_list(wmsa_low, block_list, Time, mt1x, getTime)
    d2 = mt_table_to_list(wmsa_low, block_list, Time, mt20x, getTime)
    d3 = mt_table_to_list(wmsa_high, block_list, Time, mt1x, getTime)
    d4 = mt_table_to_list(wmsa_high, block_list, Time, mt20x, getTime)
    wmsa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    wmsa_list = [i for i in range(len(wmsa_keys))]
    p1 = [d1[i] for i in wmsa_keys]
    p2 = [d2[i] for i in wmsa_keys]
    p3 = [d3[i] for i in wmsa_keys]
    p4 = [d4[i] for i in wmsa_keys]
    ax = plt.subplot(335)
    ax.errorbar(wmsa_list, p2, color = colors[4], fmt = '-o')
    ax.errorbar(wmsa_list, p4, color = colors[5], fmt = '-o')
    ax.errorbar(wmsa_list, p1, color = colors[6], fmt = '-o', linestyle = '--')
    ax.errorbar(wmsa_list, p3, color = colors[7], fmt = '-o', linestyle = '--')
    ax.set_xticks(wmsa_list)
    ax.set_xticklabels(wmsa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('log')
    ax.set_xlabel('Block Size')
    ax.set_ylim(5, 600000)
    # ax.legend(['mt20x deMEM-WMSA', 'mt20x deMEM-WMSA-H', 'mt1x deMEM-WMSA', 'mt1x deMEM-WMSA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'e', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # abPOA Time
    d1 = mt_table_to_list(abpoa_low, block_list, Time, mt1x, getTime)
    d2 = mt_table_to_list(abpoa_low, block_list, Time, mt20x, getTime)
    d3 = mt_table_to_list(abpoa_high, block_list, Time, mt1x, getTime)
    d4 = mt_table_to_list(abpoa_high, block_list, Time, mt20x, getTime)
    abpoa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    abpoa_list = [i for i in range(len(abpoa_keys))]
    p1 = [d1[i] for i in abpoa_keys]
    p2 = [d2[i] for i in abpoa_keys]
    p3 = [d3[i] for i in abpoa_keys]
    p4 = [d4[i] for i in abpoa_keys]
    ax = plt.subplot(336)
    ax.errorbar(abpoa_list, p2, color = colors[8], fmt = '-o')
    ax.errorbar(abpoa_list, p4, color = colors[9], fmt = '-o')
    ax.errorbar(abpoa_list, p1, color = colors[10], fmt = '-o', linestyle = '--')
    ax.errorbar(abpoa_list, p3, color = colors[11], fmt = '-o', linestyle = '--')
    ax.set_xticks(abpoa_list)
    ax.set_xticklabels(abpoa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('log')
    ax.set_xlabel('Block Size')
    ax.set_ylim(5, 600000)
    # ax.legend(['mt20x deMEM-abPOA', 'mt20x deMEM-abPOA-H', 'mt1x deMEM-abPOA', 'mt1x deMEM-abPOA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'f', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    '''
    Space
    '''
    # MAFFT Space
    d1 = mt_table_to_list(mafft_low, block_list, Memory, mt1x, getMemory)
    d2 = mt_table_to_list(mafft_low, block_list, Memory, mt20x, getMemory)
    d3 = mt_table_to_list(mafft_high, block_list, Memory, mt1x, getMemory)
    d4 = mt_table_to_list(mafft_high, block_list, Memory, mt20x, getMemory)
    p1 = [d1[i] for i in fftns_keys]
    p2 = [d2[i] for i in fftns_keys]
    p3 = [d3[i] for i in fftns_keys]
    p4 = [d4[i] for i in fftns_keys]
    ax = plt.subplot(337)
    ax.errorbar(fftns_list, p2, color = colors[0], fmt = '-o')
    ax.errorbar(fftns_list, p4, color = colors[1], fmt = '-o')
    ax.errorbar(fftns_list, p1, color = colors[2], linestyle = '--', fmt = '-o')
    ax.errorbar(fftns_list, p3, color = colors[3], linestyle = '--', fmt = '-o')
    ax.set_xticks(fftns_list)
    ax.set_xticklabels(fftns_keys, rotation = 45, ha = 'right')
    ax.set_xlabel('Block Size')
    ax.set_ylabel('Memory/MB')
    ax.set_yscale('log')
    ax.set_ylim(50, 140000)
    # ax.legend(['mt20x deMEM-FFTNS1', 'mt20x deMEM-FFTNS1-H', 'mt1x deMEM-FFTNS1', 'mt1x deMEM-FFTNS1-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    # ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'g', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # WMSA Space
    d1 = mt_table_to_list(wmsa_low, block_list, Memory, mt1x, getMemory)
    d2 = mt_table_to_list(wmsa_low, block_list, Memory, mt20x, getMemory)
    d3 = mt_table_to_list(wmsa_high, block_list, Memory, mt1x, getMemory)
    d4 = mt_table_to_list(wmsa_high, block_list, Memory, mt20x, getMemory)
    wmsa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    wmsa_list = [i for i in range(len(wmsa_keys))]
    p1 = [d1[i] for i in wmsa_keys]
    p2 = [d2[i] for i in wmsa_keys]
    p3 = [d3[i] for i in wmsa_keys]
    p4 = [d4[i] for i in wmsa_keys]
    ax = plt.subplot(338)
    ax.errorbar(wmsa_list, p2, color = colors[4], fmt = '-o')
    ax.errorbar(wmsa_list, p4, color = colors[5], fmt = '-o')
    ax.errorbar(wmsa_list, p1, color = colors[6], fmt = '-o', linestyle = '--')
    ax.errorbar(wmsa_list, p3, color = colors[7], fmt = '-o', linestyle = '--')
    ax.set_xticks(wmsa_list)
    ax.set_xticklabels(wmsa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('log')
    ax.set_xlabel('Block Size')
    ax.set_ylim(50, 140000)
    # ax.legend(['mt20x deMEM-WMSA', 'mt20x deMEM-WMSA-H', 'mt1x deMEM-WMSA', 'mt1x deMEM-WMSA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'h', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)

    # abPOA Space
    d1 = mt_table_to_list(abpoa_low, block_list, Memory, mt1x, getMemory)
    d2 = mt_table_to_list(abpoa_low, block_list, Memory, mt20x, getMemory)
    d3 = mt_table_to_list(abpoa_high, block_list, Memory, mt1x, getMemory)
    d4 = mt_table_to_list(abpoa_high, block_list, Memory, mt20x, getMemory)
    abpoa_keys = sorted(list(set(d1.keys()) & set(d2.keys()) & set(d3.keys()) & set(d4.keys())))
    abpoa_list = [i for i in range(len(abpoa_keys))]
    p1 = [d1[i] for i in abpoa_keys]
    p2 = [d2[i] for i in abpoa_keys]
    p3 = [d3[i] for i in abpoa_keys]
    p4 = [d4[i] for i in abpoa_keys]
    ax = plt.subplot(339)
    ax.errorbar(abpoa_list, p2, color = colors[8], fmt = '-o')
    ax.errorbar(abpoa_list, p4, color = colors[9], fmt = '-o')
    ax.errorbar(abpoa_list, p1, color = colors[10], fmt = '-o', linestyle = '--')
    ax.errorbar(abpoa_list, p3, color = colors[11], fmt = '-o', linestyle = '--')
    ax.set_xticks(abpoa_list)
    ax.set_xticklabels(abpoa_keys, rotation = 45, ha = 'right')
    ax.set_yscale('log')
    ax.set_xlabel('Block Size')
    ax.set_ylim(50, 140000)
    # ax.legend(['mt20x deMEM-abPOA', 'mt20x deMEM-abPOA-H', 'mt1x deMEM-abPOA', 'mt1x deMEM-abPOA-H'])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.text(-0.05, 1.05, 'i', transform = ax.transAxes, fontsize = 18, fontweight = 'bold', ha = 'center', va = 'center')
    # ax.grid(True)


    # plt.show()

if __name__ == '__main__':
    df = pd.read_csv('data.csv', encoding = 'gb2312')
    matplotlib_config(df)
    df = df[(df['Test_data'] == 'mt1x') | (df['Test_data'] == 'mt20x')]
    df["Maxmemory"] = df["Maxmemory"] / 1024
    del df['Q'], df['Q_view'], df['TC'], df['TC_view'], df['SP'], df['Sub_BlockSize'], df['SAMEMAlign_Version'], df['Vaild_tests'], df['Q_raw'], df['TC_raw'], df['Time_raw'], df['Memory_raw'], df['Time_view']
    pd.set_option('display.max_rows', 10)
    mt_graph(df)
    graph_final()
