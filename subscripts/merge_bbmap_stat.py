#!/bin/python3
# %%
import pandas as pd
import glob
import argparse
import re
import os

# %%
parser = argparse.ArgumentParser(
    prog='bbmap statistics merge',
    description='Merges bbmap statistics files and evaluates mapping percentage',
    epilog=' Made by glitch 	♡( ◡‿◡ )')
parser.add_argument('-i', '--input_reports', required=True, type=argparse.FileType('r'), nargs='+', help="Specify reports name with the glob expression")
parser.add_argument('-o', '--output_file_name', default='bbmap_statistics_merged', help="Specify the name of the output file")
args = parser.parse_args()
reports = args.input_reports
fname = args.output_file_name

# %%
df = pd.read_csv(reports[0], sep='\t', names=range(5)).iloc[[
    0, *range(8, 12), *range(22, 26)]].reset_index(drop=True)
df.iloc[0, 4] = re.search('\d+', df.iloc[0, 2])[0]
df.iloc[0, 2] = df.iloc[0, 1]
df = df.iloc[:, [0, 2, 4]].T.reset_index(drop=True).T
df = df.apply(pd.to_numeric, errors='coerce').fillna(df)
df_all = df.copy(deep=True)
print(f'[{1}/{len(reports)}] {reports[0].name} is processed\r')

# %%
for num, report in enumerate(reports[1:]):
    print(f'[{num+2}/{len(reports)}] {report.name} is processed\r')
    df = pd.read_csv(report, sep='\t', names=range(5)).iloc[[
        0, *range(8, 12), *range(22, 26)]].reset_index(drop=True)
    df.iloc[0, 0] = "Used:"
    df.iloc[0, 4] = re.search('\d+', df.iloc[0, 2])[0]
    df.iloc[0, 2] = df.iloc[0, 1]
    df = df.iloc[:, [0, 2, 4]].T.reset_index(drop=True).T
    df = df.apply(pd.to_numeric, errors='coerce').fillna(df)
    df_all.iloc[[0, *range(2, 5), *range(6, 9)], [1, 2]
                ] += df.iloc[[0, *range(2, 5), *range(6, 9)], [1, 2]]

# %%
total_mapped_reads = df_all.iloc[2,1]+df_all.iloc[7,1]
total_unmapped_reads = df_all.iloc[0,1] - total_mapped_reads
total_mapped_bases = df_all.iloc[2,2]+df_all.iloc[7,2]
total_unmapped_bases = df_all.iloc[0,2] - total_mapped_bases
df_all.iloc[0.5], df_all.iloc[0.6], df_all.iloc[0.7], df_all.iloc[0.8] = ['Total unmapped:', total_unmapped_reads, total_unmapped_bases],['Total mapped:', total_mapped_reads, total_mapped_bases],['Mapped/unmapped:', f'{round(total_mapped_reads/total_unmapped_reads*100, ndigits=5)}%', f'{round(total_mapped_bases/total_unmapped_bases, ndigits=5)*100}%'], ['Mapped/Used:', f'{round(total_mapped_reads/df_all.iloc[0,1]*100, ndigits=5)}%', f'{round(total_mapped_bases/df_all.iloc[0,2], ndigits=5)*100}%']
df_all = df_all.sort_index().reset_index(drop=True)

# %%
df_all = df_all.rename(columns={1:'Reads', 2:'Bases'})
df_all.iloc[0,0] = "Used:"

# %%
df_all.to_csv(f'{fname}.csv', index=False)
# view results with `column -s, -t < $output_file_name`
