#!/usr/bin/env python3

"""
Created on Tue Mar  5 08:09:15 2019

@author: Benoit
"""

#%% library loading
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy


#%% args
def usage():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--true-pos",)
    parser.add_argument("-f", "--false-pos",)
    parser.add_argument("-o", "--output-dir",)
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()

args = usage()


#%% define paths
# ~ cwd = os.path.dirname(os.path.realpath(__file__))
# ~ os.chdir(cwd)

trueP_file = args.true_pos
falseP_file = args.false_pos
plot_dir = args.output_dir


#%% seaborn global parameters and load vaf file as dataframe
sns.set(font_scale=1.4)
sns.set_style("white")                              # default theme


df_tp = pd.read_csv(trueP_file, header=0, sep='\t')


df_tp.info()

#%% vaf dot plot
x = "rdeer_VAF"
y = "depmap_VAF"

ax_tp = sns.scatterplot(data=df_tp, x=x, y=y, hue="gene", legend=False)
ax_tp.set_xlabel("Reindeer VAF")
ax_tp.set_ylabel("Depmap VAF")


### add Pearson/Spearman
pearson = df_tp[x].corr(df_tp[y])
spearman_rho, spearman_pval = scipy.stats.spearmanr(df_tp[x], df_tp[y])
r, p = scipy.stats.pearsonr(df_tp[x], df_tp[y])

t = plt.text(.05, .8,
         f"Pearson: {round(pearson, 2)}\nSpearman: {round(spearman_rho, 2)}",
         transform=ax_tp.transAxes,
         # backgroundcolor='0.90',
         )
t.set_bbox(dict(facecolor='white', alpha=0.6, linewidth=0))


plot_file = os.path.join(plot_dir, "plot-vaf-scatter.png")
fig = ax_tp.get_figure()
fig.savefig(plot_file, dpi=300, bbox_inches="tight")
plt.close()


#%% test de Wilcoxon

res = scipy.stats.wilcoxon(df_tp[x], df_tp[y])
print(res)
res = scipy.stats.ranksums(df_tp[x], df_tp[y])
print(res)
res = scipy.stats.mannwhitneyu(df_tp[x], df_tp[y], alternative='two-sided')
print(res)


#%% VAF boxplot

#### false positive file as dataframe
df_fp = pd.read_csv(falseP_file, header=0, sep='\t')

### add a status columns for both dataframe
df_tp["status"] = "true+"
df_fp["status"] = "false+"

### concatene false and true dataframe
df = pd.concat([df_tp, df_fp])

### build plot
x = "status"
y = "rdeer_VAF"
ax = sns.boxplot(data=df, x=x, y=y, hue=x, palette="husl")
ax.set_xlabel(None)
ax.set_ylabel("Reindeer VAF")

plot_file = os.path.join(plot_dir, "plot-vaf-box.png")
fig = ax.get_figure()
fig.savefig(plot_file, dpi=300, bbox_inches="tight")
plt.close()

#%% VAF violin
ax = sns.violinplot(data=df, x=x, y=y, hue=x, palette="husl")
ax.set_xlabel(None)
ax.set_ylabel("Reindeer VAF")

### add Wilcoxon test as Text
tp = df_tp['rdeer_VAF']     # VAF list if true+
fp = df_fp['rdeer_VAF']     # VAF list if false+
wilco = scipy.stats.mannwhitneyu(tp, fp)

t = plt.text(-.45, -20,
         f"Wilcoxon p-val: {wilco.pvalue:.3g}",
         )

plot_file = os.path.join(plot_dir, "plot-vaf-violin.png")
fig = ax.get_figure()
fig.savefig(plot_file, dpi=300, bbox_inches="tight")


#%% test de Wilcoxon

tp = df_tp['rdeer_VAF']
fp = df_fp['rdeer_VAF']
print(f"#true+: {len(tp)}")
print(f"#false+: {len(fp)}")
# ~ scipy.stats.mannwhitneyu(t
# ~ scipy.stats.ranksums(tp, fp)
