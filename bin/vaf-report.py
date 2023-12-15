#!/usr/bin/env python3

"""
Module doc
"""

import sys
import os
import argparse
import statistics
import pandas as pd


__appname__   = "vaf-report"
__shortdesc__ = "Short description."
__licence__   = "none"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    args = usage()
    report = ""                                 # report string, Markdown formatted
    file_info = get_files(args)                 # [ {file_path: xyz, gene: xyz} ]
    report += get_stats(args, file_info)

    output(args, report)


def get_files(args):
    '''
    get input files path and gene names
    '''
    file_info = []

    basedir, dirs, files = list(os.walk(args.vaf_dir))[0]

    for file in sorted(files):
        path = os.path.join(basedir, file)
        gene = file.split('_')[-1].rstrip(".tsv")
        file_info.append({"path":path, "gene": gene})

    return file_info


def get_stats(args, file_info):
    '''
    SOME STATISTICS
    '''

    stats = []      # [gene, TP, FP, FN, TN, presi, reca, speci, rdee_no_ref, depm_no_ref, rdee_depm_no_ref]

    for file in file_info:
        ### file as Dataframe
        df = pd.read_csv(file["path"], sep="\t")
        ### True+, False+, False-, true-
        TP = len(df[(df["rdeer_VAF"] > 0) & (df["depmap_VAF"] != "_")])
        FP = len(df[(df["rdeer_VAF"] > 0) & (df["depmap_VAF"] == "_")])
        FN = len(df[(df["rdeer_VAF"] == 0) & (df["depmap_VAF"] != "_")])
        TN = len(df[(df["rdeer_VAF"] == 0) & (df["depmap_VAF"] == "_")])
        ### Precision, Recall and Specificity
        # ~ print(file["gene"], len(TP), len(TP), len(FP))
        try:
            precision = round(TP / (TP + FP),2)         # precision
            recall = round(TP / (TP + FN),2)            # recall/Sensitivity
            specificity = round(TN / (TN + FP),2)       # specificity
        except ZeroDivisionError:
            precision = recall = specificity = 0
        ### Some other stuff
        rdee_no_ref = df["rdeer_VAF"].eq(100).sum()
        depm_no_ref = df["depmap_VAF"].eq("100").sum()
        rdee_depm_no_ref = len(df[ ( df["rdeer_VAF"] == 100 ) & ( df["depmap_VAF"] == "100" )  ])

        ### Add in dict and
        item = [ file["gene"],
                TP, FP, FN, TN,
                precision, recall, specificity,
                rdee_no_ref, depm_no_ref, rdee_depm_no_ref,
        ]
        stats.append(item)

    ### Compute Totals
    # ~ print([round(statistics.mean(i),2) for i in list(zip(*stats))[1:]])
    # ~ print([round(sum(i),2) for i in list(zip(*stats))[1:]])


    ### Reporting
    report = ""

    ### reporting stats as a table

    table  = "| gene | true+ | false+ | false- | true- | Precision | Recall | specificity | rdee_no_ref | depm_no_ref| rdee_depm_no_ref |\n"
    table += "|-----:|:-----:|:------:|:------:|:-----:|:---------:|:------:|:-----------:|:-----------:|:----------:|:----------------:|\n"
    for item in stats:
        table += f"| {'|'.join([str(i) for i in item])} |\n"
    report += f"\n{table}\n"

    report += "\n## Stats by gene:\n"

    report += "\n### Some explanations:\n"
    report += "\n**To compute the stats, we assume than:**\n\n"
    report += "- found in reindeer     : rdeer_VAF > 0\n"
    report += "- not found in reindeer : rdeer_VAF == 0\n"
    report += "- found in depmap       : depmap_VAF != '_'\n"
    report += "- not found in depmap   : depmap_VAF == '_'\n"
    report += "\n**The 3 last columns:**\n\n"
    report += "- rdee_no_ref: Reindeer abundance = 0 in REF, but positive in ALT\n"
    report += "- depm_no_ref: Depmap abundance = 0 in REF, but positive in ALT\n"
    report += "- rdee_depm_no_ref: in both Reindeer and Depmap, abundance = 0 in REF\n"


    ### write table as tsv
    outfile = f"{os.path.splitext(args.output.name)[0]}.tsv"
    header = "gene\ttrue+\tfalse+\tfalse-\ttrue-\tPrecision\tRecall\tspecificity\trdee_no_ref\tdepm_no_ref\trdee_depm_no_ref\n"
    with open(outfile, 'w') as fh:
        fh.write(header)
        for items in stats:
            fh.write('\t'.join([str(i) for i in items]) + '\n')

    return report



def output(args, report):
    """ Function doc """
    args.output.write(report)




def usage():
    parser = argparse.ArgumentParser()
    parser.add_argument("vaf_dir",
                        help="Input VAF dir",
                       )
    parser.add_argument("-o", "--output",
                        type=argparse.FileType('w'),
                        help="Ouput markdown report file",
                        required=True,
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()


#!/usr/bin/env python3

"""
Stats and plots from vaf
"""







'''
sns.set_style("white")             # theme par défaut
'''





''' #%% PLOT VAF

TSV_FILE = os.path.join(TSV_DIR, "CCLE_mutations_DepMap_22Q2_BRAF.tsv")

with open(TSV_FILE) as fh:
    header = fh.readline()
    df_vaf = pd.DataFrame(columns=['Location', 'VAF', 'Source'])
    for row in fh:
        gene, mut, srr, r_af, r_vaf, d_af, d_vaf, *rest = row.split('\t')
        ### VAF
        if int(r_vaf) > 0 and d_vaf != '_':
            df_vaf.loc[df_vaf.index.size] = (f"{mut[:16]}:{srr[3:]}", int(r_vaf), "Reindeer" )
            df_vaf.loc[df_vaf.index.size] = (f"{mut[:16]}:{srr[3:]}", int(d_vaf), "Depmap" )


# Draw a nested barplot by rdeer_VAF and depmap_VAF
ax = sns.catplot(data=df_vaf[:100],
                 kind="bar",
                 x="Location", y="VAF", hue="Source",
                 palette = ["#DB5F57", "#57D3DB",],
                 aspect=3,
                 errorbar=None,
                 )


plt.xticks(rotation=90)
plt.xlabel("mutation")

sns.move_legend(ax,
                bbox_to_anchor=(0.91,0.97),
                loc=1,     # légende dedans en haut à droite
                title = None,
                fontsize='13',
                frameon=True,
                )
# sns.set(rc={'figure.figsize':(15,4)}, style="whitegrid")
# ax.set_axis_labels("depmap VAF", "reindeer VAF")
# ax.legend.set_title("VAF Comparison")

### On enregistre dans un fichier
outfile = "plot-vaf-af.png"
ax.figure.savefig(os.path.join(outfile), dpi=300, bbox_inches='tight')  # bbox_inches pour éviter de perdre le bas de la figure
'''




''' #%% PLOT AF

TSV_FILE = os.path.join(TSV_DIR, "CCLE_mutations_DepMap_22Q2_BRAF.tsv")

with open(TSV_FILE) as fh:
    header = fh.readline()
    df_af = pd.DataFrame(columns=['Location', 'AF', 'Source'])
    for row in fh:
        gene, mut, srr, r_af, r_vaf, d_af, d_vaf, *rest = row.split('\t')
        ### VAF
        if int(r_vaf) > 0 and d_vaf != '_':
            df_af.loc[df_af.index.size] = (f"{mut[:16]}:{srr[3:]}", int(r_af.split(':')[0]), "Reindeer" )
            df_af.loc[df_af.index.size] = (f"{mut[:16]}:{srr[3:]}", int(d_af.split(':')[0]), "Depmap" )


# Draw a nested barplot by rdeer_VAF and depmap_VAF
ax = sns.catplot(data=df_af[:100],
                 kind="bar",
                 x="Location", y="AF", hue="Source",
                 palette = ["#DB5F57", "#57D3DB",],
                 aspect=3,
                 errorbar=None,
                 )


plt.xticks(rotation=90)
plt.xlabel("mutation")

sns.move_legend(ax,
                bbox_to_anchor=(0.91,0.97),
                loc=1,     # légende dedans en haut à droite
                title = None,
                fontsize='13',
                frameon=True,
                )


### On enregistre dans un fichier
outfile = "plot-af-af.png"
ax.figure.savefig(os.path.join(outfile), dpi=300, bbox_inches='tight')  # bbox_inches pour éviter de perdre le bas de la figure
'''



''' #%% PLOT True+, False+, False-

files = [ os.path.join(TSV_DIR, file) for file in next(os.walk(TSV_DIR))[2] if file.endswith(".tsv")]

df = pd.DataFrame(columns=['gene', 'number of pairs', 'type'])  # type: true_pos, false_pos, false_neg

for file in files:
    # print(os.path.basename(file))
    types = {
        "true_pos": 0,
        "false_pos": 0,
        "false_neg": 0,
    }
    with open(file) as fh:
        header = fh.readline()
        rows = fh.readlines()

        ### handle file with only header
        if not rows:
            print(file)
            gene = file.split('_')[-1].rstrip(".tsv")
            for k,val in types.items():
                df.loc[df.index.size] = (gene, val, k)
            continue

        ### build dataframe
        for row in rows:
            gene, mut, srr, r_af, r_vaf, d_af, d_vaf, *rest = row.split('\t')
            if int(r_vaf) > 0 and d_vaf != '_':
                types["true_pos"] += 1
            elif int(r_vaf) > 0 and d_vaf == '_':
                types["false_pos"] -= 1
            elif int(r_vaf) == 0 and d_vaf != '_':
                types["false_neg"] -=1
    for k,val in types.items():
        df.loc[df.index.size] = (gene, val, k)

df = df.sort_values("gene")
'''


''' #%%   Build plot wirh SEABORN

ax = sns.catplot(data=df,
                 kind="bar",
                 # stacked=True,
                 # multiple = 'stack',
                 x="gene", y="number of pairs", hue="type",
                 palette = ["#DB5F57", "#57D3DB", "#DBA857",],
                 aspect=2,
                 errorbar=None,
                )

plt.xticks(rotation=90)
plt.xlabel("Gene")

sns.move_legend(ax,
                bbox_to_anchor=(0.84,0.15),
                loc=4,     # légende dedans en haut à droite
                title = None,
                fontsize='13',
                frameon=True,
                )


### On enregistre dans un fichier
outfile = "plot-compare.png"
ax.figure.savefig(os.path.join(outfile), dpi=300, bbox_inches='tight')  # bbox_inches pour éviter de perdre le bas de la figure
'''



''' #%% Build plot with PLOTLY
### https://plotly.com/python/bar-charts/

for file in files:
    # print(os.path.basename(file))
    types = {
        "true_pos": 0,
        "false_pos": 0,
        "false_neg": 0,
    }
    data = {}
    with open(file) as fh:
        header = fh.readline()
        rows = fh.readlines()

        ### handle file with only header
        if not rows:
            print(file)
            gene = file.split('_')[-1].rstrip(".tsv")
            for k,val in types.items():
                df.loc[df.index.size] = (gene, val, k)
            continue

        ### filter rows
        for row in rows:
            gene, mut, srr, r_af, r_vaf, d_af, d_vaf, *rest = row.split('\t')

            ### filter on min abundance
            if int(r_vaf) > 0 and d_vaf != '_':


            if int(r_vaf) > 0 and d_vaf != '_':
                types["true_pos"] += 1
            elif int(r_vaf) > 0 and d_vaf == '_':
                types["false_pos"] -= 1
            elif int(r_vaf) == 0 and d_vaf != '_':
                types["false_neg"] -=1

color_discrete_map = {"false_neg": "#DB5F57", "true_pos": "#57D3DB", "false_pos":"#DBA857"}

fig = px.bar(
    df,
    x="gene",
    y="number of pairs",
    color="type",
    template= "plotly_white",
    # color_discrete_map = color_discrete_map,
    )

fig.update_layout(legend_title=None)
fig.show(renderer="svg")
fig.write_image("true_and_false.svg")
'''






#%% test zone









''' #%% COMMON

### Prepare output
os.makedirs(OUTDIR, exist_ok=True)
### Load file as dataframe
homedir = os.path.expanduser("~")
### open ssh connection on ella
if not os.listdir(f"{homedir}/sshfs/ella"):
    os.system('sshfs ella:/data ~/sshfs/ella')

pwd = f"{homedir}/sshfs/ella/nas/projects/kmer-collections/mutations/DepMap"
table = os.path.join(pwd, DEPMAP_RELEASE, OUTPUT ,"table.tsv")
table_af = os.path.join(pwd, DEPMAP_RELEASE, OUTPUT_WITH_AF ,"table.tsv")


os.chdir(pwd)
df = pd.read_csv(table, sep="\t")
df_af = pd.read_csv(table_af, sep="\t")

### avoid kmers without expression in depmap or too different (e.g. TP53)
df = df[(df['Kmer-alt>0'] > 0) & (df['Kmer-alt>0'] < 250) ]
df_af = df_af[(df_af['Kmer-alt>0'] > 0) & (df_af['Kmer-alt>0'] < 250) ]

### Renommage des colonnes d'intérêt (pour la légende)
df = df.rename(columns = {
    "Total SRR": "DepMap",
    "Kmer-alt>0 (ok)": "Kmer-alt",
    "Kmer-alt>6 (hit<100 ok)": "Refined kmer-alt",
})
df_af = df_af.rename(columns = {
    "Total SRR": "DepMap",
    "Kmer-alt>0 (ok)": "Kmer-alt",
    "Kmer-alt>6 (hit<100 ok)": "Refined kmer-alt",
})

### 2 plots at each time, with and without Allelic Frequency
dfs = {
       'df-af': df_af,
       'df': df,
       }

### print columns and rows
lcols = list(df.columns)
'''


''' #%% plot d1 & plot d2_af : pd.melt() --- sensibility


g = "Gene"
x = "DepMap"                     # Total SRR
y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)


for name, df in dfs.items():

    var_name = "Condition"
    # palette = ["royalblue", "dimgray", "tomato", ]
    palette = ["#DB5F57", "#57D3DB", ]
    markers = ["s", "o", "D"]
    tidy = df[[g,x,y]].melt(id_vars="Gene", var_name=var_name)
    # ATTENTION de ne pas réutilser un nom de colonne existant, j'ai ajouté un espace à "Gene "

    pl = sns.relplot(
            data=tidy,              # dataframe
            x=g, y="value",         # les x et les y
            hue=var_name,           # couleurs en fonction de la variable
            style=var_name,         # formes en fonction de la variable
            palette=palette,        # couleur des items
            markers=markers,        # redéfinition des formes
            edgecolor="dimgrey",    # couleur des bordures de point
            s = 100,                # taille des points
            aspect=2,               # ratio largeur/hauteur
            )
    pl.set(ylabel="mut/sample")

    sns.move_legend(pl,
                    bbox_to_anchor=(0.81,0.97), loc=9,     # légende dedans en haut à droite
                    # title=var_name, title_fontsize='14',   # titre et taille du titre de la légende
                    title = None,
                    fontsize='13',                         # taille des fontes
                    markerscale=1.5,                       # agrandir les markers
                    frameon=True,
                    )

    pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés

    ### On enregistre dans un fichier
    outfile = f"plot-sensibility-{name}.png"
    pl.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)
'''





''' #%% plot d1 & plot d2_af : pd.melt() --- specificity


g = "Gene"
x = "DepMap"                     # Total SRR
y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)
z = "Refined kmer-alt"           # specificity   # (P+R) Kmer-alt>6 (hit<100 ok)


for name, df in dfs.items():

    var_name = "Condition"
    palette = ["#DB5F57", "#57D3DB", "#DBA857"]
    markers = ["s", "o", "D"]
    tidy = df[[g,x,y,z]].melt(id_vars="Gene", var_name=var_name)
    # ATTENTION de ne pas réutilser un nom de colonne existant, j'ai ajouté un espace à "Gene "

    pl = sns.relplot(
            data=tidy,              # dataframe
            x=g, y="value",         # les x et les y
            hue=var_name,           # couleurs en fonction de la variable
            style=var_name,         # formes en fonction de la variable
            palette=palette,        # couleur des items
            edgecolor="dimgrey",    # couleur des bordures de point
            markers=markers,        # redéfinition des formes
            s = 100,                # taille des points
            aspect=2,               # ratio largeur/hauteur
    )
    pl.set(ylabel="mut/sample")


    sns.move_legend(pl,
                    bbox_to_anchor=(0.74,0.97), loc=9,     # légende dedans en haut à droite
                    # title=var_name, title_fontsize='14',   # titre et taille du titre de la légende
                    title = None,
                    fontsize='13',                         # taille des fontes
                    markerscale=1.5,                       # agrandir les markers
                    frameon=True,
    )

    pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés

    ### On enregistre dans un fichier
    outfile = f"plot-specificity-{name}.png"
    pl.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)
'''



''' #%%  barplot sans seaborn
### SANS SEABORN

# https://pythonguides.com/stacked-bar-chart-matplotlib/

marked_by = "Gene"
x = "DepMap"                     # Total SRR
y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)
z = "Refined kmer-alt"           # specificity   # (P+R) Kmer-alt>6 (hit<100 ok)

# palette = ["#DB5F57", "#57D3DB", "#DBA857"]
palette = ["#DBA857", "#57D3DB", "#DB5F57"]


for name, df in dfs.items():
    tidy = df[[marked_by,x,y,z]]


    total = df[x] + df[y] + df[z]
    rel = tidy[tidy.columns[1:]].div(total, 0)*100  # valeurs rapportées en pourcentage
    coeff = tidy[x] / 100                           # coefficents pour DepMap
    new = tidy[tidy.columns[1:]].div(coeff, 0)      # les autres colonnes sont divisées par ce coeff
    new[x] = new[x] - new[y]                        # Depmap = DepMap - Kmer-alt
    new[y] = new[y] - new[z]                        # kmer-alt = kmer-alt - Refined
    rel = new[[z, y, x]]                            # on inverse l'ordre


    rel[marked_by] = df[marked_by]

    # fig, ax = plt.subplots()

    ax = rel.plot (
        x = marked_by,                      # sur les X
        kind = 'bar',                       # graphique en barre
        stacked = True,                     # stackée
        # title = f"Percentage - {name}",
        mark_right = True,                  # je ne sais pas, ça ne change rien pour l'instant
        color = palette,
        )

    ax.legend(reverse=True, loc=4)          # ordre inversé, légende en bas à droite

    ### Pour mettre les valeurs de pourcentage (mais le 2eme plot plante)
    # for bar in ax.patches:
        # height = bar.get_height()
        # width = bar.get_width()
        # x = bar.get_x()
        # y = bar.get_y()
        # label_text = round(float(height))
        # label_x = x + width / 2
        # label_y = y + height / 2
        # ax.text(label_x, label_y, label_text, ha='center', va='center') # , rotation=90)

    plt.xticks(rotation=45, horizontalalignment='right', fontsize=8)    # noms des gènes penchés
    plt.yticks(fontsize=9)

    ### taille de la figure
    fig = plt.gcf()
    fig.set_size_inches(9, 4.5)

    ### On enregistre dans un fichier
    outfile = f"plot-stackbar-percent-{name}.png"
    ax.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)
'''




''' #%% Comparo couples ALT and REF

alt_file = os.path.join(pwd, "output/reindeer-merged")
print(os.listdir())
'''

###############################################################################
#                                OLDER
###############################################################################


''' #%% relation entre les kmers de reindeer et ceux de DepMap

### pour mettre les nomes de gene : https://stackoverflow.com/questions/58714451/add-text-to-each-point-in-a-seaborn-relplot

def lmplot(g, x, y, suffix):
    for name, data in dfs.items():
        pl = sns.lmplot(data=data, x=x, y=y, aspect=1.3)
        ax = pl.axes[0,0]
        for k,v in data[[g,x,y]].iterrows():
            ax.text(v[1], v[2] + 3, v[0], horizontalalignment='center', fontsize=5)
        outfile= f"plot-{suffix}_{name}.png"
        pl.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)

g = "Gene"

x = "Kmer-alt>6 (hit<100)"
y = "Kmer-alt>6 (hit<100 ok)"
suffix = 'a_pos'
lmplot(g, x, y, suffix)

x = "(P+R) Kmer-alt>6 (hit<100)"
y = "(P+R) Kmer-alt>6 (hit<100 ok)"
suffix = 'a_rescued'
lmplot(g, x, y, suffix)
'''


''' #%% plot 2.af : Bon avec pd.melt()

### voir https://stackoverflow.com/questions/40877135/plotting-two-columns-of-dataframe-in-seaborn
### pour pd.melt() voir https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.melt.html


def plot_2(g, x, y, suffix):
    for name, df in dfs.items() :

        var_name = "Condition"
        markers = ["s", "o"]
        tidy = df[[g,x,y]].melt(id_vars="Gene", var_name=var_name)
        # ATTENTION de ne pas réutilser un nom de colonne existant, j'ai ajouté un espace à "Gene "

        pl = sns.relplot(
                data=tidy,              # dataframe
                x=g, y="value",         # les x et les y
                hue=var_name,           # couleurs en fonction de la variable
                style=var_name,         # formes en fonction de la variable
                markers=markers,        # redéfinition des formes
                s = 80,                 # taille des points
                aspect=2,               # ratio largeur/hauteur
                )
        sns.move_legend(pl, bbox_to_anchor=(0.69,0.95), loc=9 )         # légende dedans en haut à droite

        pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés
        ### On enregistre dans un fichier
        outfile= f"plot-{suffix}_{name}.png"
        pl.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)

x = "Kmer-alt>6 (hit<100)"
y = "Kmer-alt>6 (hit<100 ok)"
suffix = 'b_pos'
plot_2(g, x, y, suffix)

x = "(P+R) Kmer-alt>6 (hit<100)"
y = "(P+R) Kmer-alt>6 (hit<100 ok)"
suffix = 'b_rescued'
plot_2(g, x, y, suffix)
'''


''' #%% plot a2 & plot a2_af : relation entre les kmers de reindeer et ceux de DepMap
### pour mettre les nomes de gene : https://stackoverflow.com/questions/58714451/add-text-to-each-point-in-a-seaborn-relplot
g = "Gene"
x = "Kmer-alt>6 pairs (hit<100)"    # changed
y = "Kmer-alt>6 (hit<100 ok)"

pl = sns.lmplot(data=df, x=x, y=y, aspect=1.3)
ax = pl.axes[0,0]
for k,v in df[[g,x,y]].iterrows():
    ax.text(v[1], v[2] + 3, v[0], horizontalalignment='center', fontsize=5)
pl.figure.savefig('plot_a2.png', dpi=300)


pl_af = sns.lmplot(data=df_af, x=x, y=y, aspect=1.3)
ax_af = pl_af.axes[0,0]
for k,v in df_af[[g,x,y]].iterrows():
    ax_af.text(v[1], v[2] + 3, v[0], horizontalalignment='center', fontsize=5)
pl_af.figure.savefig('plot_a2_af.png', dpi=300)
'''





''' #%% plot b2 & plot b2_af : pd.melt()
### voir https://stackoverflow.com/questions/40877135/plotting-two-columns-of-dataframe-in-seaborn
### pour pd.melt() voir https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.melt.html



var_name = "Condition"
markers = ["s", "o"]
tidy = df[[g,x,y]].melt(id_vars="Gene", var_name=var_name)
# ATTENTION de ne pas réutilser un nom de colonne existant, j'ai ajouté un espace à "Gene "

pl = sns.relplot(
        data=tidy,              # dataframe
        x=g, y="value",         # les x et les y
        hue=var_name,           # couleurs en fonction de la variable
        style=var_name,         # formes en fonction de la variable
        markers=markers,        # redéfinition des formes
        s = 80,                 # taille des points
        aspect=2,               # ratio largeur/hauteur
        )
sns.move_legend(pl, bbox_to_anchor=(0.69,0.95), loc=9 )         # légende dedans en haut à droite

pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés
pl.figure.savefig('plot_b2.png', dpi=300)                        # On enregistre



tidy_af = df_af[[g,x,y]].melt(id_vars="Gene", var_name=var_name)
pl_af = sns.relplot(
        data=tidy_af,           # dataframe
        x=g, y="value",         # les x et les y
        hue=var_name,           # couleurs en fonction de la variable
        style=var_name,         # formes en fonction de la variable
        markers=markers,        # redéfinition des formes
        s = 80,                 # taille des points
        aspect=2,               # ratio largeur/hauteur
        )
sns.move_legend(pl_af, bbox_to_anchor=(0.69,0.95), loc=9 )         # légende dedans en haut à droite

pl_af.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés
pl_af.figure.savefig('plot_b2_af.png', dpi=300)
'''


''' #%% plot c1 & plot c2_af : pd.melt()

g = "Gene"
x = "Total SRR"
y = "Kmer-alt>0 (ok)"
z = "(P+R) Kmer-alt>0 (ok)"

for name, df in dfs.items():

    var_name = "Condition"
    markers = ["s", "o", "D"]
    tidy = df[[g,x,y,z]].melt(id_vars="Gene", var_name=var_name)
    # ATTENTION de ne pas réutilser un nom de colonne existant, j'ai ajouté un espace à "Gene "

    pl = sns.relplot(
            data=tidy,              # dataframe
            x=g, y="value",         # les x et les y
            hue=var_name,           # couleurs en fonction de la variable
            style=var_name,         # formes en fonction de la variable
            palette=["royalblue", "dimgray", "tomato", ],
            markers=markers,        # redéfinition des formes
            s = 80,                 # taille des points
            aspect=2,               # ratio largeur/hauteur
            )
    sns.move_legend(pl, bbox_to_anchor=(0.69,0.95), loc=9 )         # légende dedans en haut à droite

    pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés
    pl.figure.savefig('plot_c1.png', dpi=300)                        # On enregistre

    ### On enregistre dans un fichier
    outfile = f"plot-c_{name}.png"
    pl.figure.savefig(os.path.join(OUTDIR, outfile), dpi=300)
'''
