#!/usr/bin/env python3


"""
Plots and stats from input table
"""

import sys
import os
import argparse
import yaml
# import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    """ Function doc """
    args = usage()
    ### open config file and configure some path
    with open("config.yaml") as fh:
        config = yaml.full_load(fh)
    ### build dataframe
    df = common(args, config)
    ### build plots
    sensibility(args, df)
    # ~ if config["filter"]["enabled"]:
    specificity(args, df)
    barplot_percent(args, df)


def common(args, config):
    sns.set_style("white")             # défault theme

    df = pd.read_csv(args.table, sep="\t")

    ### remove last row if gene == Total
    if df.loc[len(df)-1, "Gene"] == "Total":
        df.drop(len(df)-1)

    ### avoid kmers without expression in depmap or too different (e.g. TP53)
    df = df[(df['Kmer-alt>0'] > 0) & (df['Kmer-alt>0'] < 250) ]

    ### Renommage des colonnes d'intérêt (pour la légende)
    a_min = args.abund_min if config["filter"]["enabled"] else 1
    r_max = args.recur_max if config["filter"]["enabled"] else "100%"
    df = df.rename(columns = {
        "Total SRR": "DepMap",
        "Kmer-alt>0 (ok)": "Kmer-alt",
        f"Kmer-alt>{a_min} (hit<{r_max} ok)": "Refined kmer-alt",
    })

    return df


def sensibility(args, df):

    g = "Gene"
    x = "DepMap"                     # Total SRR
    y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)

    var_name = "Condition"
    palette = ["#DB5F57", "#57D3DB"]
    markers = ["s", "o"]
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
                # title=var_name, title_fontsize='14', # titre et taille du titre de la légende
                title = None,
                fontsize='13',                         # taille des fontes
                markerscale=1.5,                       # agrandir les markers
                frameon=True,
                )

    pl.set_xticklabels(rotation=45, horizontalalignment='right')    # noms des gènes penchés

    ### On enregistre dans un fichier
    outfile = f"plot-sensibility-{args.version}.png"
    pl.figure.savefig(os.path.join(args.output_dir, outfile), dpi=300)


def specificity(args, df):

    g = "Gene"
    x = "DepMap"                     # Total SRR
    y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)
    z = "Refined kmer-alt"           # specificity   # (P+R) Kmer-alt>6 (hit<100 ok)


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
    outfile = f"plot-specificity-{args.version}.png"
    pl.figure.savefig(os.path.join(args.output_dir, outfile), dpi=300)


def barplot_percent(args, df):
    ### SANS SEABORN
    # https://pythonguides.com/stacked-bar-chart-matplotlib/

    marked_by = "Gene"
    x = "DepMap"                     # Total SRR
    y = "Kmer-alt"                   # sensibility   # (P+R) Kmer-alt>0 (ok)
    z = "Refined kmer-alt"           # specificity   # (P+R) Kmer-alt>6 (hit<100 ok)

    # palette = ["#DB5F57", "#57D3DB", "#DBA857"]
    palette = ["#DBA857", "#57D3DB", "#DB5F57"]

    tidy = df[[marked_by,x,y,z]]

    total = df[x] + df[y] + df[z]
    rel = tidy[tidy.columns[1:]].div(total, 0)*100  # valeurs rapportées en pourcentage
    coeff = tidy[x] / 100                           # coefficents pour DepMap
    new = tidy[tidy.columns[1:]].div(coeff, 0)      # les autres colonnes sont divisées par ce coeff
    new[x] = new[x] - new[y]                        # Depmap = DepMap - Kmer-alt
    new[y] = new[y] - new[z]                        # kmer-alt = kmer-alt - Refined
    rel = new[[z, y, x]]                            # on inverse l'ordre

    rel[marked_by] = df[marked_by]

    ax = rel.plot (
        x = marked_by,                      # sur les X
        kind = 'bar',                       # graphique en barre
        stacked = True,                     # stackée
        # title = f"Percentage - {name}",
        mark_right = True,                  # je ne sais pas, ça ne change rien pour l'instant
        color = palette,
        )

    ax.legend(reverse=True, loc=4)          # ordre inversé, légende en bas à droite

    plt.xticks(rotation=45, horizontalalignment='right', fontsize=8)    # noms des gènes penchés
    plt.yticks(fontsize=9)

    ### taille de la figure
    fig = plt.gcf()
    fig.set_size_inches(9, 4.5)

    ### On enregistre dans un fichier
    outfile = f"plot-stackbar-percent-{args.version}.png"
    ax.figure.savefig(os.path.join(args.output_dir, outfile), dpi=300)


def usage():
    """
    Help function with argument parser.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("table",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-r", "--recur-max",
                        default=100,
                       )
    parser.add_argument("-a", "--abund-min",
                        default=6,
                       )
    parser.add_argument("-o", "--output_dir",
                        required=True,
                       )
    parser.add_argument("-v", "--version",
                        required=True,
                       )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()



if __name__ == "__main__":
    main()
