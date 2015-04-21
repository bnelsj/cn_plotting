import argparse
import pandas as pd
import cluster
import genotyper as gt
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pybedtools
import pybedtools.contrib.plotting
import sys

def make_legend(sample_legend_dims, axis, pop_data, alpha = 0.6):
    plt.sca(axis)
    axis.set_axis_off()
    legend_entries = []

    legend_data = pop_data[["label", "color"]].drop_duplicates(subset = ["label"])

    for idx, pop_label in legend_data.iterrows():
        pop_entry = mpl.patches.Patch(label=pop_label["label"], color=pop_label["color"], alpha=alpha)
        legend_entries.append(pop_entry)

    plt.legend(handles=legend_entries, bbox_to_anchor=sample_legend_dims, borderaxespad=0.)


def make_region_annotation(axis, bedstring, plot_domain, fn_regions, merge=True, strand=True):
    plt.sca(axis)
    axis.set_axis_off()
    focal_region = pybedtools.BedTool(bedstring, from_string=True)
    full_region_set = pybedtools.BedTool(fn_regions)
    region_set = full_region_set.intersect(focal_region)
    if merge and len(region_set) > 0:
        region_set = region_set.merge(c = "4,5,6", o = "distinct,max,first") 

    plot_end = plot_domain[1] - plot_domain[0]
    axis.set_xlim(*plot_domain)  
    xmin = int(plot_domain[0])
    region_labels = map(lambda x: (x[1], x[3]), region_set)
    track = pybedtools.contrib.plotting.Track(region_set, alpha=1, visibility="dense", stranded=strand, facecolor="lightgray", ybase=0.25, yheight=0.5)
    for start, label in region_labels:
        strt = max(int(start), xmin)
        axis.text(strt, 0.75, label, ha="left", va="bottom", size=10, rotation=45)
    axis.add_collection(track)

    axis.axhline(color = "k", y = 0.5, zorder=-1)

    return region_set

def plot_violins(axis, region_list, g, pop_data, colors):
    import seaborn
    plt.sca(axis)
    seaborn.set(style="whitegrid")
    sp_list = pop_data["super_pop"].unique().tolist()
    super_pop = map(lambda x: pop_data.loc[pop_data.sample == x, "super_pop"].values[0], g.indivs)
    color = map(lambda x: pop_data.loc[pop_data.sample == x, "color"].values[0], g.indivs)
    sp_map = {}
    color_list = []

    for sp in sp_list:
        sp_map[sp] = np.array([i for i, x in enumerate(super_pop) if x == sp])
        color_list.extend([color[i] for i in sp_map[sp].tolist()])
    cp_map = {}
    name_list = []
    sbn_palette = seaborn.color_palette(map(lambda x: tuple(x), colors))
    sbn_palette = seaborn.color_palette("hls", len(colors))

    cp_df_list = []

    for i, interval in enumerate(region_list):
        pop_data_tmp = pop_data.copy()
        chr, start, end, name = str(interval).split()
        X, idx_s, idx_e = g.get_gt_matrix(chr, start, end)
        cps = np.mean(X, 1)
        pop_data_tmp["Copy number"] = pop_data_tmp.apply(lambda x: cps[g.indivs.index(x.sample)], axis=1)
        pop_data_tmp["Name"] = name
        cp_df_list.append(pop_data_tmp)

    pop_copy_num = pd.concat(cp_df_list)

    with seaborn.color_palette("hls", len(colors)):
        p1 = seaborn.FacetGrid(pop_copy_num, col="Name", hue=label, ylim=(0, math.ceil(pop_copy_num["Copy number"].max())))
        p1.map(seaborn.violinplot, x="super_pop", y="Copy number", data=pop_copy_num, inner=None, bw=0.2, lw=1., color="white", facecolors="none", edgecolor="black", alpha=1, cut=0)
        #p1.map(seaborn.stripplot, "super_pop", "Copy number", label, data=pop_copy_num, jitter=True, edgecolor=sbn_palette, facecolors="none", split=False, alpha=0.5)
        #p1 = seaborn.factorplot(x="super_pop", y="Copy number", data = pop_copy_num, col="Name", kind="violin", split=False, bw=0.2, lw=1.5, facecolors="none", edgecolor="black")
        p1.map(seaborn.stripplot, x="super_pop", y="Copy number", hue=label, data=pop_copy_num, jitter=True, edgecolor=sbn_palette, facecolors="none", split=False, alpha=0.7, size=1, lw=0.4)
        p1.set_xlabels("Population").set_titles("{col_name}")
        p1.add_legend()

def make_line_plot(data, region_name, pop_data, axis, pop_color_scheme, alpha=0.6, max_cp = None, sort_col = None):
    plt.sca(axis)
    plt.xlabel("Base pairs from %s" % data.columns[0])
    plt.ylabel("WSSD copy number")
    plt.title(region_name)

    starts = map(int, data.columns)
    line_plot_axis.set_xlim(0, starts[-1] - starts[0])
    labels = line_plot_axis.get_xticks().tolist()
    new_labels = [str(int(x)) for x in labels]
    line_plot_axis.set_xticklabels(new_labels)

    line_plot_axis.tick_params(direction='out')
    line_plot_axis.xaxis.set_ticks_position('bottom')
    line_plot_axis.yaxis.set_ticks_position('left')

    # Sort samples by pop so that super_pops aren't completely drawn over
    if sort_col is not None:
        pop_data = pop_data.sort(sort_col)

    # Plot D&N and Archaics on top for hgdp settings
    if pop_color_scheme == "hgdp":
        pop_data = pop_data.loc[~pop_data["super_pop"].isin(["D&N", "ARC"])].append(pop_data.loc[pop_data["super_pop"].isin(["ARC"])]).append(pop_data.loc[pop_data["super_pop"].isin(["D&N"])])

    for idx, pop_label in pop_data.iterrows():
        indiv = pop_label.sample
        plt.plot(map(lambda x: str(x - starts[0]), starts), data.loc[indiv], label = pop_label[label], color = pop_label.color, linestyle = "-", alpha=alpha)

    plt.plot([0, starts[-1] - starts[0]], [2, 2], "r--")
    line_plot_axis.set_ylim(ymin = 0)

    if max_cp is not None:
        ymax = line_plot_axis.get_ylim()[1]
        line_plot_axis.set_ylim(ymax = min(max_cp, ymax))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gglob_dir", help="Gglob directory")
    parser.add_argument("--plot_dir", default="line_plots")
    parser.add_argument("--outfile", help="Output file name")

    parser.add_argument("--fn_fa", default="/net/eichler/vol7/home/psudmant/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta", help="reference genome fasta file (Default: %(default)s)")
    parser.add_argument("--GC_DTS", dest="fn_GC_DTS", default="/net/eichler/vol7/home/psudmant/genomes/GC_tracks/windowed_DTS/HG19/500_bp_slide_GC", help="GC tracks DTS file (Default: %(default)s")
    parser.add_argument("--DTS_contigs", dest='fn_DTS_contigs', default="/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis/windows/hg19_slide/500_bp_windows.pkl.contigs", help="Contig sizes file (Default: %(default)s)")
    parser.add_argument("--dup_tabix", dest="fn_dup_tabix", default="/net/eichler/vol7/home/psudmant/genomes/annotations/hg19/superdups/superdups.merged.bed.gz", help="Superdups tabix file (Default: %(default)s)")
    parser.add_argument("contig")
    parser.add_argument("start", type = int)
    parser.add_argument("end", type = int)
    parser.add_argument("name")
    parser.add_argument("--pop_file", help="Population manifest file of samples to include")
    parser.add_argument("--exclude_super_pops", nargs="+", help="List of super_pops to exclude (e.g. 'GORILLA CHIMPANZEE')")
    parser.add_argument("--pop_color_scheme", choices = ["hgdp", "primate", "custom"], required=True, help="Choose color scheme")
    parser.add_argument("--custom_label", help="Custom label for color scheme")
    parser.add_argument("--max_cp", type=float, help="Maximum copy number to display (ymax)")
    parser.add_argument("--genes", action="store_true", help="Add gene annotation track?")
    parser.add_argument("--regions", help="Path to bed file of region annotations")
    parser.add_argument("--violins", action="store_true", help="Add violin plot for each region")
    parser.add_argument("--pad", type=int, default=0, help="Number of bp to extend region")
    parser.add_argument("--sort_col", default = None, help="Column in pop_file for sorting samples (Default: %(default)s)")

    args = parser.parse_args()

    if args.outfile is not None:
        outfile_name =args.outfile
    else:
        outfile_name = "%s/%s_%d_%d_%s.pdf" % (args.plot_dir, args.contig, args.start, args.end, args.name)
    bedstring = "%s %d %d %s" % (args.contig, args.start, args.end, args.name)
    region_name = "%s:%d-%d %s" % (args.contig, args.start, args.end, args.name)

    pop_data = pd.read_csv(args.pop_file, sep="\t", header=0)
    pop_data.index = pop_data.sample
    
    if args.exclude_super_pops is not None:
        pop_data = pop_data.loc[~pop_data["super_pop"].isin(args.exclude_super_pops)]

    if args.pop_color_scheme == "hgdp":
        pop_data = pop_data.loc[pop_data["species_label"] == "Human"]
        label = "hgdp_label"
        n_colors = len(pop_data[label].unique().tolist())   
        colors = list(plt.cm.hsv(np.linspace(0, 0.9, n_colors-2))) + list(plt.cm.Greys(np.linspace(0.5, 1, 2)))
        #colors = ["red", "chocolate", "darkgreen", "gold", "darkcyan", "darkblue", "darkorchid"] + list(plt.cm.Greys(np.linspace(0.5, 1, 2)))
    elif args.pop_color_scheme == "primate":
        label = "species_label"
        n_colors = len(pop_data[label].unique().tolist())   
        colors = ["royalblue"] + list(plt.cm.Greys(np.linspace(1, 0.25, 4)))
    else:
        if args.custom_label is None:
            print "Must specify custom label for custom color scheme"
            sys.exit(1)
        label = args.custom_label
        n_colors = len(pop_data[label].unique().tolist())
        colors = list(plt.cm.hsv(np.linspace(0, 0.9, n_colors-2))) + list(plt.cm.Greys(np.linspace(0.5, 1, 2)))

    indivs = pop_data.sample.unique().tolist()

    g = gt.genotyper(args.contig, gglob_dir = args.gglob_dir, plot_dir = args.plot_dir, subset_indivs = indivs, fn_fa=args.fn_fa, dup_tabix = args.fn_dup_tabix, GC_inf = args.fn_GC_DTS)
    X, idx_s, idx_e = g.get_gt_matrix(args.contig, args.start - args.pad, args.end + args.pad)
    starts = g.wnd_starts[idx_s:idx_e]
    data = pd.DataFrame(X, index = g.indivs, columns = map(str, starts))

    pop_data = pop_data.loc[pop_data["sample"].isin(g.indivs)]
    unique_labels = pop_data[label].unique().tolist()
    pop_data["color"] = pop_data[label].map(lambda x: colors[unique_labels.index(x)])
    pop_data["label"] = pop_data[label]

    alpha = 0.6

    fig = plt.figure(figsize=(24,15))

    if args.genes is not None:
        genes_h = 0.1
    else:
        genes_h = 0

    if args.regions is not None:
        regions_h = 0.1
    else:
        regions_h = 0

    # Set axis dimensions
    xmin, ymin = 0.05, 0.05
    lpt_w, lpt_h = 0.8, 0.9 - genes_h - regions_h
    legend_w, legend_h = 0.1, 0.6
    xspace, yspace = 0.045, 0.05
    line_plot_dims = [xmin, ymin + genes_h + regions_h, lpt_w, lpt_h]
    sample_legend_dims = [xmin + lpt_w + xspace, ymin, legend_w, legend_h]

    # Draw legend
    legend_axis = fig.add_axes(sample_legend_dims)
    make_legend(sample_legend_dims, legend_axis, pop_data, alpha=alpha)

    # Draw line plot
    line_plot_axis = fig.add_axes(line_plot_dims)
    make_line_plot(data, region_name, pop_data, line_plot_axis, args.pop_color_scheme, alpha=0.6, max_cp = args.max_cp, sort_col = args.sort_col)

    # Draw gene annotations
    plot_domain = [starts[0], starts[-1]]
    gene_dims = [xmin, ymin + regions_h, lpt_w, genes_h - yspace]
    region_dims = [xmin, ymin, lpt_w, regions_h - yspace]

    if args.genes:
        gene_axis = fig.add_axes(gene_dims)
        make_region_annotation(gene_axis, bedstring, plot_domain, fn_regions = "/net/eichler/vol2/eee_shared/assemblies/hg19/genes/refGene.bed", merge=True, strand=True)

    if args.regions:
        region_axis = fig.add_axes(region_dims)
        region_list = make_region_annotation(region_axis, bedstring, plot_domain, fn_regions = args.regions, merge = False, strand = False)

    plt.savefig(outfile_name, facecolor="w")

    if args.violins:
        fig2, axis = plt.subplots(1, 1, figsize=(15, 15))
        plot_violins(axis, region_list, g, pop_data, colors)
        plt.savefig(outfile_name.replace(".pdf", ".violin.pdf"))
