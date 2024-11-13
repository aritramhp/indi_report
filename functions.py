import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance


def run_rhea(rhea_directory,ref_alpha,ref_phyla,ref_genera,input_dir):
    print('Running Normalization.R')
    res = subprocess.call('Rscript ' + os.path.join(rhea_directory,'1.Normalization','Normalization.R'), shell=True)
    if res:
        print('Error in Normalization.R')
    else:
        print('Running Alpha-Diversity.R')
        res = subprocess.call('Rscript ' + os.path.join(rhea_directory,'2.Alpha-Diversity','Alpha-Diversity.R'), shell=True)
        if res:
            print('Error in Alpha-Diversity.R')
        else:
            df1 = pd.read_csv(os.path.join(rhea_directory,'2.Alpha-Diversity','alpha-diversity.tab'), sep='\t')
            df2 = pd.read_csv(ref_alpha, sep='\t')
            df3 = pd.concat([df1,df2],axis='index')
            alpha_diversity = os.path.join(input_dir,'alpha-diversity_ref_current.tab')
            df3.to_csv(alpha_diversity,sep='\t',index=None)
        
        print('Running Taxonomic-Binning.R')
        res = subprocess.call('Rscript ' + os.path.join(rhea_directory,'4.Taxonomic-Binning','Taxonomic-Binning.R'), shell=True)
        if res:
            print('Error in Taxonomic-Binning.R')
        else:
            df1 = pd.read_csv(os.path.join(rhea_directory,'4.Taxonomic-Binning','Taxonomic-Binning','1.Phyla.all.tab'),sep='\t')
            df2 = pd.read_csv(ref_phyla,sep='\t')
            df3 = pd.merge(df1, df2, on='Unnamed: 0', how='outer').fillna(0)
            tax_phylum = os.path.join(input_dir,'1.Phyla.all_ref_current.tab')
            df3.to_csv(tax_phylum,sep='\t',index=None)
            df1 = pd.read_csv(os.path.join(rhea_directory,'4.Taxonomic-Binning','Taxonomic-Binning','5.Genera.all.tab'),sep='\t')
            df2 = pd.read_csv(ref_genera,sep='\t')
            df3 = pd.merge(df1, df2, on='Unnamed: 0', how='outer').fillna(0)
            tax_genera= os.path.join(input_dir,'5.Genera.all_ref_current.tab')
            df3.to_csv(tax_genera,sep='\t',index=None)
    return alpha_diversity, tax_phylum, tax_genera

def merge_with_ref(ref_alpha,ref_phyla,ref_genera,input_dir):
    df1 = pd.read_csv(os.path.join(input_dir,'alpha-diversity.tab'), sep='\t')
    df2 = pd.read_csv(ref_alpha, sep='\t')
    df3 = pd.concat([df1,df2],axis='index')
    alpha_diversity = os.path.join(input_dir,'alpha-diversity_ref_current.tab')
    df3.to_csv(os.path.join(input_dir,'alpha-diversity_ref_current.tab'),sep='\t',index=None)

    df1 = pd.read_csv(os.path.join(input_dir,'1.Phyla.all.tab'),sep='\t')
    df2 = pd.read_csv(ref_phyla,sep='\t')
    df3 = pd.merge(df1, df2, on='Unnamed: 0', how='outer').fillna(0)
    tax_phylum = os.path.join(input_dir,'1.Phyla.all_ref_current.tab')
    df3.to_csv(os.path.join(input_dir,'1.Phyla.all_ref_current.tab'),sep='\t',index=None)
    df1 = pd.read_csv(os.path.join(input_dir,'5.Genera.all.tab'),sep='\t')
    df2 = pd.read_csv(ref_genera,sep='\t')
    df3 = pd.merge(df1, df2, on='Unnamed: 0', how='outer').fillna(0)
    tax_genera = os.path.join(input_dir,'5.Genera.all_ref_current.tab')
    df3.to_csv(os.path.join(input_dir,'5.Genera.all_ref_current.tab'),sep='\t',index=None)

    return alpha_diversity, tax_phylum, tax_genera


def alpha_plot(df,parameter,my_pal,sample_alpha,yticks,sample_lbl,xlabel,figfile):
    plt.subplots(figsize=(12, 3))
    sns.set_style("whitegrid")
    circle1 = plt.Circle((0, 0), 0.2, color='r')
    ax = sns.boxplot(data=df, x=parameter, y="class",hue='class',palette=my_pal, showfliers = False, legend=False, linewidth=3, medianprops={"linestyle": "--"})
    ax.set_yticks([0,1,2])
    ax.set_yticklabels(yticks, fontsize=30)
    ax.tick_params(axis='x', labelsize=25)
    ax.axvline(sample_alpha, color='#FC7600', alpha=1, linewidth=4.5)
    ax.set(ylabel=None)
    plt.xlabel(xlabel,fontsize=30)

    sample_dates = np.array([sample_alpha])
    for i, x in enumerate(sample_dates):
        plt.text(x, -1, sample_lbl, rotation=0, verticalalignment='center',horizontalalignment='center', color='#FC7600',fontweight='semibold',in_layout=True,fontsize=30)
    plt.savefig(figfile,dpi=200,bbox_inches='tight')
    plt.close()


def alpha_diversity(input_dir, samplename, alpha_div, category, parameter):
    sample_dir = os.path.join(input_dir, samplename)
    if not os.path.exists(os.path.join(input_dir,samplename)): os.mkdir(sample_dir)
    color_library = ['#DFDED4','#F1F4F4','#899499']
    my_pal = {category[i][0]:color_library[i] for i in range(len(category))}
    alpha_list, cat_class, quartile = [], [], []
    for cat, cat_sample in category:
        df_alpha = pd.read_csv(alpha_div, sep='\t').set_index('Unnamed: 0')
        sample_alpha = df_alpha.loc[samplename,parameter]
        cat_alpha = df_alpha.loc[cat_sample,parameter].to_list()
        alpha_list.extend(cat_alpha)
        cat_class.extend([cat]*len(cat_alpha))
        median, q1, q3 = np.percentile(np.array(cat_alpha), 50), np.percentile(np.array(cat_alpha), 25), np.percentile(np.array(cat_alpha), 75)
        quartile.append([cat,'upper quartile'] if sample_alpha > q3 else [cat,'lower quartile'] if sample_alpha < q1 else [cat,'upper half'] if sample_alpha >=median else [cat,'lower half'])
    df = pd.DataFrame(zip(alpha_list,cat_class), columns=[parameter,'class'])
    figfile = os.path.join(sample_dir,'alpha-diversity_EN.jpg')
    alpha_plot(df,parameter,my_pal,sample_alpha,['Healthy young','Healthy old','Obese'],'Sample','Number of Taxa',figfile)
    figfile = os.path.join(sample_dir,'alpha-diversity_DE.jpg')
    alpha_plot(df,parameter,my_pal,sample_alpha,['Gesund, jung','Gesund, alt','Übergewichtig'],'Probe','Zahl der Taxa',figfile)
    return quartile


def stacked_bar_itol(abun_table, sname, samplelist, output_file):
    colors = {'p__Actinobacteriota':'#f5e4ea','p__Bacteroidota':'#94ac7c','p__Firmicutes':'#847cb6','p__Fusobacteriota':'#a28145',
              'p__Proteobacteria':'#e9c334','p__Verrucomicrobiota':'#7b98d1'}
    fp = open(output_file, 'w')
    fp.write('DATASET_MULTIBAR\nSEPARATOR COMMA\nDATASET_LABEL,Rel.abun.: ' + sname + '\nCOLOR,#ff0000\n')
    fp.write('DATASET_SCALE,25-25-#0000ff-1-1-1,50-50-#0000ff-1-1-1,75-75-#0000ff-1-1-1\n')
    fp.write('WIDTH,750\nMARGIN,0\nSHOW_INTERNAL,0\nHEIGHT_FACTOR,2\nBAR_SHIFT,0\nALIGN_FIELDS,0\n')
    abun_data = pd.read_table(abun_table, sep='\t').set_index('Unnamed: 0')[samplelist]
    fp.write('FIELD_LABELS')
    for i in abun_data.index:
        fp.write(','+i)
    fp.write('\nFIELD_COLORS')
    for i in abun_data.index:
        try:
            fp.write(','+colors[i])
        except KeyError:
            fp.write(',#cccccc')
    fp.write('\n')
    
    #---- Data ---#
    fp.write('\nDATA\n')
    for s in abun_data:
        if not s == sname:
            fp.write(s)
            for ab in abun_data[s]:
                fp.write(','+str(ab))
            fp.write('\n')
    fp.close()  

def mark_sample(sname, abun_table, ref_sample, output_file):
    abun_data = pd.read_table(abun_table, sep='\t').set_index('Unnamed: 0')
    # # Euclidean distance
    # dist = [np.linalg.norm(abun_data[s].to_numpy() - abun_data[sname].to_numpy()) for s in ref_sample]
    # Bray-Curtis distance
    dist = [distance.braycurtis(abun_data[s].to_list(), abun_data[sname].to_list()) for s in ref_sample]
    similar_sample = ref_sample[dist.index(min(dist))]
    with open(output_file+'_1.txt', 'w') as fp:
        fp.write('DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,'+sname+'\nCOLOR,#FC7600\nWIDTH,100\nMARGIN,0\nHEIGHT_FACTOR,3\nBAR_SHIFT,0\n')
        fp.write('DATA\n'+similar_sample+',100')
    with open(output_file+'_2.txt', 'w') as fp:
        fp.write('DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,'+sname+'\nCOLOR,#ff0000\nFIELD_LABELS,rl0\nFIELD_COLORS,#FC7600\nFIELD_SHAPES,2\n')
        fp.write('SHOW_INTERNAL,0\nMARGIN,0\nHEIGHT_FACTOR,18\nSYMBOL_SPACING,10\nDATA\n')
        fp.write(similar_sample+',1')
    return similar_sample

def clade_strip(sname, clade_file, output_file):
    clade_data = pd.read_table(clade_file, sep='\t')
    with open(output_file, 'w') as fp:
        fp.write('DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,'+sname+'\nCOLOR,#ff0000\nSTRIP_WIDTH,50\nMARGIN,0\nBORDER_WIDTH,1\nBORDER_COLOR,#000\nSHOW_INTERNAL,0')
        #---- Data ---#
        fp.write('\nDATA\n')
        for idx in clade_data.index:
            fp.write(clade_data.loc[idx,'sample']+','+clade_data.loc[idx,'clade_clr']+'\n')


# plot stacked bar
def stacked_bar_phylum(df_barplot,clr,xticks,ylabel,figfile):
    sns.set(rc={'figure.figsize':(2,10)})
    sns.set_style("whitegrid")
    ax = df_barplot.plot(kind='bar', stacked=True, color=clr,legend='reverse',linewidth=0,rot=90,width=0.8, grid=True)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(xticks)
    ax.tick_params(labelsize=18)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title='',fontsize=18)
    ax.grid(axis='x')
    ax.set_ylabel(ylabel, fontsize = 18)
    plt.savefig(figfile,dpi=200,bbox_inches='tight')
    plt.close()


def box_plot(df_fb,sample_f_by_b,xticks,sample_lbl,figfile):
    plt.subplots(figsize=(2, 7))
    # sns.set(rc={'figure.figsize':(2,7)})
    sns.set_style("whitegrid")
    my_pal = {'healthy_young':'#DFDED4', 'healthy_old':'#F1F4F4'}
    ax = sns.boxplot(data=df_fb, x="class", y="f_by_b", hue="class", palette=my_pal, showfliers = False, legend=False, linewidth=3, medianprops={"linestyle": "--"})
    ax.set_xticks([0,1])
    ax.set_xticklabels(xticks, rotation=90)
    ax.tick_params(labelsize=25)
    plt.ylabel('Bacillota/Bacteroidota',fontsize=25)
    ax.axhline(sample_f_by_b, color='#FC7600', alpha=1, linewidth=4.5)
    ax.set(xlabel=None)
    sample_dates = np.array([sample_f_by_b])
    for i, x in enumerate(sample_dates):
        plt.text(2.5, x, sample_lbl, rotation=0, verticalalignment='center',horizontalalignment='center', color='#FC7600',fontweight='semibold',in_layout=True,fontsize=25)
    plt.savefig(figfile,dpi=200,bbox_inches='tight')
    plt.close()



# plot stacked bar for genus
def stacked_bar_genus(df_barplot,top,colors17,xticks,ylabel,figfile):
    sns.set(rc={'figure.figsize':(2,10)})
    sns.set_style("whitegrid")
    ax = df_barplot.plot(kind='bar', stacked=True, color=colors17[:len(top)]+['#bcbcbc'],legend='reverse',linewidth=0,rot=90,width=0.8, grid=True)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title='',fontsize=18)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(xticks)
    ax.tick_params(labelsize=18)
    ax.grid(axis='x')
    ax.set_ylabel(ylabel, fontsize = 18)
    plt.savefig(figfile,dpi=200,bbox_inches='tight')
    plt.close()
