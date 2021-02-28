import shelve
import numpy as np
from matplotlib import use
use("Agg") # Backend for writing to files
from matplotlib.pyplot import figure, axes, legend, close
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.gridspec as gridspec
from scipy import stats
import pandas as pd

# Function to return PDF Summary PLots
def RiboScore_PDF(df,gini_dict,pathToPDF,StudyName):
    with PdfPages(pathToPDF + "/" + "RiboScore_{}.pdf".format(StudyName)) as pdf:
	fig_num = 0
	ribofile = open("RiboScore.csv","a")
	for index, row in df.iterrows():
	    fig_num += 1
	    file_id = index
	    title_1 = "RiboScore: {}".format(file_id)
	    
	    ### Store Plot Variables
	    # Gene Body Distribution
	    gbd_values = [float(i) for i in row[65:68]]
	    
	    # Gini
	    reads = gini_dict[file_id]
	    reads = np.array(reads)
	    sorted_reads = np.sort(reads)
	    total = sum(reads)
	    step_x = float(total) / float(len(reads))
	    bin_num = 20
	    step_bin = int(len(reads) / bin_num)
	    start_bin = (len(reads) % bin_num) 
	    
	    x = np.arange(start_bin+step_bin,len(reads)+1,step_bin)
	    
	    equal_y_bin = x.copy()
	    equal_y = np.array(equal_y_bin) / float(len(reads))
	    
	    total_file_y = np.cumsum(sorted_reads)
	    file_y_bin = total_file_y[start_bin+step_bin-1::step_bin]
	    file_y = np.array(file_y_bin) / float(total)
	    
	    x = np.insert(x,0,0)
	    equal_y = np.insert(equal_y,0,0)
	    file_y = np.insert(file_y,0,0)
	    
	    # Triplet Periodicity
	    tp_values = [float(i) for i in row[61:64]]
			      
	    # RUST
	    # Index candidates for A-site Peak
	    Asite_peak = [float(row[39]),float(row[40]),float(row[41])]
	    
	    # Initialise Ideal Vector
	    vec_peak = 0.5
	    ideal_vec = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
	    
	    # Insert A-site peak into Idealized vector
	    if Asite_peak[0] == max(Asite_peak) and Asite_peak[1] == max(Asite_peak):
		ideal_vec[Asite_peak.index(max(Asite_peak))+5] = vec_peak
	    elif Asite_peak[0] == max(Asite_peak) and Asite_peak[2] == max(Asite_peak):
		ideal_vec[Asite_peak.index(max(Asite_peak))+5] = vec_peak
	    else:
		ideal_vec[Asite_peak.index(max(Asite_peak))+4] = vec_peak
		
	    rust_ref = [0] * 35 + ideal_vec + [0] * 14
	    rust_ref = np.array(rust_ref)
	    
	    rust_test = [float(i) for i in row[:60]]
	    rust_test = np.array(rust_test)
	    
	    # Feature Scores
	    std_values = [float(row[68]),float(row[69]),float(row[64]),float(row[60])]
	    
	    ### Create Summary Output
	    # Initialise Figure
	    plt.figure(fig_num,figsize=(9, 12))
	    plt.suptitle(title_1,weight="bold", size=18)
	    plt.subplots_adjust(wspace=0.35, hspace=0.35)
	    
	    # Specify Figure Layout
	    gs = gridspec.GridSpec(3, 3)
	    ax1 = plt.subplot(gs[0,0])
	    ax2 = plt.subplot(gs[0,1:])
	    ax3 = plt.subplot(gs[1,0])
	    ax4 = plt.subplot(gs[1,1:])
	    ax5 = plt.subplot(gs[2,:-1])
	    ax6 = plt.subplot(gs[2,2])
	    
	    # Barplot: Gene Body Distribution
	    gbd_labels = ("5' Leader","CDS","3' Trailer")
	    gbd_pos = np.arange(len(gbd_labels))
	    barlist1 = ax1.bar(gbd_pos, gbd_values, align='center', alpha=0.5)
	    barlist1[0].set_color('b')
	    barlist1[1].set_color('g')
	    barlist1[2].set_color('r')
	    ax1.set_title("Gene Body Distribution")
	    ax1.set_xticks(gbd_pos, minor=False)
	    ax1.set_xticklabels(gbd_labels, fontdict=None, minor=False)
	    ax1.set_yticks([.2,.4,.6,.8,1],minor=False)
	    ax1.set_ylabel("Proportion Of Reads")
	    
	    # Lineplot: Footprint Distribution (Gini)
	    ax2.set_xlim([0,len(reads)])
	    ax2.set_ylim([0,1])
	    ax2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
	    ax2.plot(x, equal_y, color='b',label = "Equal Dist")
	    ax2.plot(x, file_y, color='g',label = "File Dist")
	    ax2.fill_between(x, equal_y, file_y, where= file_y <= equal_y, facecolor='palegreen', interpolate=True)
	    ax2.legend(loc='upper left',frameon=False)
	    ax2.set_title("CDS Footprint Distribution")
	    ax2.set_xlabel("Cumulative Read Density")
	    ax2.set_ylabel("Proportion Of Total Reads")
	    
	    # Barplot: Triplet Periodicity
	    tp_labels = ('Max','Mid','Min')
	    tp_pos = np.arange(len(tp_labels))
	    barlist3 = ax3.bar(tp_pos, tp_values, align='center', alpha=0.5)
	    barlist3[0].set_color('g')
	    barlist3[1].set_color('b')
	    barlist3[2].set_color('r')
	    ax3.set_title("Triplet Periodicity")
	    ax3.set_xticks(tp_pos, minor=False)
	    ax3.set_xticklabels(tp_labels, fontdict=None, minor=False)
	    ax3.set_yticks([.2,.4,.6,.8,1],minor=False)
	    ax3.set_ylabel("Sub-Codon Proportion")
		    
	    # Lineplot: RUST
	    rust_x = list(range(-40,20))
	    rust_x = np.array(rust_x)
	    ax4.axis([-40, 20, 0, 1])
	    ax4.plot(rust_x,rust_ref,color='b',label = "Ref")
	    ax4.plot(rust_x,rust_test,color='g',label = "File")
	    ax4.legend(loc='upper right')
	    ax4.set_title("Metafootprint Profile")
	    ax4.set_xlabel("Codon Position")
	    ax4.set_ylabel("Kullback Liebler Divergence")
	    
	    # Barplot: Standardised Feature Scores
	    std_labels = ('GBD','Gini','Periodicity','RUST')
	    std_pos = np.arange(len(std_labels))
	    barlist5 = ax5.bar(std_pos, std_values, align='center', alpha=0.5)
	    barlist5[0].set_color('g')
	    barlist5[1].set_color('b')
	    barlist5[2].set_color('r')
	    barlist5[3].set_color('purple')
	    ax5.set_title("Standardised Feature Scores")
	    ax5.set_xticks(std_pos, minor=False)
	    ax5.set_xticklabels(std_labels, fontdict=None, minor=False)
	    ax5.set_yticks([.2,.4,.6,.8,1],minor=True)
	    ax5.set_ylabel("Score")
	    
	    # Text Box: QC Scores
	    GBD = float(std_values[0])
	    GINI = float(std_values[1])
	    PERIOD = float(std_values[2])
	    RUST = float(std_values[3])
	    SS = round(sum([GBD,GINI,PERIOD,RUST])/4,2)
	    summary = int(SS*100)
	    
	    props = dict(boxstyle='round', facecolor='lightskyblue', alpha=0.5)
	    textstr1 = '      QC Scores    '
	    textstr2 = "\n\nGBD: {}\nGini: {}\nPeriodicity: {}\nRUST: {}\n\nSummary Score: {}".format(GBD,GINI,PERIOD,RUST,summary)
	    ax6.axis('off')
	    ax6.text(-.15, .05, textstr1 + textstr2, transform=ax6.transAxes, fontsize = 14, bbox = props)
		    
	    # Save Figure
	    pdf.savefig(fig_num,tight_layout = True)
	    plt.close(fig_num)
	    
	    # TEMP #
	    # Write out to RiboQC.csv
	    ribofile.write("{},{},{},{},{},{},{}\n".format(StudyName,file_id,GBD,GINI,PERIOD,RUST,summary))
	ribofile.close()
	print "Files Added to Pdf: {}".format(fig_num)
	close()
    return

# Functiom to return logo figues
def RiboScore_Logo(df,pathToPDF):  
    logo_num = 0
    for index, row in df.iterrows():
	file_id = index
	logo_num += 1
	labels = "GBD","GINI","TP","RUST",""
	GBD, GINI, TP, RUST = float(row[68]),float(row[69]),float(row[64]),float(row[60])
	SS = round(sum([GBD,GINI,TP,RUST])/4,2)
	measures = [round(float(GBD/4),2),round(float(GINI/4),2),round(float(TP/4),2),round(float(RUST/4),2),round(float(1-SS),2)]
	patches,texts = plt.pie(measures,colors=['green','blue','red','purple','white'],counterclock=False,startangle=90)
	my_circle=plt.Circle( (0,0), 0.95, color='white')
	p=plt.gcf()
	p.gca().add_artist(my_circle)
	gbd_patch = mpatches.Patch(color='green',label='GBD')
	gini_patch = mpatches.Patch(color='blue',label='GINI')
	tp_patch = mpatches.Patch(color='red',label='TP')
	rust_patch = mpatches.Patch(color='purple',label='RUST')
	legend = plt.legend(handles=[gbd_patch,gini_patch,tp_patch,rust_patch],frameon=False,ncol=4,loc=8)
	plt.setp(legend.get_texts(),family='serif',color='gray',alpha=0.75)
	p.gca().add_artist(legend)
	font1 = {'family':'serif','color':'gray','weight':'normal','size':32,}
	font2 = {'family':'serif','color':'black','weight':'normal','size':64,}
	plt.text(0,0,'Summary Score',va='center',ha='center', fontdict=font1,alpha=0.25)
	summary = int(SS*100)
	plt.text(0,0,summary,va='center',ha='center',fontdict=font2)
	plt.savefig(pathToPDF + "/" + 'RiboScore_{}.png'.format(file_id),tight_layout = True)
	plt.close()
    print "Logos Generated: {}".format(logo_num)
    return    

