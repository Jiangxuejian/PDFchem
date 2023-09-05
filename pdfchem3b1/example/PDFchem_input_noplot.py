# %%
import numpy as np
import math, matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
#from matplotlib import colors, cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.offsetbox import AnchoredText
import os, glob, datetime
from scipy.interpolate import griddata
# import concurrent.futures
# Suppress all warnings
import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('./')
matplotlib.rc('font',size=17)

# ------------------------------------------------------
# -------------- User defined parameters ---------------
# Check if there are AV-PDF file(s) provided by the user
# avpdf_dir =  os.path.normpath('./avpdf_input')
avpdf_dir =  os.path.normpath('../../crop_image/67pixel_pdf')
redo_calculation = 1

# -------------- User defined parameters ---------------
# ------------------------------------------------------


# %%
avpdf_files= sorted(glob.glob(os.path.join(avpdf_dir,'*.dat')), reverse=1)
#if os.path.exists(avpdf_files[0]):
if len(avpdf_files) > 0:
    print(f'User PDF data file available.')
    pdf_type = 'user'
else:
    print(f'Making a PDF based on parameters av_bar={av_bar} and width={s}.')
    pdf_type = 'model'

def screen_time():
    return datetime.datetime.now().strftime("\033[0;32m%H:%M:%S %Y-%m-%d > \33[0m")


#%% Parallel Running PDFchem
import pdfchem
def process_file(avpdf_file):
    file_index =  avpdf_file.split('/')[-1].split('.')[0][-6:]
    pdfchem_output_file = os.path.join(f'./pdfchem_output', f'output{file_index}.dat')
    if not os.path.isfile(pdfchem_output_file):
        pdfchem.main(avpdf_file, pdfchem_output_file)
    # pdf_i =  avpdf_file.split('/')[-1].split('.')[0][-6:]
    #print(f'{screen_time()}{pdf_i}: Reading PDF file "{avpdf_file}" and running PDFchem...\n')
    # pdfchem_output_file = os.path.join(f'pdfchem_output', f'output{pdf_i}.dat')
    # pdfchem
    #subprocess.call(['./pdfchem_algorithm_PDFinput', avpdf_file, pdfchem_output_file])
    #os.system(f'./pdfchem_algorithm_PDFinput {avpdf_file} output_{pdf_i}.dat')
    #shutil.copy('output.dat', pdfchem_output_file)
    return

print(f'{screen_time()}: Started.')
if not os.path.exists('./pdfchem_output'): os.mkdir('./pdfchem_output')
if redo_calculation == 1:
    for f in avpdf_files:
        process_file(f)
print(f'{screen_time()}: Finished.')
