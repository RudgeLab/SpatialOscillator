import numpy as np
import matplotlib.pyplot as plt
import cv2
import json
import os
import tempfile
from pdf2image import convert_from_path, convert_from_bytes

path_pdfs = 'pdfs/1.0_1000.0/'
files = os.listdir(path_pdfs)
files.sort()
#files
save_dir = 'pngs/1.0_1000.0/'

# pdf to png
for filename in files:
    with tempfile.TemporaryDirectory() as path:
         images_from_path = convert_from_path(path_pdfs+filename, dpi=300, output_folder=path, last_page=1, first_page =0)

    base_filename  =  os.path.splitext(os.path.basename(filename))[0] + '.png'     

    

    for page in images_from_path:
        page.save(os.path.join(save_dir, base_filename), 'PNG')
