# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 09:44:42 2024

@author: zaman.i
"""

import cv2
import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.transform import resize
import tifffile as tiff
from skimage.io import imread, imsave
import glob
from cellpose import models

def get_files_from_folders(folders, image_ext=("*.png", "*.tif"), label_ext="*_seg.npy"):
    image_files = []
    label_files = []
    
    # Loop through all provided folders
    for folder in folders:
        # Loop through each image extension and get matching files
        for ext in image_ext:
            image_files.extend(glob.glob(os.path.join(folder, ext)))
        
        # Get label files
        label_files.extend(glob.glob(os.path.join(folder, label_ext)))
    
    return image_files, label_files

parent_dir = r"/path/to/merscope_train_images"
output_dir = r"/path/to/xenium_output_train_folder"
for file in os.listdir(parent_dir):
    if file.endswith(".npy"):
        path1 = os.path.join(parent_dir,file)
        mask = np.load(os.path.join(parent_dir,file), allow_pickle = True).item()['masks']
        # Resize the mask (nearest-neighbor interpolation)
        scale_factor = 0.2125 / 0.1080
        new_shape = (int(mask.shape[0] / scale_factor), int(mask.shape[1] / scale_factor))
        resized_mask = resize(mask, new_shape, order=0, preserve_range=True)

        # Save the resized mask
        file_name = file[0:-8]
        save_path = os.path.join(output_dir, f"{file_name}_downsampled_seg.npy")
        print(f"Saving: {save_path}")
        np.save(save_path, {'masks': resized_mask})


# Scale factor: MERSCOPE (0.108 µm) → Xenium (0.2125 µm)
scale_factor = 0.2125 / 0.1080  # ≈ 1.9685

for file in os.listdir(parent_dir):
    if file.lower().endswith((".tif", ".tiff", ".png")):
        path = os.path.join(parent_dir, file)
        
        # Load the image
        img = imread(path)
        
        # Compute new shape (downsample → fewer pixels)
        new_shape = (int(img.shape[0] / scale_factor), int(img.shape[1] / scale_factor))
        
        # Resize the image using bilinear interpolation (order=1)
        resized_img = resize(
            img, 
            new_shape, 
            order=1, 
            preserve_range=True, 
            anti_aliasing=True
        ).astype(img.dtype)
        
        # Save the resized image
        file_name = os.path.splitext(file)[0]
        save_path = os.path.join(output_dir, f"{file_name}_downsampled.tif")
        imsave(save_path, resized_img)
        
        print(f"Saved: {save_path}")
        

# Set input/output paths
input_train_folder = r"/path/to/training_images_resized"
input_test_folder = r"/path/to/test_images_resized"


# Define the folder paths
train_data = [input_train_folder]
test_data = [input_test_folder]


# Get list of image files in train and test folders
train_image_files, train_label_files = get_files_from_folders(train_data)
test_image_files, test_label_files = get_files_from_folders(test_data)

# Create empty arrays for train and test images and labels
train_images = []
train_labels = []
test_images = []
test_labels = []

# Loop through train image files and append to train_images array
for image_file in train_image_files:
    # Get corresponding label file
    label_file = os.path.splitext(image_file)[0] + '_seg.npy'
    
    # Check if label file exists
    if not os.path.exists(label_file):
        print(f"Label file {label_file} not found for image file {image_file}")
        continue
    
    # Read image and label
    img = imread(image_file)
    label = np.load(label_file, allow_pickle=True).item()['masks']
    
    # Append to train_images and train_labels arrays
    train_images.append(img)
    train_labels.append(label)

# Loop through test image files and append to test_images array
for image_file in test_image_files:
    # Get corresponding label file
    label_file = os.path.splitext(image_file)[0] + '_seg.npy'
    
    # Check if label file exists
    if not os.path.exists(label_file):
        print(f"Label file {label_file} not found for image file {image_file}")
        continue
    
    # Read image and label
    img = imread(image_file)
    label = np.load(label_file, allow_pickle=True).item()['masks']
    
    # Append to test_images and test_labels arrays
    test_images.append(img)
    test_labels.append(label)
    

# Initialize model and train on combined data
model_type = 'cyto2'
model_name = f'{model_type}_Xenium_withAnnotations_extraROIs_v2'
print(f"Training model with name: {model_name}")
# Initialize Cellpose model
model = models.CellposeModel(gpu=True, model_type=model_type, net_avg=True, nchan=2)
from cellpose.io import logger_setup
logger_setup()


resample=True
net_average = True
channel_axis  = 0 

do_3D = False

#We are using our own normalisation
normalize = True

#average channels
channels = [0,1]

#trained CP model
#model_type = 'cyto2'

model = models.CellposeModel(gpu=True, model_type=model_type,net_avg=net_average,nchan=2)


model_save_path = r"/path/to/model_output"
lr = 0.1
weight_decay=1e-5
n_epochs = 100
name = model_name
model.train(train_data=train_images,train_labels=train_labels,test_data=test_images,test_labels=test_labels,channels=channels,min_train_masks=0,
            normalize=normalize,save_path=model_save_path,learning_rate=lr,weight_decay=weight_decay,n_epochs=n_epochs,model_name=name)


