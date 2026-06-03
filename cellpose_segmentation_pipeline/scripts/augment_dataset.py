import os
import glob
import shutil
import numpy as np
from skimage.io import imread, imsave
import albumentations as A


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

INPUT_FOLDER = "path/to/training_images"
OUTPUT_FOLDER = "path/to/output_augmented_folder"


# -----------------------------------------------------------------------------
# Augmentation
# -----------------------------------------------------------------------------

augmentation = A.Compose([
    [A.GaussianBlur(blur_limit=(3, 5), p=0.5), 
     A.ColorJitter(brightness=(0.8, 1.2), contrast=(0.8, 1.2), saturation=(0.8, 1.2), hue=(0, 0), p=0.5), 
     A.GaussNoise(var_limit=(10.0, 50.0), p=0.5)],
])


# -----------------------------------------------------------------------------
# Create output folder
# -----------------------------------------------------------------------------

if os.path.exists(OUTPUT_FOLDER):
    shutil.rmtree(OUTPUT_FOLDER)

os.makedirs(OUTPUT_FOLDER)


# -----------------------------------------------------------------------------
# Process images
# -----------------------------------------------------------------------------

image_files = glob.glob(os.path.join(INPUT_FOLDER, "*.png"))

print(f"Found {len(image_files)} images")

for image_file in image_files:

    image_name = os.path.splitext(os.path.basename(image_file))[0]

    mask_file = os.path.splitext(image_file)[0] + "_seg.npy"

    if not os.path.exists(mask_file):
        print(f"Skipping {image_name}: mask not found")
        continue

    # Load image
    image = imread(image_file)

    # Load Cellpose mask
    mask = np.load(
        mask_file,
        allow_pickle=True
    ).item()["masks"]

    # Apply augmentation
    augmented = augmentation(
        image=image,
        mask=mask
    )

    aug_image = augmented["image"]
    aug_mask = augmented["mask"]

    # Save outputs
    out_image = os.path.join(
        OUTPUT_FOLDER,
        f"{image_name}_aug.png"
    )

    out_mask = os.path.join(
        OUTPUT_FOLDER,
        f"{image_name}_aug_seg.npy"
    )

    imsave(out_image, aug_image)

    np.save(
        out_mask,
        {"masks": aug_mask}
    )

    print(f"Saved {image_name}")


print("Done!")