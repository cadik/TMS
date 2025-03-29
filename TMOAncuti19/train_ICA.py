import cv2
import numpy as np
from sklearn.decomposition import FastICA
import os
import glob
import yaml

patch_size = 31
num_components = 10
patches_per_image = 100
dataset_dir = '/TMS_build/images_subset'

# Save model in YAML OpenCV format.
def save_model_opencv_format(filename, icaBasis, mu, sigma):
    with open(filename, "w") as f:
        f.write("%YAML:1.0\n")
        f.write("---\n")
        
        def write_matrix(key, mat):
            rows, cols = mat.shape
            dt = "d"
            data_str = ", ".join(map(str, mat.flatten()))
            f.write(f"{key}: !!opencv-matrix\n")
            f.write(f"   rows: {rows}\n")
            f.write(f"   cols: {cols}\n")
            f.write(f"   dt: {dt}\n")
            f.write(f"   data: [ {data_str} ]\n")
        
        write_matrix("icaBasis", icaBasis)
        write_matrix("mu", mu)
        write_matrix("sigma", sigma)

# Extract patches from an image.
def extract_patches(image, patch_size, num_patches):
    patches = []
    # Convert to grayscale if needed.
    if len(image.shape) == 3:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    else:
        gray = image.copy()
    # Normalize to [0,1] and convert to float64.
    gray = gray.astype(np.float64) / 255.0
    h, w = gray.shape
    for _ in range(num_patches):
        r = np.random.randint(0, h - patch_size)
        c = np.random.randint(0, w - patch_size)
        patch = gray[r:r+patch_size, c:c+patch_size]
        patches.append(patch.flatten())
    return np.array(patches)

# Collect patches from a dataset.
all_patches = []
image_files = glob.glob(os.path.join(dataset_dir, '*.*'))
for filepath in image_files:
    image = cv2.imread(filepath)
    if image is None:
        continue
    patches = extract_patches(image, patch_size, patches_per_image)
    all_patches.append(patches)
all_patches = np.vstack(all_patches)
print("Total patches extracted:", all_patches.shape[0])

# Run FastICA.
ica = FastICA(n_components=num_components, random_state=0)
ica.fit(all_patches)
ica_basis = ica.components_

# Compute the independent coefficients for all patches.
coeffs = ica.transform(all_patches)
# Compute mu and sigma from the independent coefficients.
mu = np.mean(coeffs, axis=0).reshape(-1, 1)
sigma = np.std(coeffs, axis=0).reshape(-1, 1)

# Save the model
save_model_opencv_format("ica_model.yml", ica_basis, mu, sigma)
print("ICA model saved to ica_model.yml")
