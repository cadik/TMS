# Train ICA

This tool generates the ICA model used for the global weight map in the Ancuti19 method, which is described in:

> C. Ancuti, C. O. Ancuti, M. Feixas and M. Sbert, “Image Decolorization Based on Information Theory,” *2019 IEEE International Conference on Image Processing (ICIP)*, Taipei, 2019, pp. 3242–3246, doi:10.1109/ICIP.2019.8803485. Available at: [https://ieeexplore.ieee.org/document/8803485](https://ieeexplore.ieee.org/document/8803485)

It samples random patches from a directory of grayscale images, runs FastICA to learn a set of independent basis functions, and writes the learned basis, component means, and variances to an OpenCV-compatible YAML file (`../ica_model.yml`).

## 1. Prerequisites

Ensure the following libraries and tools are installed:

- **it++** (≥ 4.3): [IT++ documentation](https://itpp.sourceforge.net/4.3.1/)
- **OpenCV** (≥ 4.0)  
- **BLAS** & **LAPACK**  
- **pkg-config**, **make**, **g++** (C++17)

On Fedora you can install all of these uring:

```bash
sudo dnf install itpp-devel opencv-devel blas-devel lapack-devel
```

## 2. Data Preparation

Create a directory called `images_subset/` in the project root. Populate it with at least ~ 3 500 grayscale images (JPG, PNG, etc.), each ≥ 31×31 px. For example, a random subset of ImageNet images works well for training.

## 3. Build

Simply run make, which produces the train_ICA executable.

```bash
make
```

## 4. Run

Execute the training program:

```bash
./train_ICA
```
This scans the `/images_subset/` directory, extracts 31×31 patches (100 per image), runs FastICA with 10 components, and writes the `ica_model.yml` file into the parent directory (`../ica_model.yml`).

Note: Training may take several minutes.
