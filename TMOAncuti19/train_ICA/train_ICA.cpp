#include <opencv2/opencv.hpp>
#include <itpp/itbase.h>
#include <itpp/signal/fastica.h>
#include <itpp/stat/misc_stat.h>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <random>

namespace fs = std::filesystem;
using namespace itpp;
using namespace std;

static constexpr int PATCH_SIZE = 31;
static constexpr int NUM_COMPONENTS = 10;
static constexpr int PATCHES_PER_IMG = 100;
const string DATASET_DIR = "/TMS_ica_images";
const string OUTPUT_FILE = "../ica_model.yml";

// Extract patches from an image
mat extract_patches(const cv::Mat& gray, int patch_size, int count) {
    int rows = gray.rows - patch_size;
    int cols = gray.cols - patch_size;
    mat patches(count, patch_size * patch_size);

    mt19937 rng{ random_device{}() };
    uniform_int_distribution<> dr(0, rows), dc(0, cols);

    for (int i = 0; i < count; ++i) {
        int r = dr(rng), c = dc(rng);
        cv::Mat block = gray(cv::Rect(c, r, patch_size, patch_size));
        block.convertTo(block, CV_64F, 1.0 / 255.0);
        vec v(patch_size * patch_size);
        memcpy(v._data(), block.ptr<double>(), v.size() * sizeof(double));
        patches.set_row(i, v);
    }
    return patches;
}

// Save model in YAML OpenCV format
void save_yaml(const string& fname, const mat& basis, const vec& mu, const vec& sigma) {
    ofstream f{ fname };
    auto write_mat = [&](const string& key, const mat& m) {
        f << key << ": !!opencv-matrix\n"
          << "   rows: " << m.rows() << "\n"
          << "   cols: " << m.cols() << "\n"
          << "   dt: d\n"
          << "   data: [ ";
        for (int i = 0; i < m.rows(); ++i)
            for (int j = 0; j < m.cols(); ++j)
                f << m(i,j)
                  << ((i==m.rows()-1 && j==m.cols()-1) ? " ]\n" : ", ");
    };
    f << "%YAML:1.0\n---\n";
    write_mat("icaBasis", basis);

    // Write mu and sigma
    mat mu_m   = reshape(mu,   mu.size(),    1);
    mat sig_m  = reshape(sigma, sigma.size(), 1);
    write_mat("mu",    mu_m);
    write_mat("sigma", sig_m);
}

// Compute mean of each row
vec rowwise_mean(const mat& M) {
    vec out(M.rows());
    for (int i = 0; i < M.rows(); ++i)
        out[i] = mean(M.get_row(i));
    return out;
}

// Compute standard deviation of each row
vec rowwise_std(const mat& M, const vec& mu) {
    vec out(M.rows());
    for (int i = 0; i < M.rows(); ++i) {
        vec d = M.get_row(i) - mu[i];
        out[i] = sqrt(sum(sqr(d)) / M.cols());
    }
    return out;
}

int main() {
    mat all_patches;
    // Collect patches from a dataset
    for (auto& entry : fs::directory_iterator(DATASET_DIR)) {
        cv::Mat img = cv::imread(entry.path().string(), cv::IMREAD_GRAYSCALE);
        if (img.rows < PATCH_SIZE || img.cols < PATCH_SIZE) continue;
        mat P = extract_patches(img, PATCH_SIZE, PATCHES_PER_IMG);
        all_patches = all_patches.rows() ? concat_vertical(all_patches, P) : P;
    }
    mat data = all_patches.transpose();  

    // Run FastICA
    Fast_ICA ica(data);
    ica.set_nrof_independent_components(NUM_COMPONENTS);
    ica.set_approach(FICA_APPROACH_SYMM);
    ica.set_non_linearity(FICA_NONLIN_TANH);
    ica.set_max_num_iterations(200);
    ica.set_epsilon(1e-4);
    ica.separate();

    // Compute the independent coefficients for all patches
    mat basis = ica.get_separating_matrix();
    mat comps = ica.get_independent_components();
    vec mu    = rowwise_mean(comps);
    vec sigma = rowwise_std(comps, mu);

    // Save the model
    save_yaml(OUTPUT_FILE, basis, mu, sigma);
    
    return 0;
}
