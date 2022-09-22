/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Yundi Jiang, Xiao Chen, Jari Kolehmainen,                 *
 * Ioannis G. Kevrekidis, Ali Ozel, Sankaran Sundaresan,     *
 * Development of data-driven filtered drag model for        *
 * industrial-scale fluidized beds,                          *
 * Chemical Engineering Science, Volume 230, 2021.           *
 * https://doi.org/10.1016/j.ces.2020.116235                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#include "akeras_model.H"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename scalarType = double> class NNDragCorrectorJiang {
private:
  enum ModelProperty {
    ReynoldsNumber = 0,
    DimlessFilterSize = 1,
    ScaledSolidsVolumeFraction = 2,
    DimlessGradP = 3,
    DimlessSlipVelocity = 4,
    NumFeatures = DimlessSlipVelocity + 1
  };

private:
  template <typename T>
  void read_into_vector(std::vector<T> &vec, const std::string filename) {
    std::ifstream reader;
    reader.open(filename);
    if (reader.fail()) {
      std::cerr << "Could not open " << filename << "!" << std::endl;
      std::abort();
    }
    for (T &elem : vec) {
      reader >> elem;
    }
    reader.close();
  }

private:
  akeras::aKerasModel DFnnModel_;
  const int num_features_;
  std::vector<scalarType> means_;
  std::vector<scalarType> stds_;
  akeras::DataChunkFlat dataChunk_;

public:
  NNDragCorrectorJiang(const std::string _nnfile, const std::string _meansfile,
                       const std::string _stdsfile, const bool _verbose = false)
      : DFnnModel_(_nnfile, _verbose),
        num_features_(ModelProperty::NumFeatures), means_(num_features_),
        stds_(num_features_), dataChunk_(num_features_, 0.) {
    if (_verbose) {
      std::cout << "Loading means from " << _meansfile << "\nand stds from "
                << _stdsfile << std::endl;
    }
    read_into_vector(means_, _meansfile);
    read_into_vector(stds_, _stdsfile);
  }

  double
  predict_scaled_drift_flux(const scalarType &Resys,
                            const scalarType &dimless_filter_size,
                            const scalarType &scaled_solid_volume_fraction,
                            const scalarType &dimless_grad_P,
                            const scalarType &dimless_slip_velocity) {
    // to make this thread-safe, move dataChunk creation in here
    auto &features = dataChunk_.get_1d_rw();

    // collect features and normalize the features using mean and std input
    features[ReynoldsNumber] =
        (Resys - means_[ReynoldsNumber]) / stds_[ReynoldsNumber];
    features[DimlessFilterSize] =
        (dimless_filter_size - means_[DimlessFilterSize]) /
        stds_[DimlessFilterSize];
    features[ScaledSolidsVolumeFraction] =
        (scaled_solid_volume_fraction - means_[ScaledSolidsVolumeFraction]) /
        stds_[ScaledSolidsVolumeFraction];
    features[DimlessGradP] =
        (dimless_grad_P - means_[DimlessGradP]) / stds_[DimlessGradP];
    features[DimlessSlipVelocity] =
        (dimless_slip_velocity - means_[DimlessSlipVelocity]) /
        stds_[DimlessSlipVelocity];

    // compute prediction
    const auto prediction = DFnnModel_.compute_output(&dataChunk_);
    const auto scaled_dimless_drift_flux = prediction[0];

    return scaled_dimless_drift_flux;
  }
};
