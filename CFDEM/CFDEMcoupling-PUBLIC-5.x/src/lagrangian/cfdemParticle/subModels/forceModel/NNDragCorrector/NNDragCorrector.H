/*---------------------------------------------------------------------------*\
License
    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#pragma once

#include "NNDragCorrectorJiang.H"
#include "forceModel.H"

class NNDragCorrector {
private:
  const dictionary propsDict_;

  const scalar Resys_;
  const scalar fs_;
  const word pFieldName_;
  const volScalarField &p_;
  mutable volVectorField gradP_;
  const bool interpolateGradP_;
  const scalar uTerminal_;
  const scalar alphaMax_;
  const scalar rho_s_;

  const word gravityFieldName_;
  const uniformDimensionedVectorField &g_;
  mutable scalar magG_;
  mutable vector directionG_;
  mutable vector directionGbyGrhoS_;

  mutable volScalarField driftFlux_;
  mutable volScalarField dragCorrection_;

  mutable NNDragCorrectorJiang<float> nnDragJiang_;
  mutable autoPtr<interpolationCellPoint<vector>> gradPInterpolatorPtr_;

public:
  NNDragCorrector(const dictionary &dict, fvMesh const &mesh)
      : propsDict_(dict), Resys_(propsDict_.lookupOrDefault("Resys", -1.)),
        fs_(propsDict_.lookupOrDefault("filter_size", -1.)),
        pFieldName_(propsDict_.lookupOrDefault<word>("pFieldName", "p")),
        p_(mesh.lookupObject<volScalarField>(pFieldName_)),
        gradP_(fvc::grad(p_)),
        interpolateGradP_(propsDict_.lookupOrDefault("interpolation", true)),
        uTerminal_(propsDict_.lookupOrDefault<scalar>("terminal_velocity", -1)),
        alphaMax_(
            propsDict_.lookupOrDefault("MaximumSolidVolumeFraction", 0.64)),
        rho_s_(readScalar(propsDict_.lookup("solidDensity"))),
        gravityFieldName_(
            propsDict_.lookupOrDefault<word>("gravityFieldName", "g")),
        g_(mesh.lookupObject<uniformDimensionedVectorField>(gravityFieldName_)),
        magG_(0.), directionG_(vector::zero), directionGbyGrhoS_(vector::zero),
        driftFlux_(IOobject("driftFlux", mesh.time().timeName(), mesh,
                            IOobject::NO_READ, IOobject::AUTO_WRITE),
                   mesh, dimensionedScalar("driftFlux0", dimVolume/dimArea/dimTime,0.0)),
        dragCorrection_(IOobject("dragCorrection", mesh.time().timeName(), mesh,
                                 IOobject::NO_READ, IOobject::AUTO_WRITE),
                        mesh, dimensionedScalar("dragCorrection0", dimless,0.0)),
        nnDragJiang_(propsDict_.lookupOrDefault<string>(
                         "NNPath",
                         "$CFDEM_SRC_DIR/lagrangian/cfdemParticle/subModels/"
                         "forceModel/"
                         "NNDragCorrector/Jiang2020/DFkerasParameters.nnet"),
                     propsDict_.lookupOrDefault<string>(
                         "MeanPath", "$CFDEM_SRC_DIR/lagrangian/cfdemParticle/"
                                     "subModels/forceModel/"
                                     "NNDragCorrector/Jiang2020/mean.csv"),
                     propsDict_.lookupOrDefault<string>(
                         "StdPath", "$CFDEM_SRC_DIR/lagrangian/"
                                    "cfdemParticle/subModels/forceModel/"
                                    "NNDragCorrector/Jiang2020/std.csv"),
                     propsDict_.lookupOrDefault("verbose", false)) {

    if (uTerminal_ <= 0.) {
      //         uTerminal_ = calculateTerminalVelocity(rho_s_);
      FatalErrorInFunction
          << "Please specifiy \"terminal_velocity\" larger than 0."
          << abort(FatalError);
    }
  }

  ~NNDragCorrector() {}

  void preLoop() const {
    magG_ = mag(g_.value());
    directionG_ = -g_.value() / magG_;
    directionGbyGrhoS_ = directionG_ / (magG_ * rho_s_);
    gradP_ = fvc::grad(p_);
    if (interpolateGradP_) {
      gradPInterpolatorPtr_.reset(new interpolationCellPoint<vector>(gradP_));
    }
  }

  double calcDragCorrection(const label &celli, const scalar &voidfraction,
                            const vector &position, const vector &Ur,
                            const scalar &Vcell, const scalar &dp,
                            const scalar &nuf) const {
    const scalar alphaP = 1. - voidfraction;
    if (!(alphaP > 1e-2) && (alphaP < alphaMax_))
      return 1.;

    const scalar dimlessGradP =
        (interpolateGradP_
             ? gradPInterpolatorPtr_().interpolate(position, celli)
             : gradP_[celli]) &
        directionGbyGrhoS_;

    const scalar Repsys = Resys_ > 0. ? Resys_ : dp * uTerminal_ / nuf;

    const scalar filterLengthScale = Foam::cbrt(Vcell);
    if (filterLengthScale < 1.e-9)
      return 1.; // return if cell volume is tiny
    const scalar FrTerminal = uTerminal_ * uTerminal_ / (dp * magG_);
    const scalar dimlessFilterSize =
        fs_ > 0. ? fs_ : dp * Foam::cbrt(FrTerminal) / filterLengthScale;

    const scalar UrDotDirectionG = Ur & directionG_;
    const scalar dimlessUslip = UrDotDirectionG / uTerminal_;

    const scalar scaledDriftFlux = nnDragJiang_.predict_scaled_drift_flux(
        Repsys, dimlessFilterSize, alphaP / alphaMax_, dimlessGradP,
        dimlessUslip);

    const scalar driftFlux = scaledDriftFlux * alphaMax_ * uTerminal_;
    driftFlux_[celli] = driftFlux;

    const scalar magSqrUr = max(1e-6, magSqr(Ur));
    const scalar dragCorrectionRaw =
        1.0 + driftFlux * UrDotDirectionG / alphaP / magSqrUr;
    const scalar dragCorrectionBounded = min(1.5, max(1e-3, dragCorrectionRaw));
    dragCorrection_[celli] = dragCorrectionBounded;

    return dragCorrectionBounded;
  }
};
