#include "FusionEKF.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    if (!ekf_.isInitialized()) {
        ekf_.Init(measurement_pack);
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    ekf_.InitCovarianceMatrices(measurement_pack.timestamp_);

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
