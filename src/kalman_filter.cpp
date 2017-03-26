#include <iostream>
#include "kalman_filter.h"
#include "tools.h"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
    // state matrix
    x_ = VectorXd(4);
    x_ << 1, 1, 1, 1;

    // state uncertainty covariance matrix
    P_ = MatrixXd(4, 4);
    P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    // the initial transition matrix F_
    F_ = MatrixXd(4, 4);
    F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    // prediction uncertainty covariance matrix
    Q_ = MatrixXd(4, 4);
    Q_ << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;

    // measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    // measurement matrix - laser
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    //measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    // measurement matrix - radar
    Hj_ = MatrixXd(3, 4);
    Hj_ << 0, 0, 0, 0,
            1e-9, 1e-9, 0, 0,
            0, 0, 0, 0;

    // identity matrix
    I_ = MatrixXd::Identity(4, 4);

    // motion vector
    u_ = VectorXd(4);
    u_ << 0, 0, 0, 0;

    previous_timestamp_ = 0.0;
    noise_ax_ = 9;
    noise_ay_ = 9;
    is_initialized_ = false;
}

KalmanFilter::~KalmanFilter() {}

bool KalmanFilter::isInitialized() {
    return is_initialized_;
}

void KalmanFilter::Init(const MeasurementPackage &measurement_pack) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        float x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
        float y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
        x_ << x, y, 0, 0;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
}

void KalmanFilter::InitCovarianceMatrices(long long timestamp) {
    float dt = (timestamp - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = timestamp;

    F_(0, 2) = dt;
    F_(1, 3) = dt;

    float dt4 = pow(dt, 4) / 4;
    float dt3 = pow(dt, 3) / 2;
    float dt2 = pow(dt, 2);

    Q_ << dt4 * noise_ax_, 0, dt3 * noise_ax_, 0,
            0, dt4 * noise_ay_, 0, dt3 * noise_ay_,
            dt3 * noise_ax_, 0, dt2 * noise_ax_, 0,
            0, dt3 * noise_ay_, 0, dt2 * noise_ay_;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_ + u_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    MatrixXd y = z - H_laser_ * x_;
    MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
    MatrixXd K = P_ * H_laser_.transpose() * S.inverse();
    x_ = x_ + K * y;
    P_ = (I_ - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    if (z(0) == 0 && z(1) == 0 && z(2) == 0) {
        return;
    }

    Tools::CalculateJacobian(x_, Hj_);
    MatrixXd y = z - Tools::h(x_);
    MatrixXd S = Hj_ * P_ * Hj_.transpose() + R_radar_;
    MatrixXd K = P_ * Hj_.transpose() * S.inverse();
    x_ = x_ + K * y;
    P_ = (I_ - K * Hj_) * P_;
}
