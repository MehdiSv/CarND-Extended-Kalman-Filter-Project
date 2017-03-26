#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"
#include "measurement_package.h"

class KalmanFilter {
public:
    // state vector
    Eigen::VectorXd x_;

    // state covariance matrix
    Eigen::MatrixXd P_;

    // state transistion matrix
    Eigen::MatrixXd F_;

    // process covariance matrix
    Eigen::MatrixXd Q_;

    // measurement noise matrix
    Eigen::MatrixXd R_laser_;
    // measurement function matrix
    Eigen::MatrixXd H_laser_;

    // measurement noise matrix
    Eigen::MatrixXd R_radar_;
    // measurement function matrix
    Eigen::MatrixXd Hj_;

    // identity matrix
    Eigen::MatrixXd I_;

    // external motion
    Eigen::VectorXd u_;

    long long previous_timestamp_;
    int noise_ax_;
    int noise_ay_;
    bool is_initialized_;

    KalmanFilter();

    virtual ~KalmanFilter();

    bool isInitialized();

    void Init(const MeasurementPackage &measurement_pack);

    void InitCovarianceMatrices(long long timestamp);

    void Predict();

    void Update(const Eigen::VectorXd &z);

    void UpdateEKF(const Eigen::VectorXd &z);
};

#endif /* KALMAN_FILTER_H_ */
