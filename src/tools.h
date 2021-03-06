#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
public:
    /**
    * A helper method to calculate RMSE.
    */
    static Eigen::VectorXd
    CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

    /**
    * A helper method to calculate Jacobians.
    */
    static void CalculateJacobian(const Eigen::VectorXd &x_state, Eigen::MatrixXd &Hj);

    static Eigen::VectorXd h(const Eigen::VectorXd &x);
};

#endif /* TOOLS_H_ */
