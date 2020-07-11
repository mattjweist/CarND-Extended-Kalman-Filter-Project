#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}
KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Predict from Lidar
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Update Lidar
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  UpdateEither(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Update Radar
  VectorXd hx(3);
  hx << sqrt(x_(0)*x_(0)+  x_(1)*x_(1)),
  		atan2(x_(1),x_(0)),
  		(x_(0)*x_(2) + x_(1)*x_(3))/sqrt(x_(0)*x_(0)+  x_(1)*x_(1));
  VectorXd y = z - hx;
  // ensure angle is between -pi and pi
  while(y(1)<-M_PI || y(1) > M_PI) {
    if(y(1)<-M_PI){
      y(1)+=2*M_PI;
    }
    else {
      y(1)-=2*M_PI;
    }
  }
  
  UpdateEither(y);
}
  
void KalmanFilter::UpdateEither(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New Estimate from Radar
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;
}
