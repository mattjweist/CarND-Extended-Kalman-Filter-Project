#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 	1., 0., 0., 0.,
  				0., 1., 0., 0.;
  
  // set the acceleration noise components
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
  
  // intial state vector
  VectorXd x(4);
  x << 0.0, 0.0, 0.0, 0.0;

  // initial uncertainties
  MatrixXd P(4, 4);
  P << 	1, 0, 0, 0,
  		0, 1, 0, 0,
  		0, 0, 1000, 0,
  		0, 0, 0, 1000;

  // the initial transition matrix F
  MatrixXd F(4, 4);
  F << 	1, 0, 1, 0,
  		0, 1, 0, 1,
  		0, 0, 1, 0,
  		0, 0, 0, 1;
  
  // initial noise
  MatrixXd Q(4, 4);
  Q << 	0, 0, 0, 0,
  		0, 0, 0, 0,
  		0, 0, 0, 0,
  		0, 0, 0, 0;
  
  ekf_.Init(x,P,F,H_laser_,R_laser_,Q);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*** 
  Initialization
   */
  
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      ekf_.x_ << 	measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]),
      				measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]),
      				measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]),
      				measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);
      // Velocities are known, uncertainty is low
      ekf_.P_(2,2) = 1;
      ekf_.P_(3,3) = 1; 
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state
      ekf_.x_ << 	measurement_pack.raw_measurements_[0],
      				measurement_pack.raw_measurements_[1],
      				0,
      				0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
         0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
         dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
         0, dt_3/2*noise_ay_, 0, dt_2*noise_ay_;
  
  // Call predict function
  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
  	ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
  	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
