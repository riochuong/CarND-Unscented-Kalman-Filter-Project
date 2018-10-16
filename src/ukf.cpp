#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_; 
  weights_ = VectorXd(n_aug_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_aug_; i++) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_)) ; 
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  // generate sigma points
  // Augmentation x
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_ + 1) = 0;
  x_aug(n_x_ + 2) = 0;
  
  // augmentation P
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0 ,
       0 , std_yawdd_ * std_yawdd_;
  
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q;
  MatrixXd sqrt_P_aug = P_aug.llt().matrixL(); // cholevsky square matrix
  
  // get the X_sigma_aug
  MatrixXd X_sigma_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  X_sigma_aug.col(0) = x_aug;
  for (int i = 1; i < n_aug_ + 1; i++) {
      X_sigma_aug.col(i) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i - 1); 
      X_sigma_aug.col(i + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i - 1); 
  }

  MatrixXd X_sigma_pred = MatrixXd(5, 2*n_aug_+1);
  X_sigma_pred.fill(0.0);
  // Now we can do prediction from each sigma point to generate the 
  // sigma points mapping for the predicted step
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      VectorXd x_pred = VectorXd(n_x_);
      float yawd = X_sigma_aug(n_x_ - 1, i);
      float p_x = X_sigma_aug(0, i);
      float p_y = X_sigma_aug(1, i);
      float v = X_sigma_aug(2, i);
      float yaw = X_sigma_aug(3, i);
      float nu_a = X_sigma_aug(n_x_, i);
      float nu_yawdd = X_sigma_aug(n_x_ + 1, i);

      if (yawd < fabs(0.01)) {
         x_pred(0) = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
         x_pred(1) = p_y + v / yawd * (-(1) * cos(yaw + yawd * delta_t) + cos(yaw));
         
      } else {
         x_pred(0) = p_x + v * cos(yaw) * delta_t;
         x_pred(1) = p_y + v * sin(yaw) * delta_t;
      
      }
      x_pred(0) = x_pred(0) +  1/2 * (delta_t * delta_t) * cos(yaw) * nu_a;
      x_pred(1) = x_pred(1) +  1/2 * (delta_t * delta_t) * sin(yaw) * nu_a;
      x_pred(2) = v + delta_t * nu_a;
      x_pred(3) = yaw + yawd * delta_t + 1/2 * (delta_t * delta_t) * nu_yawdd;
      x_pred(4) = yawd + 0 + delta_t * nu_yawdd;
      // assign the vector to its correct column
      X_sigma_pred.col(i) = x_pred;
  }
  // calcualte new x
  x_.fill(0.0); 
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      x_ = x_ + weights_(i) * X_sigma_pred.col(i);
  }
  // calculate new covariance matrix 
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd x_diff = X_sigma_pred.col(i) - x_;      
      // normalize angle so it wil lbe between pi and -pi
      while (x_diff(3) > M_PI) x_diff(3) -=  2 * M_PI;
      while (x_diff(3) < M_PI) x_diff(3) += 2 *M_PI;
      P_ = P_ + weights_(i) * (x_diff * x_diff.transpose());
  }
  Xsig_pred_ = X_sigma_pred; 
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // convert sigma points to radar 
  //
  MatrixXd R = MatrixXd (3,3);
  R << std_radr_ * std_radr_, 0 , 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  for (int i; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_sig = VectorXd(3);
     double p_x = Xsig_pred_(0, i);
     double p_y = Xsig_pred_(1, i);
     double yaw = Xsig_pred_(3, i);
     double v_x = Xsig_pred_(2, i) * cos(yaw);
     double v_y = Xsig_pred_(2, i) * sin(yaw);
     z_sig(0) = sqrt((p_x *p_x) + (p_y * p_y));
     z_sig(1) = atan2(p_y, p_x);
     z_sig(2) = ((p_x * v_x) + (p_y * v_y)) / (sqrt(p_x*p_x + p_y*p_y));
     while (z_sig(1) > M_PI) z_sig(1) -= 2 * M_PI;
     while (z_sig(1) < M_PI) z_sig(1) += 2 * M_PI;
     z_pred += weights_(i) * z_sig;
     Zsig.col(i) = z_sig;
  }

  // Calculate Innovation S
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
      while(z_diff(1) < M_PI) z_diff(1) += 2 * M_PI;
      S = S + weights_(i) * z_diff * z_diff.transpose(); 
  }
  // Add noise
  S += R;
  // Correlation matrix between State and measurement
  MatrixXd T = MatrixXd(n_x_, 3);
  T.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     VectorXd z_diff = Zsig - z_pred;

     // normalize angle
     while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
     while(x_diff(3) < M_PI) x_diff(3) += 2 * M_PI;
     while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
     while(z_diff(1) < M_PI) z_diff(1) += 2 * M_PI;
     T = T + weights_(i) * x_diff * z_diff.transpose(); 
  }
  MatrixXd K = T * S.inverse();
  // need to normalsize angle diff here
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
  while(z_diff(1) < M_PI) z_diff(1) += 2 * M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

}
