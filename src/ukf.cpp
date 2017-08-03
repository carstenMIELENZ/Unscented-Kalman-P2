#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_       = 2.0;  

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3; 

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;  // dont got first measurements

  n_x_     = 5         ;// state dimension
  n_aug_   = 7         ;// aug dimension
  lambda_  = 3 - n_aug_;// spreading factor

  // Signma points matrix 
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  weights_   = VectorXd(2*n_aug_+1);       
  
  //set weights
  weights_(0)   = lambda_/3.0;
  for (int i=1; i<2*n_aug_+1; i++)  //2n+1 weights
    weights_(i) = 1/6.0;
  
  // covariance
  P_ <<    1.0,    0.0,    0.0,    0.0,    0.0,
           0.0,    1.0,    0.0,    0.0,    0.0,
           0.0,    0.0,    1.0,    0.0,    0.0,
           0.0,    0.0,    0.0,    1.0,    0.0,
           0.0,    0.0,    0.0,    0.0,    1.0; 

  // NIS
  NIS_radar_num_     =  0.0;
  NIS_radar_num_5_   =  0.0;
  NIS_lidar_num_     =  0.0;
  NIS_lidar_num_5_   =  0.0;
}

UKF::~UKF() { }


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
  
  //-------------------------------------------------------------------------------
  //
  // code entered July-23-2017
  //
  //-------------------------------------------------------------------------------

  // Initialization
  
  if (!is_initialized_) {
    
     std::cout << std::endl << " NOTE (ProcessMeasurement): init ukf" << endl;
    
    // set start time
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
      double px = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
      double py = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
      x_ << px, py, 0.0, 0.0, 0.0;
    
      is_initialized_ = true;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER ) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      x_ << px, py, 0.0, 0.0, 0.0;
      
      //set covariance
      P_(0,0) = std_laspx_ * std_laspx_;  
      P_(1,1) = std_laspy_ * std_laspy_;  
      
      is_initialized_ = true;
    }
   
    return;
  }
  
  
  // Select between RADAR and LASER assumption either RADAR or LASER exists
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {  
  
     // Set time step
     double delta_t  = (meas_package.timestamp_ - time_us_)/1000000.0;
            time_us_ =  meas_package.timestamp_;        
  
     // Predict 
     Prediction(delta_t);
     // Update
     UpdateRadar(meas_package); 
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
  
     // Set time step
     double delta_t  = (meas_package.timestamp_ - time_us_)/1000000.0;
            time_us_ =  meas_package.timestamp_;        
 
     // Predict 
     Prediction(delta_t);
     // Update
     UpdateLidar(meas_package); 
  }
 
  return;
    
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

  //*************************************************************************
  // A) GENERATION OF SIGMA POINTS  
  //*************************************************************************

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  
  //create augmented mean state * NOTE: (5),(6) is 0 because mean = 0
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_)               = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1)           = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  const double SQRL    = sqrt(3.0);  // = 3 - n_aug_ + n_aug_ 
  Xsig_aug.col(0)      = x_aug;
  for (int i=0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + SQRL * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - SQRL * L.col(i);
  }

  //************************************************************************
  // B) PREDICTION SIGMA POINTS  
  //************************************************************************

  //predict sigma points
  for (int i=0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin(yaw+yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p  = v_p  + nu_a*delta_t;

    yaw_p  = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //************************************************************************
  // C) PREDICTED MEAN & COVARIANCE  
  //************************************************************************
 
  //get mean x
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  //predicted state mean
  for (int i=0; i < 2 * n_aug_ + 1; i++)  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);

  //get mean covariance
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  //predicted state covariance matrix
  for (int i=0; i < 2*n_aug_+1; i++)  //iterate over sigma points
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //set prediction 
  x_ = x;
  P_ = P;
  
  return;
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

  //************************************************************************
  // A) PREDICT LIDAR MEAN & COVARIANCE  
  //************************************************************************

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(2, 2*n_aug_+1);
  Zsig.row(0) = Xsig_pred_.row(0); 
  Zsig.row(1) = Xsig_pred_.row(1); 

  //mean predicted measurement
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) 
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)  //2n+1 simga points
  { //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //noise matrix lidar  
  MatrixXd R = MatrixXd(2,2);
  R <<    std_laspx_ * std_laspx_, 0,
           0, std_laspy_ * std_laspy_;
  S = S + R;

  //************************************************************************
  // B) UPDATE  
  //************************************************************************

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);
  Tc.fill(0.0);

  //calculate cross correlation matrix 
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // measurement
  VectorXd z = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();  

  //************************************************************************
  // CALC LIDAR NIS  
  //************************************************************************

  //double E = y.transpose() * S.inverse() * y;
  double       E = z_diff.transpose() * S.inverse() * z_diff;
  NIS_lidar_num_ += 1.0;
  if (E > NIS_LIDAR_5_) 
    NIS_lidar_num_5_ += 1.0;

  //cout << "NIS lidar: # of samples above " << NIS_LIDAR_5_ << " = " << NIS_lidar_num_5_ / NIS_lidar_num_ * 100.0 << "%" << endl; 

  return;
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

  //************************************************************************
  // A) PREDICT RADAR MEAN & COVARIANCE  
  //************************************************************************

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2*n_aug_+1);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)   //2n+1 simga points
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // check for zero
    double spxy = sqrt(p_x*p_x + p_y*p_y);
    if (spxy < 0.0001) 
      cout << "ERRORi (UpdateRadar): divison by zero" << endl;
    
    // measurement model
    Zsig(0,i) = spxy;                         //r
    Zsig(1,i) = atan2(p_y,p_x);               //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / spxy;    //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(3,3);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;


  //************************************************************************
  // B) UPDATE  
  //************************************************************************

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0.0);

  //calculate cross correlation matrix 
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // measurement
  VectorXd z = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();  

  //************************************************************************
  // CALC RADAR NIS  
  //************************************************************************

  double       E = z_diff.transpose() * S.inverse() * z_diff;
  NIS_radar_num_ += 1.0;
  if (E > NIS_RADAR_5_) 
    NIS_radar_num_5_ += 1.0;

  //cout << "NIS radar: # of samples above " << NIS_RADAR_5_ << " = " << NIS_radar_num_5_ / NIS_radar_num_ * 100.0 << "%" << endl; 

  return;
}

