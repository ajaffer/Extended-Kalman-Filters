#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //Process noise
  Q_ = MatrixXd(4, 4);
  
  //state covariance matrix P
  P_ = MatrixXd(4, 4); 
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  // //measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // cout << "** ProcessMeasurement" << endl;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    // cout << "first measurement " << endl;
    VectorXd x_ = VectorXd(4);
    // ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // cout << "Radar" << endl; 
      VectorXd raw = tools.PolarToCartesian(measurement_pack.raw_measurements_);
      
      //set the state with the initial location and zero velocity
      x_ << raw[0], raw[1], 0, 0;
      ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state, set the state with the initial location and zero velocity
      */
      // cout << "Laser" << endl; 
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
    }    
    
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    // cout << "Init Done" << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float noise_ax = 9;
  float noise_ay = 9;

  //compute the time elapsed between the current and previous measurements
  //- Time is measured in seconds.
  // cout << "measurement_pack.timestamp_ " << measurement_pack.timestamp_ << endl;
  // cout << "(measurement_pack.timestamp_ - previous_timestamp_)" << (measurement_pack.timestamp_ - previous_timestamp_) << endl;
  float dt = ((measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0);	//dt - expressed in seconds
  // cout << "dt: " << dt << endl;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // cout << "F prev: " << ekf_.F_ << endl;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // cout << "F: " << ekf_.F_ << endl;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
  dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
  0, dt_3/2*noise_ay, 0, dt_2*noise_ay;   

  // cout << "Q: " << ekf_.Q_ << endl;

  ekf_.Predict();
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // VectorXd rawCartesian = tools.PolarToCartesian(measurement_pack.raw_measurements_);
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
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
  // cout << "|| ProcessMeasurement" << endl;
}
