#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
	
	// cout << "KF Initialized " << endl;
	// cout << "x = " << x_ << endl;
	// cout << "P = " << P_ << endl;
	// cout << "F = " << F_ << endl;
	// cout << "H = " << H_ << endl;
	// cout << "R = " << R_ << endl;
	// cout << "Q = " << Q_ << endl;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
	*/
	// cout << "** Predict" << endl;
	x_ = F_ * x_;
	// cout << "x: " << x_ << endl;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
	// cout << "P: " << P_ << endl;
	// cout << "|| Predict" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	// cout << "** Update()"  << endl;
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;  

	// cout << "x: " << x_ << endl;
	// cout << "P: " << P_ << endl;
	// cout << "|| Update()"  << endl;
}

const double pi = 3.14159265358979323846;
const double two_pi = 6.28318530718;

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

	// cout << "** UpdateEKF()"  << endl;
	
	VectorXd z_pred = tools.CartesianToPolar(x_);
	// cout << "z: " << z << endl;
	VectorXd y = z - z_pred;
	
	double phi = y[1];

	if (phi < (-1 * pi)) {
		// cout << "[warning] phi was less than pi: " << phi << endl;
		phi += two_pi;
	} else if (phi > (pi)) {
		// cout << "[warning] phi was greater than pi: " << phi << endl;
		phi -= two_pi;
	}
	y[1] = phi;

	  

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;    

  // MatrixXd Hj = Tools::CalculateJacobian(x_);

  // VectorXd h = Tools:CartesianToPolar(x_);
	// VectorXd y = z - h;
	// MatrixXd Ht = Hj.transpose();
	// MatrixXd S = Hj * P_ * Ht + R_;
	// MatrixXd Si = S.inverse();
	// MatrixXd PHt = P_ * Ht;
	// MatrixXd K = PHt * Si;

	// //new estimate
	// x_ = x_ + (K * y);
	// long x_size = x_.size();
	// MatrixXd I = MatrixXd::Identity(x_size, x_size);
	// P_ = (I - K * Hj) * P_;  

	// cout << "|| UpdateEKF()"  << endl;
}
