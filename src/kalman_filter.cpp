#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
	// cout << "** Predict" << endl;
	x_ = F_ * x_;
	// cout << "x: " << x_ << endl;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
	// cout << "P: " << P_ << endl;
	// cout << "|| Predict" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
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

void KalmanFilter::UpdateEKF(const VectorXd &z) {
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

	// cout << "|| UpdateEKF()"  << endl;
}
