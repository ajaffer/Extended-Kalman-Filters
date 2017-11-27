#include <iostream>
#include <math.h>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
                                // // cout << "tools-0" << endl;
                                // VectorXd RMSE(4);
                                // // cout << "tools-1" << endl;
                                // RMSE << 1,1,1,1;
                                // // cout << "tools-2" << endl;
                                // return RMSE;
		// cout << "** CalculateRMSE" << endl;

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	// cout << "rmse: " << rmse << endl;
	// cout << "|| CalculateRMSE" << endl;

	return rmse;  
}

VectorXd Tools::PolarToCartesian(const VectorXd& raw_measurements_) {
	VectorXd x(2);
	// cout << "** PolarToCartesian()" << endl;

  //recover state parameters
	float rho = raw_measurements_(0);
	float phi = raw_measurements_(1);
	float rhoDot = raw_measurements_(2);

	// cout << "polar (rho, phi, rhoDot) (" << rho << 
	// 	", " << phi << 
	// 	", " << rhoDot << ")" << endl;

	//compute the Cartesian vector
	x << rho * cos(phi),
			 rho * sin(phi);
	
	// cout << "cartesian (x,y) (" << x[0] << "," << x[1] << ")" << endl;		 

	// cout << "|| PolarToCartesian()" << endl;
	
	return x;  
}

VectorXd Tools::CartesianToPolar(const VectorXd& x_state) {
	// cout << "** CartesianToPolar()" << endl;
	
	VectorXd h(3);

  //recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = sqrt(px*px+py*py);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CartesianToPolar () - Error - Division by Zero" << endl;
		return h;
	}

	//compute the Jacobian matrix
	float phi = atan2(py, px);
	h << c1,
       phi,
       ((px*vx+py*vy) / c1);

	// cout << "h: " << h << endl;		 
  // cout << "|| CartesianToPolar()" << endl;
			 
	return h;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	// cout << "** CalculateJacobian()" << endl;
	// cout << "x_state: " << x_state << endl;

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	// cout << "|| CalculateJacobian()" << endl;
			
	return Hj;  
}
