#include "kalman_filter.h"
#include <iostream>

#define PI 3.14159265
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
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
    
}

//
//inline void calculate_(MatrixXd ,){
//    
//}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd K = P_ * Ht * S.inverse();
    
    x_ = x_ + K * y;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    
    
    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    
    float rho = sqrt(px*px + py*py);
    if( rho == 0.)
        return;

    float phi = 0;
    float rho_dot = 0;
    
    // avoid division by zero
    if(fabs(px) < 0.0001){
        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << std::endl;
    }
    else {
        phi = atan2(py, px);
    }
    
    // avoid division by zero
    if (rho < 0.0001) {
        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
    }
    else {
        rho_dot = (px*vx + py*vy) / rho;
    }
    
    VectorXd h = VectorXd(3);
    h << rho, phi, rho_dot; // For radar H * x becomes h(x)
    
    VectorXd y = z - h; // Using h instead of Jacobian Hj_ here!
    
//    while (y(1)>PI) {
//        y(1) -= 2 * PI;
//    }
//    while (y(1)<-PI) {
//        y(1) += 2 * PI;
//    }
    
    if (y(1) > M_PI)
        y(1) = fmod((y(1) - M_PI), 2*M_PI) - M_PI;
    if (y(1) < -M_PI)
        y(1) = fmod((y(1) + M_PI), 2*M_PI) + M_PI;
    
    /* TODO check
     If angle > PI then
     angle = ((angle - Pi) % 2*Pi) - Pi
     
     If angle < -PI then
     angle = ((angle + Pi) % 2*Pi) + Pi
     */
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd K = P_ * Ht * S.inverse();
    
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * H_) * P_;
}
