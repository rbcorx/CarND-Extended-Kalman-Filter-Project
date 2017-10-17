#include <iostream>
#include "tools.h"

#define EPS 0.0001 // A very small number
#define EPS2 0.0000001

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse = VectorXd(4);
    rmse << 0,0,0,0;
    
    if (!estimations.size() || (estimations.size() != ground_truth.size()))
        return rmse;
    
    for (int i=0; i<estimations.size(); ++i){
        VectorXd temp = estimations[i] - ground_truth[i];
        temp = temp.array() * temp.array();
        rmse += temp;
    }
    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */     

    
    MatrixXd Hj(3,4);
    
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //check division by zero
    if( px == 0 && py == 0 )
    {
        cout << "Error:  division by zero in CalculateJacobian" << endl;
        return Hj;
    }
    
    
    float rho = sqrt( px*px + py* py );
    float rho2 = rho*rho;
    float rho3 = rho2*rho;
    Hj <<                 px/rho,                    py/rho,      0,      0,
    -py/rho2,                   px/rho2,      0,      0,
    py*( vx*py - vy*px )/rho3, px*( vy*px - vx*py )/rho3, px/rho, py/rho;
    
    return Hj;
}
