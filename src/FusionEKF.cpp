#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

#define EPS 0.0001 // A very small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    
    MatrixXd F_ = MatrixXd(4, 4);
    F_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    
    MatrixXd Q_ = MatrixXd(4, 4);
    
    VectorXd x_ = VectorXd(4);
    
    H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;
    
    MatrixXd P_ = MatrixXd(4, 4);
    VectorXd x_in = VectorXd(4);
    
    ekf_.Init(x_in, P_, F_, H_laser_, R_laser_, Q_);
    
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


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
    cout << "EKF: " << endl;
//    ekf_.x_ = VectorXd(4);
//    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        cout << "Initializing by RADAR data\n";
        double rho = measurement_pack.raw_measurements_(0);
        double phi = measurement_pack.raw_measurements_(1);
        double rhod = measurement_pack.raw_measurements_(2);

        double sigrho = R_radar_(0, 0);
        double sigphi = R_radar_(1, 1);
        double sigrhod = R_radar_(2, 2);
        
        double cosp = cos(phi);
        double sinp = sin(phi);
        
        double px = rho * cosp;
        double py = rho * sinp;
        
        // improve vx, vy approcximations. currently underestimated as perpendicular component is ignored
        double vx = 0.0;//rhod * cosp;
        double vy = 0.0;//rhod * sinp;
        
        double sigvx = sigrhod * sigphi * sigphi /2 + 1000;
        double sigvy = sigrhod * sigphi * sigphi /2 + 1000;
        double sigpx = sigrho * sigphi * sigphi /2;
        double sigpy = sigrho * sigphi * sigphi /2;
        
        ekf_.x_ << px, py, vx, vy;
//        ekf_.P_ << sigpx, 0, 0, 0,
//                    0, sigpy, 0, 0,
//                    0, 0, sigvx, 0,
//                    0, 0, 0, sigvy;
        
//        cout << "INITIAL x_ = " << ekf_.x_ << endl;
//        cout << "INTIAL P_ = " << ekf_.P_ << endl;
        
        
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        cout << "Initializing by LASER data\n";
        ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
//        ekf_.P_ << R_laser_(0, 0), R_laser_(0, 1), 0, 0,
//                    R_laser_(1, 0), R_laser_(1, 1), 0, 0,
//                    0, 0, 1000, 0,
//                    0, 0, 0, 1000;
//        cout << "INITIAL x_ = " << ekf_.x_ << endl;
//        cout << "INTIAL P_ = " << ekf_.P_ << endl;
    }
      ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;
      
      
//    // TODO Edit
//      if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
//          ekf_.x_(0) = EPS;
//          ekf_.x_(1) = EPS;
//      }
//      
      
      
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
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
    static double noise_ax = 9.0;
    static double noise_ay = 9.0;
    
    double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    if( dt > 0. ) {
        
        
        // NEW CODE //
        float dt2 = dt*dt;
        float dt3 = dt2*dt;
        float dt4 = dt3*dt;
        float dt4over4 = dt4/4.;
        float dt3over2 = dt3/2.;
        ekf_.Q_ << dt4over4*noise_ax,                 0, dt3over2*noise_ax,                0,
			     0, dt4over4*noise_ay,                 0, dt3over2*noise_ay,
        dt3over2*noise_ax,                 0,      dt2*noise_ax,                 0,
			     0, dt3over2*noise_ay,                 0,      dt2*noise_ay;
        
    
    // setting F matrix
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    ekf_.Predict();
    cout << "POST PREDICT x_ = " << ekf_.x_ << endl;
    cout << "POST PREDICT P_ = " << ekf_.P_ << endl;
}

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
      cout << "\n\n RADAR \n\n";
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
      cout << "\n\n LASER \n\n";
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl << "\n\n ENDS \n\n";
}
