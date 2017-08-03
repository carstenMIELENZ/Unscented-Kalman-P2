#include <iostream>
#include "tools.h"

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
  //-------------------------------------------------------------------------------
  //
  // code reuse from Extended kalman-filter 
  //
  //-------------------------------------------------------------------------------

  // initialize return vector
  VectorXd rmse(4);
           rmse << 0,0,0,0;

  // perform checks

  if (estimations.size() == 0 || (estimations.size() > ground_truth.size()) ) {

    cout << "ERROR (CalculateRMSE): vector size not matching or zero, returning zero" << endl;
    return rmse;

  }

  // calculate
  // 1. residual

  for (int i=0; i < estimations.size(); ++i) {

     VectorXd buf = estimations[i] - ground_truth[i];
              buf = buf.array() * buf.array();
              // update residual
              rmse += buf;
  }

  // 2. mean residual
  rmse = rmse / estimations.size();
  // 3. square root (wurzel)
  rmse = rmse.array().sqrt();

  //std::cout << "RMSE" << std::endl << rmse << std::endl; 

  return rmse;
  
}
