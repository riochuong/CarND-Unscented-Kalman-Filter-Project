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
    * Calculate the RMSE here.
  */
  VectorXd results(estimations[0].size());
  for (int i = 0; i < results.size(); i++) {
     results(i) = 0;
  }
  for (unsigned int i = 0; i < estimations.size(); i++){
     VectorXd est = estimations[i];
     VectorXd gt = ground_truth[i];
     results = results + (est - gt).cwiseProduct(est - gt) ;
  }
  results = results / estimations.size();

   for (int i = 0; i < results.size(); i++) {
     results(i) = sqrt(results(i));
  }
  return results;

}
