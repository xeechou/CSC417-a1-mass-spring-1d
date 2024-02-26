#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template <typename FORCE, typename STIFFNESS> inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force, STIFFNESS &stiffness) {
    Eigen::MatrixXd K;

    stiffness(K, q, qdot);
    //This is crucial, don't know why they flip it.
    double k = -K(0, 0);

    //the long equation for coeff.
    double coeff = (1.0 + dt * dt * k / mass);

    //the formula to solve the implicit Euler update
    Eigen::VectorXd right = dt * qdot + q;
    Eigen::VectorXd q_old = q;
    q                     = (dt * qdot + q_old) / coeff;
    qdot                  = (qdot - (dt * k / mass) * q_old) / coeff;
}
