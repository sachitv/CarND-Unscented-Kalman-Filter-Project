#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include "measurement_package.h"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd( 5 );

	// initial covariance matrix
	P_ = MatrixXd( 5, 5 );

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 2; //originally 30

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = M_PI / 2.0; //originally 30

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	// Sigma points
	Xsig_pred_ = MatrixXd( 5, 2 * n_aug_ + 1 );

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/

	// State Dimensions
	n_x_ = 5;

	// Augmented State Dimensions
	n_aug_ = 7;

	//Sigma point spreading factor
	lambda_ = 3 - n_aug_;

	// weights_
	weights_ = VectorXd( 2 * n_aug_ + 1 );

	//measurement covariance matrix - laser
	R_Laser = MatrixXd(2, 2);
	R_Laser <<    	0.0225, 0,
					0, 0.0225;

	m_H = MatrixXd(2, 5);
	m_H << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0;


	is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement( MeasurementPackage meas_package )
{
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/

	if ( !is_initialized_ )
	{
		//Let's initialize the covariance matrix
		P_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

		if ( meas_package.sensor_type_ == MeasurementPackage::RADAR )
		{
			if ( use_radar_ )
			{
				double const rho = meas_package.raw_measurements_[ 0 ];
				double const phi = meas_package.raw_measurements_[ 1 ];
				double const rhodot = meas_package.raw_measurements_[ 2 ];

				double const px = rho * cos( phi );
				double const py = rho * sin( phi );
				double const v = rhodot;

				x_.fill( 0 );
				x_ << px, py, v, 0, 0;

				is_initialized_ = true;
			}
		}
		else
		{
			if ( use_laser_ )
			{
				x_.fill( 0 );
				x_( 0 ) = meas_package.raw_measurements_( 0 );
				x_( 1 ) = meas_package.raw_measurements_( 1 );

				is_initialized_ = true;
			}

		}

		time_us_ = meas_package.timestamp_;

	}
	else
	{
		double const deltaTime = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;

		if(deltaTime > 0.f)
		{
			Prediction( deltaTime );
		}

		if ( meas_package.sensor_type_ == MeasurementPackage::RADAR )
		{
			if ( use_radar_ )
			{
				UpdateRadar( meas_package );
			}
		}
		else
		{
			if ( use_laser_ )
			{
				UpdateLidar( meas_package );
			}
		}
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction( double delta_t )
{
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	//create augmented mean vector
	VectorXd x_aug = VectorXd( 7 );
	x_aug.fill( 0 );

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd( 7, 7 );
	P_aug.fill( 0 );

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd( n_aug_, 2 * n_aug_ + 1 );
	Xsig_aug.fill( 0 );

	//create augmented mean state
	x_aug.head( 5 ) = x_;
	//create augmented covariance matrix
	P_aug.topLeftCorner( 5, 5 ) = P_;
	P_aug( 5, 5 ) = std_a_ * std_a_;
	P_aug( 6, 6 ) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();
	MatrixXd sqrtRt = sqrt( lambda_ + n_aug_ ) * A;

	//create augmented sigma points
	Xsig_aug.col( 0 ) = x_aug;

	//set remaining sigma points
	for ( int i = 0; i < n_aug_; i++ )
	{
		Xsig_aug.col( i + 1 ) = x_aug + sqrtRt.col( i );
		Xsig_aug.col( i + 1 + n_aug_ ) = x_aug - sqrtRt.col( i );
	}

	//predict sigma points
	Xsig_pred_.fill(0);
	for ( int i = 0; i < 2 * n_aug_ + 1; ++i )
	{
		Xsig_pred_( 0, i ) = Xsig_aug( 0, i );
		Xsig_pred_( 1, i ) = Xsig_aug( 1, i );
		Xsig_pred_( 2, i ) = Xsig_aug( 2, i );
		Xsig_pred_( 3, i ) = Xsig_aug( 3, i );
		Xsig_pred_( 4, i ) = Xsig_aug( 4, i );

		double const vk = Xsig_pred_( 2, i );
		double const yawK = Xsig_pred_( 3, i );
		double const yawKdot = Xsig_pred_( 4, i );

		double const nu_a = Xsig_aug( 5, i );
		double const yawdotdot = Xsig_aug( 6, i );

		//Generate second term (noise)
		VectorXd secondTerm( 5 );
		secondTerm( 0 ) = 0.5 * delta_t * delta_t * cos( yawK ) * nu_a;
		secondTerm( 1 ) = 0.5 * delta_t * delta_t * sin( yawK ) * nu_a;
		secondTerm( 2 ) = delta_t * nu_a;
		secondTerm( 3 ) = 0.5 * delta_t * delta_t * yawdotdot;
		secondTerm( 4 ) = delta_t * yawdotdot;

		//avoid division by zero and generate first term
		VectorXd firstTerm( 5 );
		if ( fabs( yawKdot ) > 0.001 )
		{
			firstTerm( 0 ) = vk / yawKdot * (sin( yawK + yawKdot * delta_t ) - sin( yawK ));
			firstTerm( 1 ) = vk / yawKdot * (-cos( yawK + yawKdot * delta_t ) + cos( yawK ));
			firstTerm( 2 ) = 0;
			firstTerm( 3 ) = yawKdot * delta_t;
			firstTerm( 4 ) = 0;
		}
		else
		{
			firstTerm( 0 ) = vk * cos( yawK ) * delta_t;
			firstTerm( 1 ) = vk * sin( yawK ) * delta_t;
			firstTerm( 2 ) = 0;
			firstTerm( 3 ) = yawKdot * delta_t;
			firstTerm( 4 ) = 0;
		}

		//write predicted sigma points into right column
		Xsig_pred_( 0, i ) += firstTerm( 0 ) + secondTerm( 0 );
		Xsig_pred_( 1, i ) += firstTerm( 1 ) + secondTerm( 1 );
		Xsig_pred_( 2, i ) += firstTerm( 2 ) + secondTerm( 2 );
		Xsig_pred_( 3, i ) += firstTerm( 3 ) + secondTerm( 3 );
		Xsig_pred_( 4, i ) += firstTerm( 4 ) + secondTerm( 4 );
	}

	//DEBUG
	std::cout << "Xsig_pred_" << std::endl;
	for ( int i = 0; i < 5; ++i )
	{
		for ( int j = 0; j < 2 * n_aug_ + 1; ++j )
		{
			std::cout << Xsig_pred_( i, j ) << "\t\t\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
	//END DEBUG

	//set weights_
	weights_( 0 ) = lambda_ / (lambda_ + n_aug_);

	for ( int i = 1; i < 2 * n_aug_ + 1; ++i )
	{
		weights_( i ) = 1 / (2 * (lambda_ + n_aug_));
	}

	//predict state mean
	x_.fill( 0 );
	for ( int i = 0; i < 2 * n_aug_ + 1; ++i )
	{
		x_ += ( weights_( i ) * Xsig_pred_.col( i ) );
	}

	//DEBUG
	std::cout << "predicted state mean\n" << x_ << std::endl;
	std::cout << std::endl << std::endl;
	//END DEBUG

	P_.fill( 0 );
	//predict state covariance matrix
	for ( int i = 0; i < 2 * n_aug_ + 1; ++i )
	{
		// state difference
		VectorXd x_diff = Xsig_pred_.col( i ) - x_;

		//angle normalization
		x_diff( 3 ) = remainder( x_diff( 3 ), 2.0 * M_PI );

		P_ = P_ + weights_( i ) * x_diff * x_diff.transpose();
	}

	//DEBUG
	std::cout << "predicted state covariance\n" << P_ << std::endl;
	std::cout << std::endl << std::endl;
	//END DEBUG
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar( MeasurementPackage meas_package )
{
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/

	VectorXd z = meas_package.raw_measurements_;

	MatrixXd const y = z - m_H * x_;

	MatrixXd const HT(m_H.transpose());

	MatrixXd const S = m_H * P_ * HT + R_Laser;
	MatrixXd const K = P_ * HT * S.inverse();


	x_ = x_ + K * y;

	long const x_size = x_.size();
	MatrixXd const I(MatrixXd::Identity(x_size, x_size));

	P_ = (I - K * m_H) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar( MeasurementPackage meas_package )
{
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/

	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd( n_z, 2 * n_aug_ + 1 );
	Zsig.fill(0);

	//transform sigma points into measurement space
	for ( int i = 0; i < 2 * n_aug_ + 1; i++ )
	{  //2n+1 sigma points

		// extract values for better readability
		double p_x = Xsig_pred_( 0, i );
		double p_y = Xsig_pred_( 1, i );
		double v = Xsig_pred_( 2, i );
		double yaw = Xsig_pred_( 3, i );

		double v1 = cos( yaw ) * v;
		double v2 = sin( yaw ) * v;

		// measurement model
		Zsig( 0, i ) = sqrt( p_x * p_x + p_y * p_y );                        		//r
		Zsig( 1, i ) = atan2( p_y, p_x );                                 			//phi

		//Avoid divide by zero
		if ( Zsig( 0, i ) < 0.01f )
		{
			Zsig( 2, i ) = 0;                             							//r_dot
		}
		else
		{
			Zsig( 2, i ) = (p_x * v1 + p_y * v2) / sqrt( p_x * p_x + p_y * p_y );   //r_dot
		}

	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd( n_z );
	z_pred.fill( 0.0 );
	for ( int i = 0; i < 2 * n_aug_ + 1; i++ )
	{
		z_pred = z_pred + weights_( i ) * Zsig.col( i );
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd( n_z, n_z );
	S.fill( 0.0 );
	for ( int i = 0; i < 2 * n_aug_ + 1; i++ )
	{  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col( i ) - z_pred;

		//angle normalization
		z_diff( 1 ) = remainder( z_diff( 1 ), 2.0 * M_PI );

		S = S + weights_( i ) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd( n_z, n_z );
	R << std_radr_ * std_radr_, 0, 0,
			0, std_radphi_ * std_radphi_, 0,
			0, 0, std_radrd_ * std_radrd_;
	S = S + R;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd( n_x_, n_z );

	//calculate cross correlation matrix
	Tc.fill( 0 );
	for ( int i = 0; i < 2 * n_aug_ + 1; ++i )
	{
		VectorXd z_diff = Zsig.col( i ) - z_pred;

		//Normalization
		z_diff( 1 ) = remainder( z_diff( 1 ), 2.0 * M_PI );

		// state difference
		VectorXd x_diff = Xsig_pred_.col( i ) - x_;

		// Normalization
		x_diff( 3 ) = remainder( x_diff( 3 ), 2.0 * M_PI );

		Tc = Tc + weights_( i ) * x_diff * z_diff.transpose();
	}

	//DEBUG OUTPUT
	std::cout<<"Cross Correlation Matrix : "<<std::endl;
	std::cout<<Tc<<std::endl<<std::endl;
	//DEBUG OUTPUT

	//calculate Kalman gain K;
	MatrixXd const K = Tc * S.inverse();

	VectorXd const z = meas_package.raw_measurements_;

	//update state mean and covariance matrix
	VectorXd z_diff = (z - z_pred);

	// Normalize
	z_diff( 1 ) = remainder( z_diff( 1 ), 2.0 * M_PI );

	x_ = x_ + K * z_diff;
	P_ = P_ - (K * S * K.transpose());

	//Calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

	//Debug
	std::cout << "X" << std::endl;
	for ( int i = 0; i < 5; ++i )
	{
		std::cout << x_( i ) << std::endl;
	}
	std::cout << std::endl << std::endl;
	// Debug

	//Debug
	std::cout << "P" << std::endl;
	std::cout<<P_<<std::endl;
	std::cout << std::endl << std::endl;
	// Debug
}
