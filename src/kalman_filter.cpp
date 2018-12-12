#include <ros/ros.h>
#include "geometry_msgs/Twist.h"
#include "std_msgs/Float32MultiArray.h"
#include <armadillo>

using namespace arma;

ros::Publisher vwHat;
ros::Publisher torque;

ros::Subscriber vw;
ros::Subscriber motor_power;

geometry_msgs::Twist torque_hat;
geometry_msgs::Twist vw_hat;

int LOOP_RATE;

double mc = 100.0;
double mw = 7.5;
double m = mc + 2* mw;
double r = 0.165;
double l = 0.285;
double d = 0.25;
double Ic = 0.5 * mc * pow(d, 2);
double Im = mw * pow(r,2);
double Iw = mw * pow(l,2);
double I = (Ic + mc * pow(d,2) + 2 * mw * pow(l,2) + 2 * Im);

double K1 = -32.3;
double K2 = -200.8;
double B1 = 195.2;
double K3 = -32.3;
double K4 = -200.8;
double B2 = 195.2;
//double K1 = -34.3;
//double K2 = -246.8;
//double B1 = 180.2;
//double K3 = -30.3;
//double K4 = -163.1;
//double B2 = 213.5;
//double K1;
//double K3;

mat Ak(4,4,fill::zeros);
mat Bk(4,2,fill::zeros);
mat Ck(2,4,fill::zeros);
mat Dk(2,2,fill::zeros);

mat Nk(4,4,fill::zeros);
mat Qk(4,4,fill::zeros);
mat Rk(2,2,fill::zeros);

//vec x(4,fill::zeros);
vec xp_hat(4,fill::zeros);
vec xpp_hat(4,fill::zeros);
mat P(4,4,fill::zeros);
mat Pp(4,4,fill::zeros);
mat Ppp(4,4,fill::zeros);

mat Kk(4,2,fill::zeros);
vec U(2,fill::zeros);
vec Y(2,fill::zeros);
vec ytilde(4,fill::zeros);


// set upp static matrixes
void init()
{
	torque_hat.linear.x = 0;	
	torque_hat.linear.y = 0;	
	torque_hat.linear.z = 0;	
	
	torque_hat.angular.x = 0;	
	torque_hat.angular.y = 0;	
	torque_hat.angular.z = 0;	
	
	vw_hat.linear.x = 0;
	vw_hat.linear.y = 0;
	vw_hat.linear.z = 0;
	
	vw_hat.angular.x = 0;
	vw_hat.angular.y = 0;
	vw_hat.angular.z = 0;

	//K1 = (k1 - 1)/LOOP_RATE;
	//K3 = (k3 - 1)/LOOP_RATE;
	
	Bk(2,0) = B1/LOOP_RATE;
	Bk(3,1) = B2/LOOP_RATE;

	Ck(0,0) = 1;
	Ck(1,1) = 1;

	Nk(0,0) = 1.0;
	Nk(1,1) = 1.0;
	Nk(2,2) = 1.0;
	Nk(3,3) = 1.0;
	Nk = 10 * Nk;	

	Qk(0,0) = 1.0;
	Qk(1,1) = 1.0;
	Qk(2,2) = 1.0;
	Qk(3,3) = 1.0;

	Rk(0,0) = 100.0;
	Rk(1,1) = 100.0;
	Rk(0,1) = 0.0;
	Rk(1,0) = 0.0;

}


void motorPowerCallback(const geometry_msgs::Twist::ConstPtr& msg)
{
	U(1) = msg->linear.x / 100; // u_l normolised
	U(0) = msg->linear.y / 100; // u_r
}

void vwCallback(const geometry_msgs::Twist::ConstPtr& msg)
{
	Y(0) = msg->linear.x; // v
	Y(1) = msg->angular.z; // w
	// ROS_INFO("vwCalback v in: %f", x(0));
}

void systemSetup()
{
	// config for Ak matrix
	Ak(0,0) = 0.0;
	Ak(0,1) = (-2 * d * mc * Y(1))/(m + (2 * Iw)/pow(r,2));
	Ak(0,2) = 1/(r * (m + (2 * Iw)/pow(r,2)));
	Ak(0,3) = 1/(r * (m + (2 * Iw)/pow(r,2)));
	
	Ak(1,0) = (d * mc * Y(1))/(I + (2 * Iw * pow(l,2))/pow(r,2));
	Ak(1,1) = (d * mc * Y(0))/(I + (2 * Iw * pow(l,2))/pow(r,2));
	Ak(1,2) = l/(r * (I + (2 * Iw * pow(l,2))/pow(r,2)));
	Ak(1,3) = -l/(r * (I + (2 * Iw * pow(l,2))/pow(r,2)));
	
	Ak(2,0) = K2;
	Ak(2,1) = K2 * l;
	Ak(2,2) = K1;
	Ak(2,3) = 0.0;
	
	Ak(3,0) = K4;
	Ak(3,1) = -K4 * l;
	Ak(3,2) = 0.0;
	Ak(3,3) = K3;
	
	mat I(4,4,fill::eye);
	Ak = (Ak/LOOP_RATE) + I;
	// Ak config end
	
}

void prediction()
{
	xp_hat = Ak * xp_hat + Bk * U;
	Pp = Ak * P * trans(Ak) + Nk * Qk * trans(Nk);
}

void update()
{
	Kk = Pp * trans(Ck) * inv(Ck * Pp * trans(Ck) + Rk);
	ytilde = Y - Ck * xp_hat;
	xpp_hat = xp_hat + Kk * ytilde;
	Ppp = Pp - Kk * Ck * Pp;
	

	// set for next iteration
	xp_hat = xpp_hat;
	P = Ppp;
}

void publich()
{
	vw_hat.linear.x = (float) xp_hat(0); // v
	vw_hat.angular.z = (float) xp_hat(1); // w

	torque_hat.linear.x = (float) xp_hat(3); // tl
	torque_hat.linear.y = (float) xp_hat(2); // tr

	vwHat.publish(vw_hat);
	torque.publish(torque_hat);
}

// looping trigerd at desierd frequens
void loopCallback(const ros::TimerEvent&)
{
	systemSetup();		// kalman filter, systemSetup, prediction, update
	prediction();
	update();
	publich();
}


int main(int argc, char** argv)
{
	ros::init(argc, argv, "kalmanfilter");
	ros::NodeHandle nh;
	ros::NodeHandle nh_("~");
	
	nh_.param<int>("loop_rate", LOOP_RATE, 100);

	vwHat = nh.advertise<geometry_msgs::Twist>("vw_hat", 5);
	torque = nh.advertise<geometry_msgs::Twist>("torque", 5);
	motor_power = nh.subscribe<geometry_msgs::Twist>("motor_power", 1, motorPowerCallback);
	vw = nh.subscribe<geometry_msgs::Twist>("vw_estimate", 1, vwCallback);

	ros::Timer loop = nh.createTimer(ros::Duration(0.01), loopCallback);	
	
	init();
	ros::spin();

	return 0;
}

