#include <ros/ros.h>
#include "geometry_msgs/Twist.h"
#include "std_msgs/Float32MultiArray.h"
#include <math.h>
#include <armadillo>
#include <iostream>

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
double I = Ic + mc * pow(d,2) + 2 * mw * pow(l,2) + 2 * Im;

double k1 = 0.22;
double K2 = 5.5;
double B1 = -4.2;
double k3 = 0.22;
double K4 = 5.5;
double B2 = -4.2;
double K1;
double K3;

mat Ak(4,4,fill::zeros);
mat Bk(4,2,fill::zeros);
mat Ck(4,4,fill::eye);
mat Dk(2,2,fill::zeros);

mat Nk(4,4,fill::zeros);
mat Qk(4,4,fill::zeros);
mat Rk(4,4,fill::zeros);

vec x(4,fill::zeros);
vec xp_hat(4,fill::zeros);
vec xpp_hat(4,fill::zeros);
mat P(4,4,fill::zeros);
mat Pp(4,4,fill::zeros);
mat Ppp(4,4,fill::zeros);

mat Kk(4,4,fill::zeros);
vec U(2,fill::zeros);
vec Y(4,fill::zeros);
vec ytilde(4,fill::zeros);

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

	K1 = (k1 - 1)/LOOP_RATE;
	K3 = (k3 - 1)/LOOP_RATE;
	
	Bk(2,0) = B1;
	Bk(3,1) = B2;

	Nk(0,0) = 100.0;
	Nk(1,1) = 100.0;
	Nk(2,2) = 0.01;
	Nk(3,3) = 0.01;

	Qk(0,0) = 0.001;
	Qk(1,1) = 0.001;
	Qk(2,2) = 1.0;
	Qk(3,3) = 1.0;

	Rk(0,0) = 0.0;
	Rk(1,1) = 0.0;
	Rk(2,2) = 0.0;
	Rk(3,3) = 0.0;

}


void motorPowerCallback(const geometry_msgs::Twist::ConstPtr& msg)
{
	U(1) = msg->linear.x; // u_l
	U(0) = msg->linear.y; // u_r
}

void vwCallback(const std_msgs::Float32MultiArray::ConstPtr& msg)
{
	x(0) = msg->data[0]; // v
	x(1) = msg->data[1]; // w
}

void systemSetup()
{
//	Ak = {{0.0, (-2 * d * mc * x(1))/(m + (2 * Iw)/pow(r,2)), 1/(r * (m + (2 * Iw)/pow(r,2))),  1/(r * (m + (2 * Iw)/pow(r,2)))},
//		  {(d * mc * x(1))/(I + (2 * Iw * pow(l,2))/pow(r,2)), (d * mc * x(0))/(I + (2 * Iw * pow(l,2))/pow(r,2)), l/(r * (I + (2 * Iw * pow(l,2))/pow(r,2))), -l/(r * (I + (2 * Iw * pow(l,2))/pow(r,2)))},
//		  {K2, K2 * l, K1 , 0},
//		  {K4, -K4 * l, 0 , K3}};

	Ak(0,0) = 0.0;
	Ak(0,1) = (-2 * d * mc * x(1))/(m + (2 * Iw)/pow(r,2));
	Ak(0,2) = 1/(r * (m + (2 * Iw)/pow(r,2)));
	Ak(0,3) = 1/(r * (m + (2 * Iw)/pow(r,2)));
	
	Ak(1,0) = (d * mc * x(1))/(I + (2 * Iw * pow(l,2))/pow(r,2));
	Ak(1,1) = (d * mc * x(0))/(I + (2 * Iw * pow(l,2))/pow(r,2));
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
	
	xp_hat = xpp_hat;
	Pp = Ppp;
}

void publich()
{
	vw_hat.linear.x = xp_hat(0);
	vw_hat.angular.z = xp_hat(1);

	torque_hat.linear.x = xp_hat(3);
	torque_hat.linear.y = xp_hat(2);

	vwHat.publish(vw_hat);
	torque.publish(torque_hat);
}

int main(int argc, char** argv)
{
	ros::init(argc, argv, "kalmanfilter");
	ros::NodeHandle nh;
	ros::NodeHandle nh_("~");
	
	nh_.param<int>("loop_rate", LOOP_RATE, 20);

	vwHat = nh.advertise<geometry_msgs::Twist>("vw_hat", 5);
	torque = nh.advertise<geometry_msgs::Twist>("torque", 5);
	motor_power = nh.subscribe<geometry_msgs::Twist>("motor_power", 1, motorPowerCallback);
	vw = nh.subscribe<std_msgs::Float32MultiArray>("vw_estimate", 1, vwCallback);
	
	ros::Rate loop_rate(LOOP_RATE);
	init();
	while(ros::ok())
	{
		ros::spinOnce();
		publich();
		loop_rate.sleep();
	}
	
	return 0;
}

