clear all
close all
clc

load Eulerdata

roll = roll * (pi/180);
pitch = pitch * (pi/180);
yaw = yaw * (pi/180);


dt = 20e-3; %50hz

[m,n]=size(xacc);
t1 = ((0:(m-1))*dt)';
t_fin = (m-1)*dt;

%%get the attitude angles by integrating the angular velocity from gyroscopes
sampleperiod=dt;
Gyroscope=[rollrate pitchrate yawrate];
gyro=zeros(length(t1),3);
gyro_integral = zeros(length(t1),3);
for t = 2:length(t1)
    gyro_integral(t,:) = gyro_integral(t-1,:) + Gyroscope(t,:)*sampleperiod;
end


%%get the attitude angles by orthogonally decomposing the acceleration
Accelerometer=[xacc yacc zacc];
Magnetometer=[magx magy magz];
acc=zeros(length(t1),3);
acc_angle = zeros(length(t1),3);
for t = 1:length(t1)
    acc_angle(t,1) = atan2(Accelerometer(t,2),Accelerometer(t,3)) * (180/pi);    
    acc_angle(t,2) = asin(-Accelerometer(t,1)/norm(Accelerometer(t,:))) * (180/pi);    
    mx = Magnetometer(t,1);
    my = Magnetometer(t,2);
    mz = Magnetometer(t,3);
    at1 = Accelerometer(t,1);
    at2 = Accelerometer(t,2);
    cosr =  cos(acc_angle(t,1)*(pi/180));
    sinr =  sin(acc_angle(t,1)*(pi/180));
    cosp =  cos(acc_angle(t,2)*(pi/180));
    sinp =  sin(acc_angle(t,2)*(pi/180));    
    acc_angle(t,3) = atan2( -(my*cosr -mz*sinr),(mx*cosp+my*sinp*sinr + mz*sinp*cosr) ) * (180/pi);   
end

%%classic linear complementary filtering
tic
tao = 0.08;
CCF_roll = classical_complementary_filters('tao', tao ,'Ts',dt);
CCF_pitch = classical_complementary_filters('tao', tao ,'Ts',dt);
CCF_yaw = classical_complementary_filters('tao', tao ,'Ts',dt);
ccf_euler_roll = zeros(length(t1),1);
ccf_euler_pitch = zeros(length(t1),1);
ccf_euler_yaw = zeros(length(t1),1);
for t = 1:length(t1)
    CCF_roll.UpdateIMU( Gyroscope(t,1)* (pi/180) , acc_angle(t,1)*(pi/180) ) ;
    ccf_euler_roll(t) = CCF_roll.angle * (180/pi);
    
    CCF_pitch.UpdateIMU( Gyroscope(t,2)* (pi/180) , acc_angle(t,2)*(pi/180) ) ;
    ccf_euler_pitch(t) = CCF_pitch.angle * (180/pi);
    
    CCF_yaw.UpdateIMU( Gyroscope(t,3)* (pi/180) , acc_angle(t,3)*(pi/180) ) ;
    ccf_euler_yaw(t) = CCF_yaw.angle * (180/pi);
end
toc
figure('Name', 'classic linear complementary filtering');
subplot(311)
plot(t1,ccf_euler_roll,t1,roll*(180/pi),'r');title('roll comparison');xlabel('Time(s)');grid;
legend('ccf','ground truth')
subplot(312)
plot(t1,ccf_euler_pitch,t1,pitch*(180/pi),'r');title('pitch comparison');xlabel('Time(s)');grid;
legend('ccf','ground truth')
subplot(313)
plot(t1,ccf_euler_yaw-mean(ccf_euler_yaw),t1,yaw*(180/pi),'r');title('yaw comparison');xlabel('Time(s)');grid;
legend('ccf','ground truth')

mean(abs(ccf_euler_roll*pi/180-roll))
mean(abs(ccf_euler_pitch*pi/180-pitch))
mean(abs((ccf_euler_yaw-mean(ccf_euler_yaw))*pi/180-yaw))