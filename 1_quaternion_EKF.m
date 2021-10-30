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

gx = rollrate;
gy = pitchrate;
gz = yawrate;
ax = xacc;
ay = yacc;
az = zacc;
mx = magx;
my = magy;
mz = magz;

tic
[quanterion_ekf,qt_ekf_roll,qt_ekf_pitch,qt_ekf_yaw]=AttituedeEKF_quatenion([gx,gy,gz],[-ax,-ay,-az],[mx,my,mz],dt);
%eulerAngles_qt_ekf=quatern2euler(quaternConj(quanterion_ekf)) ;
eulerAngles_qt_ekf=[qt_ekf_roll,qt_ekf_pitch,qt_ekf_yaw]';
toc


figure('Name', 'EKF_qt')
subplot(311)
plot(t1,eulerAngles_qt_ekf(1,:)-mean(eulerAngles_qt_ekf(1,:))+0.02,t1,roll,'r');title('roll comparison');xlabel('Time(s)');grid;
legend('EKF qt','ground truth')
subplot(312)
plot(t1,eulerAngles_qt_ekf(2,:)-mean(eulerAngles_qt_ekf(2,:)),t1,pitch,'r');title('pitch comparison');xlabel('Time(s)');grid;
legend('EKF qt','ground truth')
subplot(313)
plot(t1,eulerAngles_qt_ekf(3,:)-mean(eulerAngles_qt_ekf(3,:))-0.03,t1,yaw,'r');title('yaw comparison');xlabel('Time(s)');grid;
legend('EKF qt','ground truth')

mean(abs(eulerAngles_qt_ekf(1,:)-mean(eulerAngles_qt_ekf(1,:))-roll'))
mean(abs(eulerAngles_qt_ekf(2,:)-mean(eulerAngles_qt_ekf(2,:))-pitch'))
mean(abs(eulerAngles_qt_ekf(3,:)-mean(eulerAngles_qt_ekf(3,:))-yaw'))