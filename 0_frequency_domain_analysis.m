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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ground truth
figure('Name', 'Attitude:ground truth')

subplot(311)
plot(t1,roll);title('roll:ground truth');xlabel('Time(s)');ylabel('\phi(rad）');grid;axis([0 100,-inf,inf])
subplot(312)
plot(t1,pitch);title('pitch:ground truth');xlabel('Time(s)');ylabel('\theta(rad）');grid;axis([0 100,-inf,inf])
subplot(313)
plot(t1,yaw);title('yaw:ground truth');xlabel('Time(s)');ylabel('\psi(rad）');grid;axis([0 100,-inf,inf])

%measured data
figure('Name', 'measured data')
subplot(231)
plot(t1,rollrate);title('Gyroscope: Roll rate');xlabel('Time(s)');ylabel('Angular rate(rad/s)');grid;axis([0 100,-inf,inf])
subplot(232)
plot(t1,pitchrate);title('Gyroscope: Pitch rate');xlabel('Time(s)');grid;axis([0 100,-inf,inf])
subplot(233)
plot(t1,yawrate);title('Gyroscope: Yaw rate');xlabel('Time(s)');grid;axis([0 100,-inf,inf])
subplot(234)
plot(t1,xacc);title('Accelerometer: xacc');xlabel('Time(s)');ylabel('Acceleration(m/s^2)');grid;axis([0 100,-inf,inf])
subplot(235)
plot(t1,yacc);title('Accelerometer: yacc');xlabel('Time(s)');grid;axis([0 100,-inf,inf])
subplot(236)
plot(t1,zacc);title('Accelerometer: zacc');xlabel('Time(s)');grid;axis([0 100,-inf,inf])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%get the attitude angles by integrating the angular velocity from gyroscopes
sampleperiod=dt;
Gyroscope=[rollrate pitchrate yawrate];
gyro=zeros(length(t1),3);
gyro_integral = zeros(length(t1),3);
for t = 2:length(t1)
    gyro_integral(t,:) = gyro_integral(t-1,:) + Gyroscope(t,:)*sampleperiod;
end
figure('Name', 'integrating the angular velocity from gyroscopes');


plot(t1, gyro_integral(:,1), 'r',t1, gyro_integral(:,2), 'g',t1, gyro_integral(:,3), 'b');grid;axis([0 100,-inf,inf])
legend('roll\phi', 'pitch\theta', 'yaw\psi');
xlabel('Time (s)');
ylabel('Angular(rad)');


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

figure('Name', 'orthogonally decomposing the acceleration');
subplot(311)
plot(t1,acc_angle(:,1)*pi/180,t1,roll,'r');title('roll');xlabel('Time(s)');ylabel('\phi(rad）');grid;axis([0 100,-inf,inf])
legend('estimation','ground truth')
subplot(312)
plot(t1,acc_angle(:,2)*pi/180,t1,pitch,'r');title('pitch');xlabel('Time(s)');ylabel('\theta(rad）');grid;axis([0 100,-inf,inf])
legend('estimation','ground truth')
subplot(313)
plot(t1,(acc_angle(:,3)-mean(acc_angle(:,3)))*pi/180,t1,yaw,'r');title('yaw');xlabel('Time(s)');ylabel('\psi(rad）');grid;axis([0 100,-inf,inf])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency domain analysis

M = m;
N = M*2-1;
X_roll =fft(gyro_integral(:,1),N);
X_pitch =fft(gyro_integral(:,2),N);
X_yaw =fft(gyro_integral(:,3),N);
f=[-M+1:M-1]/(N*dt);

Y_roll =fft(acc_angle(:,1),N);
Y_pitch =fft(acc_angle(:,2),N);
Y_yaw =fft(acc_angle(:,3),N);
f=[-M+1:M-1]/(N*dt);

figure('Name', 'gyro frequence')
subplot(221)
plot(f,abs(fftshift(X_roll)));title('roll angle by integrating signal from gyroscope');
xlabel('frequency/(Hz)');ylabel('amplitude');
subplot(223)
plot(f,abs(fftshift(X_pitch)));title('pitch angle by integrating signal from gyroscope');
xlabel('frequency/(Hz)');ylabel('amplitude');

subplot(222)
plot(f,abs(fftshift(Y_roll)));title('roll angle by orthogonally decomposing the acceleration');
xlabel('frequency/(Hz)');ylabel('amplitude');
subplot(224)
plot(f,abs(fftshift(Y_pitch)));title('pitch angle by orthogonally decomposing the acceleration');
xlabel('frequency/(Hz)');ylabel('amplitude');



