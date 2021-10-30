clear all
close all
clc

load Eulerdata

roll = roll * (pi/180);
pitch = pitch * (pi/180);
yaw = yaw * (pi/180);


dt = 20e-3; %50hz记录

[m,n]=size(xacc);
t1 = ((0:(m-1))*dt)';
t_fin = (m-1)*dt;

% mahony complementary filtering
gx = rollrate;
gy = pitchrate;
gz = yawrate;
ax = xacc;
ay = yacc;
az = -zacc;
mx = magx;
my = magy;
mz = -magz;


tic
halfT = 0.5 * dt;
Kp = 0.5;
Ki = 5;

exInt = zeros(m,1);							 
eyInt = zeros(m,1);
ezInt = zeros(m,1); 

q0 = zeros(m,1);
q1 = zeros(m,1);
q2 = zeros(m,1);
q3 = zeros(m,1);
q0q0 = zeros(m,1);
q0q1 = zeros(m,1);
q0q2 = zeros(m,1);
q0q3 = zeros(m,1);
q1q1 = zeros(m,1);
q1q2 = zeros(m,1);
q1q3 = zeros(m,1);
q2q2 = zeros(m,1);
q2q3 = zeros(m,1);
q3q3 = zeros(m,1);

q0(1) = 1;

rollmahony = zeros(m,1);
pitchmahony = zeros(m,1);
yawmahony = zeros(m,1);

az = -az;

for i=2:m
  
    q0q0(i)=q0(i-1)*q0(i-1);
    q0q1(i)=q0(i-1)*q1(i-1);
    q0q2(i)=q0(i-1)*q2(i-1);
    q0q3(i)=q0(i-1)*q3(i-1);
    q1q1(i)=q1(i-1)*q1(i-1);
    q1q2(i)=q1(i-1)*q2(i-1);
    q1q3(i)=q1(i-1)*q3(i-1);
    q2q2(i)=q2(i-1)*q2(i-1);
    q2q3(i)=q2(i-1)*q3(i-1);
    q3q3(i)=q3(i-1)*q3(i-1);

    norm1 = sqrt(mx(i)*mx(i)+my(i)*my(i)+mz(i)*mz(i));
    mx(i)=mx(i)/norm1;
    my(i)=my(i)/norm1;
    mz(i)=mz(i)/norm1;

    hx = 2.0 * (mx(i) * (0.5 - q2q2(i) - q3q3(i)) + my(i) * (q1q2(i) - q0q3(i)) + mz(i) * (q1q3(i) + q0q2(i)));  
    hy = 2.0 * (mx(i) * (q1q2(i) + q0q3(i)) + my(i) * (0.5 - q1q1(i) - q3q3(i)) + mz(i) * (q2q3(i) - q0q1(i)));  
    bx = sqrt(hx * hx + hy * hy);  
    ez_ef = -hy * bx;
    
    mag_ex = 2.0 * (q1q3(i) - q0q2(i)) * ez_ef;
    mag_ey = 2.0 * (q2q3(i) + q0q1(i)) * ez_ef;
    mag_ez = (1.0 - 2.0 * q1q1(i) - 2.0 * q2q2(i)) * ez_ef;
    
    norm1=sqrt(ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i));%把加计的三维向量转成单位向量
    ax(i)=ax(i)/norm1;
    ay(i)=ay(i)/norm1;
    az(i)=az(i)/norm1;
    
    halfex = (ay(i)*(1.0 - 2.0*q1q1(i) - 2.0*q2q2(i)) - az(i)*(2.0*(q2q3(i)+q0q1(i))));                          					
    halfey = (az(i)*(2.0*(q1q3(i)-q0q2(i))) - ax(i)*(1.0 - 2.0*q1q1(i) - 2.0*q2q2(i)));
    halfez = (ax(i)*(2.0*(q2q3(i)+q0q1(i))) - ay(i)*(2.0*(q1q3(i)-q0q2(i))));
    
    ex = halfex;% + mag_ex;
    ey = halfey;% + mag_ey;
    ez = halfez;% + mag_ez;
    
    exInt(i) = exInt(i) + ex* Ki* dt;								 
    eyInt(i) = eyInt(i) + ey* Ki* dt;
    ezInt(i) = ezInt(i) + ez* Ki* dt; 

    gx(i) = (gx(i) + Kp*ex + exInt(i))*halfT;					   							
    gy(i) = (gy(i) + Kp*ey + eyInt(i))*halfT;
    gz(i) = (gz(i) + Kp*ez + ezInt(i))*halfT;

    q0(i) = q0(i-1) + (-q1(i-1)*gx(i)  - q2(i-1)*gy(i) - q3(i-1)*gz(i));%四元数微分方程
    q1(i) = q1(i-1) + ( q0(i-1)*gx(i) + q2(i-1)*gz(i) - q3(i-1)*gy(i));
    q2(i) = q2(i-1) + ( q0(i-1)*gy(i)  - q1(i-1)*gz(i) + q3(i-1)*gx(i));
    q3(i) = q3(i-1) + ( q0(i-1)*gz(i) + q1(i-1)*gy(i)  - q2(i-1)*gx(i));

    norm1 = sqrt(q0(i)*q0(i) + q1(i)*q1(i) + q2(i)*q2(i) + q3(i)*q3(i));%四元数规范化
    q0(i) = q0(i) / norm1;
    q1(i) = q1(i) / norm1;
    q2(i) = q2(i) / norm1;
    q3(i) = q3(i) / norm1;

    rollmahony(i)   = atan2(2 * q2(i) * q3(i) + 2 * q0(i) * q1(i), -2 * q1(i) * q1(i) - 2 *q2(i)* q2(i) + 1); %roll
    pitchmahony(i)  = (asin(2 * q0(i) * q2(i) - 2 * q1(i)* q3(i))); % pitch
    yawmahony(i)    = atan2(2 * q1(i) * q2(i) + 2 * q0(i) * q3(i) , -2 *q2(i) * q2(i) - 2 * q3(i) * q3(i) +1);
end

toc
figure('Name', 'mahony')
subplot(311)
plot(t1,rollmahony,t1,roll,'r');title('roll comparison');xlabel('Time(s)');grid;
legend('mahony','ground truth')
subplot(312)
plot(t1,pitchmahony,t1,pitch,'r');title('pitch comparison');xlabel('Time(s)');grid;
legend('mahony','ground truth')
subplot(313)
plot(t1,yawmahony,t1,yaw,'r');title('yaw comparison');xlabel('Time(s)');grid;
legend('mahony','ground truth')


mean(abs(rollmahony-roll))
mean(abs(pitchmahony-pitch))
mean(abs(yawmahony-yaw))