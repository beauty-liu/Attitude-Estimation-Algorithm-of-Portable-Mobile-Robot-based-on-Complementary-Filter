function [ quatention,Roll,Pitch,Yaw ] = AttituedeEKF_quatenion(gyro,acc,mag,dt)
%EKF 4state
[m n]=size(gyro);
quatention=zeros(m,4);

p=gyro(:,1);
q=gyro(:,2);
r=gyro(:,3);

ax=acc(:,1);
ay=acc(:,2);
az=acc(:,3);

mx=mag(:,1);
my=mag(:,2);
mz=mag(:,3);

%ACC�Ĳ�������
var_ax=5;
var_ay=5;
var_az=5;
%mag�Ĳ�������
var_mx=5;
var_my=5;
var_mz=5;
%״̬�����Ĺ�������
P1=[2,2,2,2];
P=diag(P1);

x_e=[1;0;0;0];
Pitch = zeros(m,1);
Roll = zeros(m,1);
Yaw = zeros(m,1);
 %״̬ת�ƾ���
F = eye(4,4);
%�������������
K = zeros(4,3);
%�۲����
H = zeros(3,4);
Hdt=0.5*dt;
Q1 = 0.0001;%Q�����е���Ҳ�Ǿ���ֵ�õ�����Ҫ�ڿ���q�еĲ���
Q2 = 0.0009;

R=zeros(6,6);
R(1,1)=var_ax;
R(2,2)=var_ay;
R(3,3)=var_az;
R(4,4)=var_mx;
R(5,5)=var_my;
R(6,6)=var_mz;

% R=zeros(3,3);
% R(1,1)=var_ax;
% R(2,2)=var_ay;
% R(3,3)=var_az;

for i=1:m
    norm = 1.0/sqrt(x_e(1,1)*x_e(1,1)+x_e(2,1)*x_e(2,1)+x_e(3,1)*x_e(3,1)+x_e(4,1)*x_e(4,1));
    for n=1:4
        x_e(n,1)=x_e(n,1)*norm;
    end
    
    %״̬ת�ƾ���
    F = [              1,  -Hdt*(p(i)),  -Hdt*(q(i)),  -Hdt*(r(i));
        Hdt*(p(i)),                1,  Hdt*(r(i)), -Hdt*(q(i));
        Hdt*(q(i)), -Hdt*(r(i)),                1,  Hdt*(p(i));
        Hdt*(r(i)),  Hdt*(q(i)), -Hdt*(p(i)),                1;];
 
    dTemp(1)=x_e(1,1)-Hdt*(p(i))*x_e(2,1)-Hdt*(q(i))*x_e(3,1)-Hdt*(r(i))*x_e(4,1);
    dTemp(2)=x_e(2,1)+Hdt*(p(i))*x_e(1,1)-Hdt*(q(i))*x_e(4,1)+Hdt*(r(i))*x_e(3,1);
    dTemp(3)=x_e(3,1)+Hdt*(p(i))*x_e(4,1)+Hdt*(q(i))*x_e(1,1)-Hdt*(r(i))*x_e(2,1);
    dTemp(4)=x_e(4,1)-Hdt*(p(i))*x_e(3,1)+Hdt*(q(i))*x_e(2,1)+Hdt*(r(i))*x_e(1,1);

    for j=1:4
        x_e(j,1)=dTemp(j);
    end
    %Э����һ��Ԥ�����
    P=F*P*F';
    %����������
    P(1,1)=P(1,1)+Q1;
    P(2,2)=P(2,2)+Q1;
    P(3,3)=P(3,3)+Q1;
    P(4,4)=P(4,4)+Q1;

    
    y(1,1)= -2.0*(x_e(2,1)*x_e(4,1)-x_e(1,1)*x_e(3,1));
    y(2,1)= -2.0*(x_e(1,1)*x_e(2,1)+x_e(3,1)*x_e(4,1));
    y(3,1)= -(x_e(1,1)*x_e(1,1)-x_e(2,1)*x_e(2,1)-x_e(3,1)*x_e(3,1)+x_e(4,1)*x_e(4,1));
    %������Acc����
    norm=sqrt(ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i));%�ѼӼƵ���ά����ת�ɵ�λ����
    ax(i)=ax(i)/norm;
    ay(i)=ay(i)/norm;
    az(i)=az(i)/norm;
    y_s(1,1) = ax(i);
    y_s(2,1) = ay(i);
    y_s(3,1) = az(i);
    
    norm = sqrt(mx(i)*mx(i)+my(i)*my(i)+mz(i)*mz(i));
    mx(i)=mx(i)/norm;
    my(i)=my(i)/norm;
    mz(i)=mz(i)/norm;
    y_s(4,1) = mx(i);
    y_s(5,1) = my(i);
    y_s(6,1) = mz(i);

    hx = 2.0 * (mx(i) * (0.5 - x_e(3,1)*x_e(3,1) - x_e(4,1)*x_e(4,1)) + my(i) * (x_e(2,1)*x_e(3,1) - x_e(1,1)*x_e(4,1)) + mz(i) * (x_e(2,1)*x_e(4,1) + x_e(1,1)*x_e(3,1)));  
    hy = 2.0 * (mx(i) * (x_e(2,1)*x_e(3,1) + x_e(1,1)*x_e(4,1)) + my(i) * (0.5 - x_e(2,1)*x_e(2,1) - x_e(4,1)*x_e(4,1)) + mz(i) * (x_e(3,1)*x_e(4,1) - x_e(1,1)*x_e(2,1)));  
    bx = sqrt(hx * hx + hy * hy);  
    bz = 2.0 * (mx(i) * (x_e(2,1)*x_e(4,1) - x_e(1,1)*x_e(3,1)) + my(i) * (x_e(3,1)*x_e(4,1) + x_e(1,1)*x_e(2,1)) + mz(i) * (0.5 - x_e(2,1)*x_e(2,1) - x_e(3,1)*x_e(3,1)));
     
    y(4,1)=bx*(x_e(1,1)*x_e(1,1) + x_e(2,1)*x_e(2,1) - x_e(3,1)*x_e(3,1) - x_e(4,1)*x_e(4,1)) + 2*bz*(x_e(2,1)*x_e(4,1) - x_e(1,1)*x_e(3,1));
    y(5,1)=2*bx*(x_e(2,1)*x_e(3,1) - x_e(1,1)*x_e(4,1)) + 2*bz*(x_e(3,1)*x_e(4,1) + x_e(1,1)*x_e(2,1));
    y(6,1)=2*bx*(x_e(2,1)*x_e(4,1) + x_e(1,1)*x_e(3,1)) + bz*(x_e(1,1)*x_e(1,1) - x_e(2,1)*x_e(2,1) - x_e(3,1)*x_e(3,1) + x_e(4,1)*x_e(4,1));

        %�۲�������
%         H = [ 2.0*x_e(3,1), -2.0*x_e(4,1),  2.0*x_e(1,1), -2.0*x_e(2,1);
%             -2.0*x_e(2,1), -2.0*x_e(1,1), -2.0*x_e(4,1), -2.0*x_e(3,1);
%             -2.0*x_e(1,1),  2.0*x_e(2,1),  2.0*x_e(3,1), -2.0*x_e(4,1);]; 
        
        H = [ 2*x_e(3,1), -2*x_e(4,1),  2*x_e(1,1), -2*x_e(2,1);
            -2*x_e(2,1), -2*x_e(1,1), -2*x_e(4,1), -2*x_e(3,1);
            -2*x_e(1,1),  2*x_e(2,1),  2*x_e(3,1), -2*x_e(4,1);
            2*( x_e(1,1)*bx - x_e(3,1)*bz), 2*( x_e(2,1)*bx + x_e(4,1)*bz), 2*(-x_e(3,1)*bx - x_e(1,1)*bz), 2*(-x_e(4,1)*bx + x_e(2,1)*bz);
            2*(-x_e(4,1)*bx + x_e(2,1)*bz), 2*( x_e(3,1)*bx + x_e(1,1)*bz), 2*( x_e(2,1)*bx + x_e(4,1)*bz), 2*(-x_e(1,1)*bx + x_e(3,1)*bz);
            2*( x_e(3,1)*bx + x_e(1,1)*bz), 2*( x_e(4,1)*bx - x_e(2,1)*bz), 2*( x_e(1,1)*bx - x_e(3,1)*bz), 2*( x_e(2,1)*bx + x_e(4,1)*bz)]; 
        
        %���㿨��������
        K=P*H'*pinv(H*P*H'+R);
        %״̬����
        x_e=x_e+K*(y_s-y);
        %x_e=x_e;
        P=P-K*H*P;

    
        norm = 1.0/sqrt(x_e(1,1)*x_e(1,1)+x_e(2,1)*x_e(2,1)+x_e(3,1)*x_e(3,1)+x_e(4,1)*x_e(4,1));
        for n=1:4
            x_e(n,1)=x_e(n,1)*norm;
        end
        
        quatention(i,:)=x_e';
        %��Ԫ��תŷ����
        Pitch(i) = asin(2*(x_e(1,1)*x_e(3,1)- x_e(2,1)*x_e(4,1)));
        Roll(i) = atan2(2*(x_e(1,1)*x_e(2,1)+x_e(3,1)*x_e(4,1)),1-2*(x_e(2,1)*x_e(2,1)+x_e(3,1)*x_e(3,1)));
        Yaw(i) = atan2(2*(x_e(2,1)*x_e(3,1)+x_e(1,1)*x_e(4,1)),1-2*(x_e(3,1)*x_e(3,1)+x_e(4,1)*x_e(4,1)));

end   


end

