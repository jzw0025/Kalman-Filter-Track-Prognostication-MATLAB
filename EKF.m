function result=EKF(x,z,Fms,R)
% to Calculate the Kalman Filter
% inputs: x --- time data
%            z --- trackable data
%            Fms --- process noise tuning parameter
%            R  --- measurement noise tuning parameter
%            
% outputs:
%            result --- output estimated value
%

TS=1;
I=[1,0,0;0,1,0;0,0,1]; % unit matrix
Xp=ones(3,10);Xe=ones(3,10);%initiate the two matrix
Xe(1,1)=1;
H=[1,0,0];
Ht=[1;0;0];
P=0.1*I;% the initiate value for the covariance matrix
result = zeros(1,x);
for t=1:length(x);
    %%%  CALCULATING JACCOBIAN  %%%
    f11=Xe(3,t);% b
    f12=1;
    f13=Xe(1,t);% x
    f21=Xe(3,t)*Xe(3,t);% b^2
    f22=Xe(3,t); % b
    f23=2*Xe(3,t)*Xe(1,t);% 2*b*x
    
    syms T;
    F=[f11,f12,f13;f21,f22,f23;0,0,0]; % System Dynamic Matrix;
    Fm=[1+f11*TS,f12*TS,f13*TS;f21*TS,1+f22*TS,f23*TS;0,0,1]; % Fundamental Matrix;
    Fmt=[1+f11*TS,f21*TS,0;f12*TS,1+f22*TS,0;f13*TS,f23*TS,1];% Transpose Fundamental Matrix(TS);
    Fm1=[1+f11*T,f12*T,f13*T;f21*T,1+f22*T,f23*T;0,0,1];% Transpose Fundamental Matrix (T);
    Fmt2=[1+f11*T,f21*T,0;f12*T,1+f22*T,0;f13*T,f23*T,1];% Transpose Fundamental Matrix (T);
    
    Q=Fms*[0,0,0;0,0,0;0,0,1]; 
    Y=Fm1*Q*Fmt2;
    QK2=int(Y,T,0,TS);
    QK3=double(QK2);
    
    Xp(3,t)=Xe(3,t)*Xe(2,t);      %get DDx at the time t;
    Xp(1,t+1)=Xe(1,t)+TS*Xp(2,t); %get x at the time t+1;
    Xp(2,t+1)=Xe(2,t)+TS*Xp(3,t); %get Dx at the time t+1;
    
    M=Fm*P*Fmt+QK3;  % calculate for the M matrix
    K=M*Ht/(H*M*Ht+R);% calculate for the Kalman Gain
    
    Xe(1,t+1)=Xp(1,t+1)+K(1,1)*(z(line,t)-Xp(1,t+1));% the next time x estimation
    Xe(2,t+1)=Xp(2,t+1)+K(2,1)*(z(line,t)-Xp(1,t+1));% the next time Dx estimation
    Xe(3,t+1)=Xe(3,t)+K(3,1)*(z(line,t)-Xp(1,t+1));  % the next time coefficient 'b' estimation
    
    P=(I-K*H)*M; % calculate for the Convariance Matrix after the Update;
    
    result(:,t) = Xe(:,t); 

end

