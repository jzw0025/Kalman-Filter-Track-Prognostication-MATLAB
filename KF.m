TLAB KF IMPLEMENTATION %%%%%%%%%%%%%%%%%

function [result] = KalmanFiltering(x,z,Fms,R)
% to Calculate the Kalman Filter
% inputs: x --- time data
%            z --- trackable data
%            Fms --- process noise tuning parameter
%            R  --- measurement noise tuning parameter
%            
% outputs:
%            result --- output estimated value
%
 
TS=1;                        % sampling time, the time between the measurements (change);
I=[1,0,0;0,1,0;0,0,1];  % the unit matrix definition;
f11=0;                       %the F11 in the system dynamic matrix;
f12=1;                       %the F12 in the system dynamic matrix;
f13=0;                       %the F13 in the system dynamic matrix;
f21=0;                       %the F21 in the system dynamic matrix;
f22=0;                       %the F22 in the system dynamic matrix;
f23=1;                       %the F23 in the system dynamic matrix;
syms T;
F=[f11,f12,f13;f21,f22,f23;0,0,0];                                          % System Dynamic Matrix;
Fm=[1+f11*TS,f12*TS,f13*TS;f21*TS,1+f22*TS,f23*TS;0,0,1]; % Fundamental Matrix;
Fmt=[1+f11*TS,f21*TS,0;f12*TS,1+f22*TS,0;f13*TS,f23*TS,1];% Transpose Fundamental Matrix(TS);
Fm1=[1+f11*T,f12*T,f13*T;f21*T,1+f22*T,f23*T;0,0,1];         % Transpose Fundamental Matrix (T);
Fmt2=[1+f11*T,f21*T,0;f12*T,1+f22*T,0;f13*T,f23*T,1];        % Transpose Fundamental Matrix (T);
Q=Fms*[0,0,0;0,0,0;0,0,1]; Y=Fm1*Q*Fmt2;                           % Calculate the Y for the intergration;
QK2=int(Y,T,0,TS);                                                                % intergrating from the time "0"to"sample time TS";
QK=double(QK2);                                                                  % double the intergrating result;
Xe=[0,0,0]';                                                                           % Initiate Matrix, the first measurement(change)
Fmt=[1,0,0;TS,1,0;0.5*TS*TS,TS,1];                                        %Transport Fundamental Matrix 
H=[1,0,0];
Ht=[1;0;0];
P=I*10^-6;                                                                            % Initiate Convariance Matrix
result = zeros(3, lengtth(x));
gain = zeros(3, lengtth(x));
res = zeros(1, length(x));
for t1=1:length(x);% the process for calculating the M,P,K;
    Xp=Fm*Xe;
    M=Fm*P*Fmt+QK;
    K=M*Ht/(H*M*Ht+R);
    Xe=Xp+K*(z(1,t1)-H*Xp);
    P=(I-K*H)*M;         %recursive KF calculation using Matrix
    result(:,t1) = Xe;
    gain(:,t1) = K;        % Kalman Gain
    res(:,t1) = z(1,t1)-H*Xp;  % residual
    V1=result(1,t1);     %  KF tracking value
    V2=result(2,t1);     %  KF tracking velocity
    V3=result(3,t1);     %  KF tracking acceleration
end
