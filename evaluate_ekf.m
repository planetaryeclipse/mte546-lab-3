function [x_ests,P_ests] = evaluate_ekf(x0,P0,z_actual,f,h,F,H,Q,R)
%EVALUATE_EKF evaluates the ekf with starting conditions x0 and P0 against
%measurement data z_actual using functions f, h, F, H (of form alpha(x))
%with constant process and measurement noise covariance of Q and R

n = size(x0,1);
m = size(z_actual,2);

x_ests = zeros(n,m);
P_ests = zeros(n,n,m);

x_ests(:,1) = x0;
P_ests(:,:,1) = P0;

for k=2:m
    x_prev = x_ests(:,k-1);
    P_prev = P_ests(:,:,k-1);

    xk_predict = f(x_prev);
    zk_predict = h(x_prev);

    Fk = F(x_prev);
    Pk_predict = Fk*P_prev*Fk' + Q;

    Hk = H(xk_predict);
    Kk = Pk_predict*Hk'*(Hk*Pk_predict*Hk'+R)^-1;

    x_ests(:,k) = xk_predict+Kk*(z_actual(:,k)-zk_predict);
    P_ests(:,:,k) = (eye(n)-Kk*Hk)*Pk_predict;
end
end