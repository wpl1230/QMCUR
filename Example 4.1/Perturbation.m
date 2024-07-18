clc; clear;
clc;clear;
M = 500; N = M; 
R = 50;
AA=lowrankQ_randn_qsvd(M,N,R);
% AA=lowrankQ_randn(M,N,R);
A1 = AA(:,1:N); A3 =  AA(:,N+1:2*N);
A2 = AA(:,2*N+1:3*N);A4 = AA(:,3*N+1:4*N);
A = quaternion(A1,A2,A3,A4);
%% Construct Simulated data
% A = randq(M,R) * randq(R,N);
% % E1 = rand(M,N);E2 = rand(M,N);E3 = rand(M,N);E4 = rand(M,N);
% M1 = randn(M,R);M2 = randn(M,R);M3 = randn(M,R);M4 = randn(M,R);
% N1 = randn(R,N);N2 = randn(R,N);N3 = randn(R,N);N4 = randn(R,N);
% MM =quaternion(M1,M2,M3,M4);NN =quaternion(N1,N2,N3,N4);
% A = MM * NN;
% E1 =0.1 * randn(M,N);E2 = 0.1 *randn(M,N);E3 =0.1 * randn(M,N);E4 =0.1 * randn(M,N);
% E = quaternion(E1,E2,E3,E4);
% E = randq(M,N); 
normE = [];
E1 = randn(M,N);E2 =randn(M,N);
E3 =randn(M,N);E4 = randn(M,N);
EE = quaternion(E1,E2,E3,E4);
Error_all_qsvd = [];  time_all_qsvd = [];
error1 = []; error2 = []; error = [];Error_all_qcur = [];  time_all_qcur  = [];
for degree = [1e-6 0.000005 1e-5 0.00005 0.0001 0.0005 0.001 0.005  0.01 0.05 0.1]
E = degree * EE;
Ahat = A +   E;
% A = randq(M,N);
K =50; step = 10;
normE = [normE, norm(E)];
for k =50%10:step :K
tStart = tic;
 % cc = k ; rr = k;
 cc = ceil(k * log(k)) ; rr = ceil(k * log(k)) ;
    rows = randsample(M, cc);
    cols = randsample(N, rr);
    C= Ahat(:,cols);
    R = Ahat(rows,:);
    % Update L = C * pinv_U * R
    [UC, SC, VC ] = Q_SVD(C );
    pinvC = VC  * pinv(SC) * UC';
       [UR, SR, VR ] = Q_SVD(R);
    pinvR = VR  * pinv(SR) * UR';
    U = pinvC * A * pinvR;
   B_QCUR = C * U * R;
          time_QCUR = toc(tStart);
              time_all_qcur = [time_all_qcur;time_QCUR];
   RE_QCUR= norm(A - B_QCUR) / norm(A);
    Error_all_qcur =[Error_all_qcur;   RE_QCUR];
[Wa,Sa,Va] = Q_SVD(A);
Vak = Va(:,1:k);
Wak = Wa(:,1:k);
f1 =  norm(A- B_QCUR) ;
error1 = [error1,f1];
f2 = norm(E) * (norm(Wak(rows,:)) + norm(Vak(cols,:)) +2);
error2 = [error2,f2];
dis = f2 - f1;
error  = [error,dis];
end
% plot(error1,'b'); hold on
% plot(error2,'r')
% 
% Ahat = A +  degree * E; k = 50;
 % for k = 50%:step:K
     tStart = tic;
   [Uq,Sigma,Vq]  = Q_SVD(Ahat);
   Vk = Vq(:,1:k);
   B_QSVD = Uq(:,1:k) * Sigma(1:k, 1:k) *Vk';
      time_QSVD = toc(tStart);
   RE_QSVD= norm(A- B_QSVD) / norm(A);
      Error_all_qsvd =[Error_all_qsvd; RE_QSVD];
 time_all_qsvd =[ time_all_qsvd;  time_QSVD];
 % end
end
 imname=['Permutation_noise','.mat'];
save(imname,'Ahat','A','normE','error','error1','error2', ...
    'Error_all_qcur','time_all_qcur','Error_all_qsvd','time_all_qsvd');