clc;clear;
M = 500; N = M; 
k = 10;
Error_all_qsvd = [];  time_all_qsvd = [];
error1 = []; error2 = []; error = [];Error_all_qcur = [];  time_all_qcur  = [];
   Error_all_cur_length =[];     time_all_cur_length  = [];
    Error_all_qcur =[];     time_all_qcur  = [];
%% Construct Simulated data
for M = [50 100 150 200 250 300 350 400 450 500]
    N = M; 
AA=lowrankQ_randn_qsvd(M,N,k);
% AA=lowrankQ_randn(M,N,R);
A1 = AA(:,1:N); A3 =  AA(:,N+1:2*N);
A2 = AA(:,2*N+1:3*N);A4 = AA(:,3*N+1:4*N);
A = quaternion(A1,A2,A3,A4);
E1 = randn(M,N);E2 =randn(M,N);
E3 =randn(M,N);E4 = randn(M,N);
EE = quaternion(E1,E2,E3,E4);
E = 0 * EE;
A_noise = A+E ;
% A_noise = A + E; step = 10;
%% Calculate Q_SVD approximation
     tStart = tic;
   [Uq,Sigma,Vq]  = Q_SVD(A_noise);
   Vk = Vq(:,1:k);
   B_QSVD = Uq(:,1:k) * Sigma(1:k, 1:k) *Vk';
      time_QSVD = toc(tStart);
   RE_QSVD= norm(A(:) - B_QSVD(:)) / norm(A(:));
      Error_all_qsvd =[Error_all_qsvd; RE_QSVD];
 time_all_qsvd =[ time_all_qsvd;  time_QSVD];

 %% Calculate Q_CUR approximation
   tStart = tic;
 % cc = k ; rr = k;
 cc = ceil(k * log(k)) ; rr = ceil(k * log(k)) ;
    rows = randsample(M, cc);
    cols = randsample(N, rr);
        D_cols = A_noise(:,cols);
        D_rows = A_noise(rows,:);
        C = D_cols;  R =D_rows;
    % Update L = C * pinv_U * R
    [UC, SC, VC ] = Q_SVD(D_cols );
    pinvC = VC  * pinv(SC) * UC';
       [UR, SR, VR ] = Q_SVD(D_rows);
    pinvR = VR  * pinv(SR) * UR';
    U = pinvC * A * pinvR;
   B_QCUR = C * U * R;
          time_QCUR = toc(tStart);
   RE_QCUR= norm(A(:) - B_QCUR(:)) / norm(A(:));
       Error_all_qcur = [Error_all_qcur;RE_QCUR];
    time_all_qcur = [time_all_qcur;time_QCUR];
%% CUR_length
          tStart = tic;
         c = ceil(k * log(k)) ; r=  ceil(k * log(k)) ; 
          % [C, U, R]  =CUR_leverage(A, k,c,r);
 % [C,R, U]  = qcur_deim(A, k);
    normX = norm(A_noise(:))^2;
    %% Update L = C * pinv_U * R
  [C,cols ] = col_length(A_noise,cc,normX);
  [R,rows] = row_length(A_noise,rr,normX);
    [UC, SC, VC ] = Q_SVD(C);
    pinvC = VC  * pinv(SC) * UC';
    [UR, SR, VR ] = Q_SVD(R);
    pinvR = VR  * pinv(SR) * UR';
    U = pinvC * A * pinvR;
   B_cur_length = C * U * R;
   time_cur_length = toc(tStart);
   RE_cur_length= norm(A(:) - B_cur_length(:)) / norm(A(:));
          Error_all_cur_length = [Error_all_cur_length;RE_cur_length];
    time_all_cur_length = [time_all_cur_length;time_cur_length];
    %    Error_all_cur_length = [Error_all_cur_length;RE_cur_length];
    % time_all_cur_length = [time_all_cur_length;time_cur_length];
end
plot(Error_all_qcur,'c-*','LineWidth',2); hold on
hold on
plot(Error_all_cur_length,'-m^','LineWidth',2); hold on
hold on
plot(Error_all_qsvd,'-bo','LineWidth',2); hold on

figure
plot(time_all_qcur,':c*','LineWidth',2); hold on
hold on
plot(time_all_cur_length,':m^','LineWidth',2); hold on
hold on
plot(time_all_qsvd,':bo','LineWidth',2); hold on
    imname=['Permutation_noise_length_sigma_no','.mat'];
save(imname,'A_noise','A', ...
    'Error_all_qcur','time_all_qcur',...
    'Error_all_cur_length','time_all_cur_length',...
    'Error_all_qsvd','time_all_qsvd');

