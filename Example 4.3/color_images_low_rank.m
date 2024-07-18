clc;clear;
% M = 500; N = M; 
% R = 50;
%% Construct Simulated data
% AA=lowrankQ_randn_qsvd(M,N,R);
% % AA=lowrankQ_randn(M,N,R);
% A1 = AA(:,1:N); A3 =  AA(:,N+1:2*N);
% A2 = AA(:,2*N+1:3*N);A4 = AA(:,3*N+1:4*N);
% A = quaternion(A1,A2,A3,A4);
addpath('E:\Paper\Quaternion\Data\dataset-main\dataset27');
Im = im2double(imread('10.png'));
A = quaternion(Im(:,:,1),Im(:,:,2),Im(:,:,3)); % Original
[M,N] = size(A);
K =120; step = 10;
%% Calculate Q_SVD approximation
Error_all_qsvd = [];  time_all_qsvd = []; B_all_QSVD = cell(1,12);
 for k = 10:step:K
     tStart = tic;
   [Uq,Sigma,Vq]  = Q_SVD(A);
   Vk = Vq(:,1:k);
   B_QSVD = Uq(:,1:k) * Sigma(1:k, 1:k) *Vk';
   B_QSVD  = vector( B_QSVD );
      time_QSVD = toc(tStart);
        i = k/10;
          B_all_QSVD{1, i}=B_QSVD;
   RE_QSVD= norm(A(:) - B_QSVD(:)) / norm(A(:));
   Error_all_qsvd =[Error_all_qsvd; RE_QSVD];
 time_all_qsvd =[ time_all_qsvd;  time_QSVD];
 end

    Error_all_qcur =[];     time_all_qcur  = [];B_all_cur = cell(1,12);
   for k = 10:step:K
   tStart = tic;
 % cc = k ; rr = k;
 cc = ceil(k * log(k)) ; rr = ceil(k * log(k)) ;
    rows = randsample(M, cc);
    cols = randsample(N, rr);
        D_cols = A(:,cols);
        D_rows = A(rows,:);
        C = D_cols;  R =D_rows;
    % Update L = C * pinv_U * R
    [UC, SC, VC ] = Q_SVD(D_cols );
    pinvC = VC  * pinv(SC) * UC';
       [UR, SR, VR ] = Q_SVD(D_rows);
    pinvR = VR  * pinv(SR) * UR';
    U = pinvC * A * pinvR;
   B_QCUR = C * U * R;  B_QCUR  = vector( B_QCUR);
          time_QCUR = toc(tStart);
          i = k/10;
          B_all_cur{1, i} =B_QCUR;
   RE_QCUR= norm(A(:) - B_QCUR(:)) / norm(A(:));
       Error_all_qcur = [Error_all_qcur;RE_QCUR];
    time_all_qcur = [time_all_qcur;time_QCUR];
   end

   Error_all_cur_length =[];     time_all_cur_length  = [];B_all_curLength=cell(1,12);
   for k =  10:step:K
          tStart = tic;
         cc =ceil(k * log(k)) ; rr= ceil(k * log(k)) ; 
          % [C, U, R]  =CUR_leverage(A, k,c,r);
 % [C,R, U]  = qcur_deim(A, k);
    normX = norm(A(:))^2;
    %% Update L = C * pinv_U * R
  [C] = col_length(A,cc,normX);
  [R] = row_length(A,rr,normX);
    [UC, SC, VC ] = Q_SVD(C);
    pinvC = VC  * pinv(SC) * UC';
       [UR, SR, VR ] = Q_SVD(R);
    pinvR = VR  * pinv(SR) * UR';
    U = pinvC * A * pinvR;
   B_cur_length = C * U * R;B_cur_length = vector(B_cur_length);
            time_cur_length= toc(tStart);
               i = k/10;
            B_all_curLength {1, i}=B_cur_length;
   RE_cur_length= norm(A(:) - B_cur_length(:)) / norm(A(:));
       Error_all_cur_length = [Error_all_cur_length;RE_cur_length];
    time_all_cur_length = [time_all_cur_length;time_cur_length];
   end
    imname=['color_images_10','.mat'];
save(imname,'A', 'B_all_QSVD','B_all_cur','B_all_curLength',...
    'Error_all_qcur','time_all_qcur',...
    'Error_all_cur_length','time_all_cur_length',...
    'Error_all_qsvd','time_all_qsvd');
