function exercise3()
  
  #using point from slide here
  #x1 = [44.7; -103.6; 47.4; -152.2; -153.3; -149.4];
  #y1 = [-142.4; -146.6; -150.1; 59.4; -96.9; 52.7];
  #z1 = [258.7; 154.4; 59.8; 245.2; 151.3; 46.9];
  #x2 = [18.5; 99.1; 13.8; 242.1; 151.1; 243.1];
  #y2 = [46.8; 146.5; 221.8; 52.5; 147.1; 224.5];
  
  #exercise 3.3
  x1 = [4; 12; 0;  0; 4; 12];
  y1 = [4; 12; 4; 12; 0;  0];
  z1 = [0;  0; 4; 12; 4; 12];
  
  x2 = [755; 1028; 557; 291; 665; 688];
  y2 = [455; 297; 455; 305; 623; 967];
  
  T1 = coordinates4x4(x1, y1, z1)
  T2 = coordinates3x3(x2, y2)
  
  A = designMatrix(T1, T2, x1, y1, z1, x2, y2)
  
  P_tilde = solveEquation(A)
  
  #reverse conditioning
  P = reverseConditioning(P_tilde, T1, T2)
  
  #exerceise 3.4
  [C, R, K] = rqDecomposition(P)
  
function T = coordinates3x3(in_x, in_y)
  #get mean of translation
  t_x = mean(in_x);
  t_y = mean(in_y);
  
  #translate points to the origin
  new_x = abs(in_x - t_x);
  new_y = abs(in_y - t_y);
  
  #get mean of the scale
  s_x = mean(new_x);
  s_y = mean(new_y);
  
  #matrix multiplied with points
  M1 = [1/s_x 0    0;
        0    1/s_y 0;
        0     0    1];
       
  M2 = [1 0 -t_x;
        0 1 -t_y;
        0 0   1];
        
  T = M1 * M2;
  
function T = coordinates4x4(in_x, in_y, in_z)
  #get mean of translation
  t_x = mean(in_x);
  t_y = mean(in_y);
  t_z = mean(in_z);
  
  #translate points to the origin
  new_x = abs(in_x - t_x);
  new_y = abs(in_y - t_y);
  new_z = abs(in_z - t_z);
  
  #get mean of the scale
  s_x = mean(new_x);
  s_y = mean(new_y);
  s_z = mean(new_z);
  
  #matrix multiplied with points
  M1 = [1/s_x 0    0     0;
         0   1/s_y 0     0;
         0    0   1/s_z  0;
         0    0    0     1];
       
  M2 = [1 0 0 -t_x;
        0 1 0 -t_y;
        0 0 1 -t_z;
        0 0 0   1];
        
  T = M1 * M2;
  
function A = designMatrix(T1, T2, in_x1, in_y1, in_z1, in_x2, in_y2)
  A = zeros(rows(in_x1)*2, 12);
  for var = 1:rows(in_x1)
    #get x tilde = T * x
    out1 = T1 * [in_x1(var,1); in_y1(var,1); in_z1(var,1); 1]
    out2 = T2 * [in_x2(var,1); in_y2(var,1); 1]
    
    #use computed points to get the designmatrix
    #put the results in one big matrix
    A([(2*var)-1, 2*var], :) = designMatrixPart(out1(1,1), out1(2,1), out1(3,1), out2(1,1), out2(2,1));
  end

#this function computes one part of the design matrix  
function A = designMatrixPart(in_x1, in_y1, in_z1, in_x2, in_y2)
  u = in_x2 * [in_x1, in_y1, in_z1, 1];
  u = normalizePoint(u);
  
  v = in_y2 * [in_x1, in_y1, in_z1, 1];
  v = normalizePoint(v);
  
  w = -1 * [in_x1, in_y1, in_z1, 1];
  w = normalizePoint(w);
  
  #matrix composition like in slides
  A = [w 0 0 0 0 u;
       0 0 0 0 w v];
  
function p = normalizePoint(p_in)
  p = p_in / abs(p_in(1,3));

function P_tilde = solveEquation(A)
  [U, S, V] = svd (A);
  S_min = min(diag(S));
  
  for var = 1:rows(S)
    if S(var,var) == S_min
      row = var;
    end
  end
  
  #reshape the column of V
  P_tilde = reshape(V(:, row), 3, 4);
  
function P = reverseConditioning(P_tilde, T1, T2)
  P = inv(T2) * P_tilde * T1;
  
function [C, R, K] = rqDecomposition(P)
  M = [ P(1, 1:3);
        P(2, 1:3);
        P(3, 1:3)];
        
  m3 = [M(3, :)];
  
  det = det(M);
  if(det > 0)
    delta = 1/norm(m3);
  else
    delta = -1/norm(m3);
  endif
  
  [R, K] = qr (inv(M));
  R = inv(R);
  K = inv(K);
  
  #svd to find PC = 0
  [S, V, D] = svd(M);
  V_min = min(diag(V));
  
  for var = 1:rows(V)
    if V(var,var) == V_min
      row = var;
    end
  end     
  
  #reshape the column of D
  C = reshape(D(:, row), 3, 1);