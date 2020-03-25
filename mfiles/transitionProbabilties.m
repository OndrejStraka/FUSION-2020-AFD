P1 = [100 5   5   1;
     5   100 1   5;
     5   1   100 5;
     1   5   5   100];
P2 = [100 5   5   0;
     2   100 0   0;
     2   0   100 0;
     1   20  20  1];
P3 = [100 5   5   1;
      2   100 1   2;
      2   1   100 2;
      1   20  20  100];
%P = [100  20  20    1;
      %20 100   1   20;
      %20   1 100   20;
       %1  20   20 100];
%P = [1 5 2 6;
     %4 1 1 2;
     %2 0 3 1;
     %3 4 4 1]/10;
P = P3;
% normalizace
P = P./sum(P,1);

% hledani P1 a P2 pomoci Bayese
[V,D] = eig(P);
Pi = V(:,abs(diag(D)-1)<eps);
Pi = Pi/sum(Pi);

P1_Bayes = kron(eye(2),ones(1,2))*(P.*Pi')*kron(eye(2),ones(2,1));
P1_Bayes = P1_Bayes./sum(P1_Bayes,1);
P2_Bayes = kron(ones(1,2),eye(2))*(P.*Pi')*kron(ones(2,1),eye(2));
P2_Bayes = P2_Bayes./sum(P2_Bayes,1);
diff_Bayes = P - kron(P1_Bayes,P2_Bayes)
norm_diff_Bayes = norm(diff_Bayes,'fro')

return

% hledani P1 a P2 neomezenou minimalizaci |P - (P1 kron P2)|_F
P11 = P(1:2,1:2);
P12 = P(1:2,3:4);
P21 = P(3:4,1:2);
P22 = P(3:4,3:4);
permP = [P11(:)';P21(:)';P12(:)';P22(:)'];
%[U,S,V] = svd(permP);
%P1_min = reshape(-U(:,1)*sqrt(S(1,1)),2,2);
%P2_min = reshape(-V(:,1)*sqrt(S(1,1)),2,2);
% normalizace
%P1_min = P1_min./sum(P1_min,1);
%P2_min = P2_min./sum(P2_min,1);
%P1_min = P1_min/sum(P1_min(:,1));
%P2_min = P2_min*sum(P1_min(:,1));

%min_error = norm(P - kron(P1_min,P2_min),'fro')
%norm(permP - P1_min(:)*P2_min(:)','fro')
%
% hledani P1 a P2 omezenou minimalizaci |P - (P1 kron P2)|_F
objective = @(x) norm(permP - [x(1);1-x(1);x(2);1-x(2)]*[x(3) 1-x(3) x(4) 1-x(4)],'fro');
x0 = [0.5 0.5 0.5 0.5];
lb = [0 0 0 0];
ub = [1 1 1 1];
x = fmincon(objective,x0,[],[],[],[],lb,ub);
P1_mincon = [x(1) x(2);1-x(1) 1-x(2)]
P2_mincon = [x(3) x(4);1-x(3) 1-x(4)]
diff_mincon = P - kron(P1_mincon,P2_mincon)
norm_diff_mincon = norm(diff_mincon,'fro')
