%%%%%%%%%%%%%%%%
% FINAL REPORT %
%%%%%%%%%%%%%%%%
% In our dataset the order of units is according to the name of cities
% (178) and near the name there is the value for the index of  quality if
% life 

%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY ANALISYS %
%%%%%%%%%%%%%%%%%%%%%%%%

% load the matrix X (178 x 9)
QLCit24;
% sample 60 observations from your date of birth
% fix the random seed (23/12/2002)
rng(231202);
% permute the rows
rowp = randperm(178);
Xp = Q(rowp, :);
X = Xp(1:60,:);
u60 = ones(60,1);
% mean of variables
mean(X);
% range of variables
min(X);
max(X);

%%%%%%%%%%%%%%
% EXERCISE 1 %
%%%%%%%%%%%%%%

% 1. STANDARDIZE DATA
Jc = eye(60)-(1/60)*ones(60);
Sc = 1/60*X'*Jc*X;
D = diag(diag(Sc).^0.5);
Z = Jc * X * D^-1;
% variance and covarianze matrix of Z or correlation matrix of X
Sz = 1/60*Z'*Z;
heatmap(Sz);
colormap('summer');

% 2. COMPUTE THE EUCLIDIAN DISTANCE BETWEEN UNITS
D = squareform(pdist(Z, 'euclidean'));
% heatmap of the matrix of distances 
heatmap(D);
colormap('summer');
% nearest and farthest cities
% nearest cities
D_no_diag = D + diag(Inf(size(D, 1), 1)); 
[min_dist, idx_min] = min(D_no_diag(:));
[row_min, col_min] = ind2sub(size(D), idx_min); 
fprintf('The nearest cities are %d and %d with a distance of %.2f\n', row_min, col_min, min_dist);
% farthest cities
[max_dist, idx_max] = max(D(:)); % Trovare il valore massimo e l'indice
[row_max, col_max] = ind2sub(size(D), idx_max); % Convertire l'indice in coordinate di riga e colonna
fprintf('The farthest cities are %d and %d with a distance of %.2f\n', row_max, col_max, max_dist);

% 3. COMPUTE THE WSPP
% 2 clusters 
[U2,b2, a2, f2,iter2]=WSPP(D, 2, 50);
[pf2,Dw2,Db2] = psF(Z,U2);
% 4 clusters
[U4,b4, a4, f4, iter4]=WSPP(D, 4, 50);
[pf4,Dw4,Db4] = psF(Z,U4);
% 3 clusters: chosen
[U3,b3, a3, f3,iter3]=WSPP(D, 3, 50);
[pf3,Dw3,Db3] = psF(Z,U3);
Pe = b3*(u60*u60' - U3*U3') + a3*(U3*U3' - eye(60));
% heatmap of membership matrix U3
heatmap(U3)
colormap('summer');
% how many cities in each cluster
unitCounts1 = sum(U3);
% which are the cities for each cluster
find(U3(:,1));
find(U3(:,2));
find(U3(:,3));
% find range of each cluster
indices1 = find(U3(:,1)== 1);
Y1 = X(indices1, :);
min(Y1(:,1))
max(Y1(:,1))
indices2 = find(U3(:,2)== 1);
Y2 = X(indices2, :);
min(Y2(:,1))
max(Y2(:,1))
indices3 = find(U3(:,3)== 1);
Y3 = X(indices3, :);
min(Y3(:,1))
max(Y3(:,1))
% silhouette plot
[~, clusterAssignments1] = max(U3, [], 2); 
unitCounts1 = histcounts(clusterAssignments1, 1:(3+1));
figure;
silhouette(D,clusterAssignments1);
title('Silhouett plot');
xlabel('Silhouette index');
ylabel('Cluster');
mean(silhouette(Z,clusterAssignments1));


% 4. COMPUTE THE WSP
% 2 clusters
[U2,Db2, Dw2, f2,iter2]=WSP(D, 2, 50);
[pf2,Dw2f,Db2f] = psF(Z,U2);
% 3 clusters
[U3,Db3, Dw3, f3,iter3]=WSP(D, 3, 50);
[pf3,Dw3f,Db3f] = psF(Z,U3);
% 4 clusters
[U4,Db4, Dw4, f4,iter4]=WSP(D, 4, 50);
[pf4,Dw4f,Db4f] = psF(Z,U4);
% heatmap of membership matrix U3
heatmap(U3)
colormap('summer');
% how many cities in each cluster
unitCounts1 = sum(U3);
% which are the cities for each cluster
find(U3(:,1));
find(U3(:,2));
find(U3(:,3));
% silhouette plot
[~, clusterAssignments1] = max(U3, [], 2); 
unitCounts1 = histcounts(clusterAssignments1, 1:(3+1));
figure;
silhouette(D,clusterAssignments1);
title('Silhouett plot');
xlabel('Silhouette index');
ylabel('Cluster');
mean(silhouette(Z,clusterAssignments1));


% 5. COMPUTE THE PD 
% 2 clusters 
[U2pd,Db2pd, Dw2pd, f2pd,iter2pd]=PD(D, 2, 50);
[pf2pd,Dw2pd,Db2pd] = psF(Z,U2pd);
% 3 clusters
[U3pd,Db3pd, Dw3pd, f3pd,iter3pd]=PD(D, 3, 50);
[pf3pd,Dw3pd,Db3pd] = psF(Z,U3pd);
% 4 clusters
[U4pd,Db4pd, Dw4pd, f4pd,iter4pd]=PD(D, 4, 50);
[pf4pd,Dw4pd,Db4pd] = psF(Z,U4pd);
% heatmap of membership matrix U3
heatmap(U3pd)
colormap('summer');
% how many cities in each cluster
unitCounts1 = sum(U3pd);
% which are the cities for each cluster
find(U3pd(:,1));
find(U3pd(:,2));
find(U3pd(:,3));
% silhouette plot
[~, clusterAssignments1] = max(U3pd, [], 2); 
unitCounts1 = histcounts(clusterAssignments1, 1:(3+1));
figure;
silhouette(D,clusterAssignments1);
title('Silhouett plot');
xlabel('Silhouette index');
ylabel('Cluster');
mean(silhouette(Z,clusterAssignments1));
% verify ultrametricity --> if it is 0 the matrix is ultrametric 
Q3 = U3pd*Db3pd*U3pd'+U3pd*Dw3pd*U3pd'-diag(diag(U3pd*Dw3pd*U3pd'));
Verultrametrica(Q3);

%%%%%%%%%%%%%%
% EXERCISE 2 %
%%%%%%%%%%%%%%

% 2. COMPUTE PCA TO DETERMINE THE NUMBER OF PRINCIPAL COMPONENTS 
% cambia A con T
[A, L] = eigs(Sz, 9);
diag(L);
% in L we have 3 components according to the Kaiser rule since we have 3
% large eigenvalues, greater than 1
% so we suppose to have 3 components
A3 = A(:,1:3);
A3r = rotatefactors(A3);
Y = Z * A3r;
% scatter plot
figure;

% Componente 1 vs Componente 2
subplot(1,3,1);
scatter(Y(:,1), Y(:,2), 50, 'filled');
xlabel('Component 1');
ylabel('Component 2');
title('Component 1 vs Component 2');
grid on;

% Component 1 vs Component 3
subplot(1,3,2);
scatter(Y(:,1), Y(:,3), 50, 'filled');
xlabel('Component 1');
ylabel('Component 3');
title('Component 1 vs Component 3');
grid on;

% Component 2 vs Component 3
subplot(1,3,3);
scatter(Y(:,2), Y(:,3), 50, 'filled');
xlabel('Component 2');
ylabel('Component 3');
title('Component 2 vs Component 3');
grid on;


% 3. COMPUTE K-MEANS ON THE NUMBER OF COMPONENT IDENTIFIED IN THE STEP 1
[loopOtt,UOtt,fOtt,iterOtt] = kmeansVICHI(Y,2,50);
[pf,Dw,Db] = psF(Z,UOtt);
% silhouette plot 
[~, clusterAssignments1] = max(UOtt, [], 2); 
unitCounts1 = histcounts(clusterAssignments1, 1:(2+1));
figure;
silhouette(Z,clusterAssignments1);
title('Silhouett plot');
xlabel('Silhouette index');
ylabel('Cluster');
mean(silhouette(Z,clusterAssignments1));

% 4. IDENTIFY THE BEST K WITH pF
K = 10;
pfKm = zeros(K,3);
for k = 2:K
    [~,UOtt,~,~]=kmeansVICHI(Y,k,50);
    [pF,Dw,Db]=[psF,Dw,Db];
    pfKm(k,:)=[psF,Dw,Db];
end
% k = 2
[loopOtt2,UOtt2,fOtt2,iterOtt2]=kmeansVICHI(Y,2,100);
[pfOtt2,Dw2,Db2]=psF(Y,UOtt2);
% k = 3
[loopOtt3,UOtt3,fOtt3,iterOtt3]=kmeansVICHI(Y,3,100);
[pfOtt3,Dw3,Db3]=psF(Y,UOtt3);
% k = 4
[loopOtt4,UOtt4,fOtt4,iterOtt4]=kmeansVICHI(Y,4,100);
[pfOtt4,Dw4,Db4]=psF(Y,UOtt4);
% k = 5
[loopOtt5,UOtt5,fOtt5,iterOtt5]=kmeansVICHI(Y,5,100);
[pfOtt5,Dw5,Db5]=psF(Y,UOtt5);
% k = 6
[loopOtt6,UOtt6,fOtt6,iterOtt6]=kmeansVICHI(Y,6,100);
[pfOtt6,Dw6,Db6]=psF(Y,UOtt6);
% k = 7
[loopOtt7,UOtt7,fOtt7,iterOtt7]=kmeansVICHI(Y,7,100);
[pfOtt7,Dw7,Db7]=psF(Y,UOtt7);
% k = 8
[loopOtt8,UOtt8,fOtt8,iterOtt8]=kmeansVICHI(Y,8,100);
[pfOtt8,Dw8,Db8]=psF(Y,UOtt8);
% k = 9
[loopOtt9,UOtt9,fOtt9,iterOtt9]=kmeansVICHI(Y,9,100);
[pfOtt9,Dw9,Db9]=psF(Y,UOtt9);
% k = 10
[loopOtt10,UOtt10,fOtt10,iterOtt10]=kmeansVICHI(Y,10,100);
[pfOtt10,Dw10,Db10]=psF(Y,UOtt10);
% store it in a vector
c = [pfOtt2, pfOtt3, pfOtt4, pfOtt5, pfOtt6, pfOtt7, pfOtt8, pfOtt9, pfOtt10];
figure;
plot(2:10, c, '-o');
title('Screeplot');
xlabel('Number of Clusters');
ylabel('pseudoF value');
% 2 clusters
heatmap(UOtt2);
colormap('summer');
sum(UOtt2);
find(UOtt2(:,1));
find(UOtt2(:,2));
% 3 clusters
heatmap(UOtt3);
colormap('summer');
sum(UOtt3);
find(UOtt3(:,1));
find(UOtt3(:,2));
find(UOtt3(:,3));

% 5. COMPUTE REDUCED K-MEAN;
% k = 2 BETTER
[Urkm,Arkm, Yrkm,frkm,inrkm]=REDKM(Z, 2, 3, 50);
pFrkm = psF(Yrkm,Urkm);
% confusion matrix
U2'*Urkm;
% k = 3
[Urkm,Arkm, Yrkm,frkm,inrkm]=REDKM(Z, 3, 3, 50);
pFrkm = psF(Yrkm,Urkm);
% heatmap for U
heatmap(Urkm);
colormap('summer');

% 6. COMPUTE FACTORIAL K-MEAN;
[Ufkm,Afkm, Yfkm,ffkm,infkm]=FKM(Z, 2, 3, 20);
pFfkm = psF(Yfkm,Ufkm);
% pseudoF better with k = 2
% confusion matrix
U2'*Ufkm;
heatmap(Ufkm);
colormap('summer');
% find cities in clusters
sum(Ufkm);
find(Ufkm(:,1))
find(Ufkm(:,2))
% rotation of factors
rotatefactors(Afkm);

% 7. COMPUTE CDPCA
[Vcdpca,Ucdpca,Acdpca, Ycdpca,fcdpca,incdpca]=CDPCA(Z, 2, 3, 50);
heatmap(Ucdpca);
colormap('summer');
heatmap(Vcdpca);
colormap('summer');

% 8. COMPUTE DKM
% I choose k = 2 from the application of KM
% I choose q = 2 from the applicationo of KM on columns of the Z
% choice of q
[loopd,UOttd,fOttd,iterOttd]=kmeansVICHI(Z',2,50);
[pfOttd,Dwd,Dbd]=psF(Z',UOttd);
[Vdkm,Udkm,Ymdkm,fdkm,indkm] = DKM(Z, 2, 2, 50);
heatmap(Udkm);
colormap('summer');
heatmap(Vdkm);
colormap('summer');

% 9. COMPARE THE RESULTS BY USING THE CONFUSION MATRIX (CONTINGENCY TABLE)
% BETWEEN 4 AND 5, 6 AND 7
% confusion matrix 4 and 5
UOtt*Urkm;
% confusion matrix 6 and 7
Ufkm'*Ucdpca;

% 10. COMPUTE THE EXPLAINED VARIANCE FROM 4 TO 8
% 4
Db/(sum(var(Y)));
