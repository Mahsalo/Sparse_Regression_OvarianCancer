%%Final code
%%%%%% Concordance index %%%%%%%%%
%%reading TCGA agilent training an Excel file in MATLAB %%%%
training_tcga=xlsread('TCGA_time.xlsx');
clear
%%%%%%switching NA values to "zero"
[m,n]=find(training_tcga==0);
%%13104 genes and 578 samples
%%eliminating the samples that have NA recurrence time (m'), every column
%%is a patient
training_exprs= xlsread('TCGA_training.xlsx');
copy=training_exprs;
[m1,n1]=size(training_exprs);

% load tarining_data_ready_workspace

load('training_data_ready_workspace.mat')
copy(:,m)=[];%%eliminate columns with NA/0 time recurrence value,dimension of copy is (feature * sample)
Aprime=copy;
[u1,u2]=size(Aprime);

y=training_tcga;%%time recurrence of a tumor
yy=y;
yy(m)=[];

% Sorting the survival values "yy" and matrix gene expression "A"
%%accordingly 
A=zeros(u1,u2);
[yy,oldindex]=sort(yy,'ascend');
for u3=1:length(yy)
A(:,u3)=Aprime(:,oldindex(u3));
end

% Find common genes for all datasets (training+testing) 

[num1,GSE30161_genes]=xlsread('GSE30161_genes.xlsx');
[num2,GSE17260_genes]=xlsread('GSE17260_genes.xlsx');
[num3,GSE9891_genes]=xlsread('GSE9891_genes.xlsx');
[num4,TCGA_genenames]=xlsread('TCGA_genenames.xlsx');
[m1,n1]=size(GSE30161_genes);
[m2,n2]=size(GSE17260_genes);
[m3,n3]=size(GSE9891_genes);
[p1,p2]=ismember(GSE30161_genes,GSE9891_genes);
[p3,p4]=find(p1==1);
common1=GSE30161_genes(p1,1);
[p5,p6]=ismember(common1,GSE17260_genes);
[p7,p8]=find(p5==1);
common2=common1(p7,1);
[p9,p10]=ismember(TCGA_genenames,common2);
[p11,p12]=find(p9==1);
[p13,p14]=find(p9==0);
common_genes=TCGA_genenames(p9,1); %%there are 12249 common genes  
A(p13,:)=[];
[j1,j2]=size(A);%%12249*522 j2=sample , j1=feature
A=A';%%sample*feature or 522*12249
lambda=10;
mu=0.2;
data=A;
survival=log(yy);
[absolute_error nentries1 selected1  z1 weight_all predicted_survival randomization]= MENrfe2(data,survival)
[absolute_error nentries1 selected1  z1  bias weight_all predicted_survival randomization]= MENrfe(data,survival,lambda,mu);

%% Optimization

[normalmatrix] = normalize(A,j2);
lambda=0.1;
mu=0.2;
[w1 b] = optimizer(normalmatrix,yy,j1,lambda,mu)%%%lagrange
w1=optimizer2(normalmatrix,yy,12249,mu)%%%CLOT
w1=optimizer3(normalmatrix,yy,12249)%%Lasso

% Sorting w1 to see nonzero elements on the graph
%instead of optimizing again we can load a saved workspace
load('lambdayeksadom.mat')
zt=w1;
w1sorted=sort(abs(w1),'ascend');
figure;
plot(w1sorted,'*');
xlabel('Indices of the output');
ylabel('The output of the optimizer');
title('The histogram of the optimizer sorted output, lambda=0.01, mu=0.2');

figure;
plot(zt,'*');
xlabel('Indices of the output');
ylabel('The output of the optimizer');
title('The histogram of the optimizer sorted output, lambda=0.01, mu=0.2');

% Iterations for RFE
y1_predicted=normalmatrix*w1+b;
firstmse=(1/length(yy))*sum((yy-y1_predicted).^2);
[j2,j1]=size(normalmatrix);
denom=j2*(j2-1)/2;
su=0;

for u5=1:j2
    for u6=u5+1:j2
      if y1_predicted(u6)>y1_predicted(u5)
          su=su+1;
      end
    end
end

TCGA_con=su/denom

%%

myMSEi=zeros(10,1);%%save MSE error here
dis=zeros(100,1);%%save the distance between consecutive MSEs here
iter=1;%%number of iterrations
TCGA_con=zeros(100,1);%%%save TCGA concordance here
con9891=zeros(100,1);
con30161=zeros(100,1);
con17260=zeros(100,1);
nonzeronum=zeros(100,1);%%save the number of nonzero elements here
newA=normalmatrix;%%new gene expression matrices after faeture elimination
weight=zeros(12249,100);
com_genes=cell(12249,100);
comcom=common_genes;
bb=zeros(100,1);
bb(1,1)=b;
weight(:,1)=w1;

while dis(iter,1)<1
iter
zt=w1;
mean1=mean(abs(zt));
[pp1,pp2]=find(abs(zt)<=10^(-3));%%zt has dimensions like features*1,pp1 is the index for features
i=2;
while length(pp1) < 1
    i=i-1;
    [pp1,pp2]=find(abs(zt)<=mean1*10^(-i));
end
zt(pp1,1)=0;
nonzeronum(iter,1)=length(find(zt~=0)); 
%%concordance of training
[j2,j1]=size(newA);%%522*12249
yy_predicted=newA*zt+bb(iter,1);
clear b;
denom=j2*(j2-1)/2;
su=0;

for u5=1:j2
    for u6=u5+1:j2
      if yy_predicted(u6)>yy_predicted(u5)
          su=su+1;
      end
    end
end

TCGA_con(iter,1)=su/denom; 
%%Make a criterion for the performance of the method (COST FUNCTION)
myMSEi(iter,1)=(1/length(yy))*sum((yy-yy_predicted).^2);
zt(pp1,:)=[];

comcom(pp1,:)=[];%%common genes

absz=abs(zt);
[ztsort,oldindex1]=sort(absz,'descend');
removal=round(0.1*length(zt));
oldindex1((length(oldindex1)-removal):length(oldindex1),:)=[];
[u11,u22]=size(newA);
newA2=zeros(u11,length(oldindex1));
size(newA2)
for u3=1:length(oldindex1)
    newA2(:,u3)=newA(:,oldindex1(u3));
end
weight(1:length(zt),iter+1)=zt;
lambda=0.01;
mu=0.5;
[w,b] = optimizer(newA2,yy,length(oldindex1),lambda,mu);%%%lagrange
bb(iter+1,1)=b;
if iter==1
    dis(iter,1)=myMSEi(iter,1);
else
    dis(iter,:)=abs(myMSEi(iter,1)-myMSEi(iter-1,1));
end
w1=w;
clear pp1 pp2;
newA=newA2;
if dis(iter,1)>1
    dis(iter+1,1)=2;
end

iter=iter+1;
end
% [l1,l2]=find(dis>1);
% dis(l1,:)=[];
% [l1,l2]=find(dis==0);
% dis(l1,:)=[];
% [M1,N1]=size(dis);
% weight(:,M1+2:end)=[];
% [M2,N2]=size(weight);

%% First iterations

% zt=weight(1,1);
% b=bb(1,1);

%% Testing data GSE9891 
test_9891time=xlsread('9891time.xlsx');
%%%%%%%switching NA values to "zero"
[m11,n11]=find(test_9891time==0);
test_9891exp= xlsread('9891genes.xlsx');
copy1=test_9891exp;
copy1(:,m11)=[];%%19816*273
[u11,u21]=size(copy1);
%%Normalization-zscore on each 
y1=test_9891time;%%time recurrence of a tumor
yy1=y1;
yy1(m11)=[];%time recurrence matrix without zero or NaN values
%%delete unwanted genes (not common)

[o1,o2]=ismember(common_genes,GSE9891_genes);%% o2 shows the common_genes indices in GSE9891
gx_1=copy1(o2,:);%%instead of gx1
gx_1=normalize(gx_1',u21);
A1=zeros(length(zt),u21);
[yy1,oldindex1]=sort(yy1,'ascend');
gx_1=gx_1';
for u3=1:u21
A1(:,u3)=gx_1(:,oldindex1(u3));
end

yy1_predicted=A1'*zt+b(1,1);%%273
denom1=u21*(u21-1)/2;
su1=0;
yyy1=yy1_predicted;
for u5=1:u21
     for u6=u5+1:u21
       if yy1_predicted(u6)>yy1_predicted(u5)
           su1=su1+1;
       end
     end
 end
 total_concordance_GSE9891=su1/denom1
     
%% Testing data GSE17260

test_17260time=xlsread('17260time.xlsx');
% %%%%%%%switching NA values to "zero"
[m22,n22]=find(test_17260time==0);
test_17260exp= xlsread('17260genes.xlsx');
copy2=test_17260exp;
copy2(:,m22)=[];%%19816*273
[u11,u21]=size(copy2);
y2=test_17260time;%%time recurrence of a tumor
yy2=y2;
yy2(m22)=[];%time recurrence matrix without zero or NaN values
%%delete unwanted genes (not common)

[o1,o2]=ismember(common_genes,GSE17260_genes);%% o2 shows the common_genes indices in GSE9891
gx_2=copy2(o2,:);%%instead of gx2
gx_2=normalize(gx_2',u21);
A2=zeros(length(zt),u21);
[yy2,oldindex2]=sort(yy2,'ascend');
gx_2=gx_2';
for u3=1:u21
A2(:,u3)=gx_2(:,oldindex2(u3));
end 
%yy2_predicted=A2'*weight(:,iter)+bb(iter,1);%%273
 yy2_predicted=A2'*zt+b(1,1);%%273
denom2=u21*(u21-1)/2;
su1=0;
yyy2=yy2_predicted;
for u5=1:u21
    for u6=u5+1:u21
       if yy2_predicted(u6)>yy2_predicted(u5)
           su1=su1+1;
       end
     end
 end
total_concordance_GSE17260=su1/denom2  

%% Testing data GSE30161

test_30161time=xlsread('30161time.xlsx');
%%%%%%%switching NA values to "zero"
[m22,n22]=find(test_30161time==0);
test_30161exp= xlsread('30161genes.xlsx');
copy3=test_30161exp;
copy3(:,m22)=[];%%19816*273
[u11,u21]=size(copy3);
y3=test_30161time;%%time recurrence of a tumor
yy3=y3;
yy3(m22)=[];%time recurrence matrix without zero or NaN values
%%delete unwanted genes (not common)
%[o1,o2]=ismember(common_genes,GSE30161_genes);%% o2 shows the common_genes indices in GSE9891
[o1,o2]=ismember(common_genes,GSE30161_genes);%% o2 shows the common_genes indices in GSE9891
gx_3=copy3(o2,:);%%instead of gx3
gx_3=normalize(gx_3',u21);
A3=zeros(length(zt),u21);
[yy3,oldindex3]=sort(yy3,'ascend');
gx_3=gx_3';
for u3=1:u21
A3(:,u3)=gx_3(:,oldindex3(u3));
end
yy3_predicted=A3'*zt+b(1,1);%%273
denom3=u21*(u21-1)/2;
su1=0;
yyy3=yy3_predicted;
for u5=1:u21
    for u6=u5+1:u21
      if yy3_predicted(u6)>yy3_predicted(u5)
          su1=su1+1;
      end
    end
end
 total_concordance_GSE30161=su1/denom3  
