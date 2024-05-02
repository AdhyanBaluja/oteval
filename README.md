% Simplex Method
%max z=2x1+5X2
%x1+4x2<=24
%3x1+1x2<=21
%x1+x2<=9
clc
clear all
format short
Noofvariables=2;
C=[2 5];
a=[1 4; 3 1; 1 1]
b=[24; 21; 9]
s=eye(size(a,1))
A=[a s b]
cost=zeros(1,size(A,2))
cost(1:Noofvariables)=C
bv= Noofvariables+1:1:size(A,2)-1
zjcj=cost(bv)*A-cost
zcj=[zjcj; A]
simptable=array2table(zcj);
simptable.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','s_3','sol'}
RUN=true;
while RUN
if any(zjcj<0); %check for (most) negative value
    fprintf(' the current BFS is not optimal \n')
   zc=zjcj(1:end-1);
   [Enter_val, pvt_col]= min(zc) 
   if all(A(:,pvt_col)<=0)
    error('LPP is Unbounded all enteries are <=0 in column %d',pvt_col);
   else
       sol=A(:,end)
       column=A(:,pvt_col)
       for i=1:size(A,1)
         if column(i)>0
            ratio(i)= sol(i)./column(i)
         else
            ratio(i)=inf
         end
       end
       [leaving_val, pvt_row]=min(ratio)
   end
bv(pvt_row)=pvt_col
pvt_key=A(pvt_row, pvt_col)
A(pvt_row,:)=A(pvt_row,:)./pvt_key
for i=1:size(A,1)
    if i~=pvt_row
        A(i,:)=A(i,:)-A(i, pvt_col).*A(pvt_row,:)
    end
end
    zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:)
    zcj=[zjcj;A]
    table=array2table(zcj)
    table.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','s_3','sol'}
else
    RUN=false;
    fprintf('The current BFS is optimal \n')
end
end


%BigM method
format short
clear all
clc
% Cost=[-4 -5 0 0 -1000 -1000 0]
% A=[3 1 1 0 0 0 27; 3 2 0 -1 1 0 3; 5 5 0 0 0 1 60]
% % BV=[3 5 6]
Cost=[-2 -1 0 0 -10000 -10000 0]
A=[3 1 0 0 1 0 3; 4 3 -1 0 0 1 6 ;1 2 0 1 0 0 3]
BV=[5 6 4]
 
ZjCj=Cost(BV)*A-Cost
 zcj=[Cost;ZjCj;A];
    bigmtable=array2table(zcj);
    bigmtable.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','A_1','A_2','sol'}
 
RUN= true;
while RUN
    ZC=ZjCj(1:end-1)
    if any(ZC<0)
        fprintf('  The current BFS is not optimal\n')
        [ent_col,pvt_col]=min(ZC)
        fprintf('Entering Col =%d \n' , pvt_col);
        sol=A(:,end)
        Column=A(:,pvt_col)
        if Column<=0
            error('LPP is unbounded');
        else
            for i=1:size(A,1)
                if Column(i)>0
                    ratio(i)=sol(i)./Column(i)
                else
                    ratio(i)=inf
                end
            end
            [MinRatio,pvt_row]=min(ratio)
            fprintf('leaving Row=%d \n', pvt_row);
        end
        BV(pvt_row)=pvt_col;
        pvt_key=A(pvt_row,pvt_col);
        A(pvt_row,:)=A(pvt_row,:)./ pvt_key;
        for i=1:size(A,1)
            if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
            end
        end
        ZjCj=ZjCj-ZjCj(pvt_col).*A(pvt_row,:)
        ZCj=[ZjCj;A]
        TABLE=array2table(ZCj);
        TABLE.Properties.VariableNames(1:size(ZCj,2))={'x_1','x_2','s_1','s_2','A_1','A_2','sol'}
    else
        RUN=false;
        fprintf('  Current BFS is Optimal \n');
    end
end


%Least Cost Method
clc
clear all
format short
% Matlab Code of Least Cost Method (LCM)
% Input Information
%% Input Phase
Cost=[11 20 7 8; 21 16 10 12; 8 12 18 9]
A=[50 40 70]
B=[30 25 35 40]
%% To check unbalanced/balanced Problem
if sum(A)==sum(B)
    fprintf('Given Transportation Problem is Balanced \n')
else
   fprintf('Given Transportation Problem is Unbalanced \n') 
   if sum(A)<sum(B)
       Cost(end+1,:)=zeros(1,size(B,2))
       A(end+1)=sum(B)-sum(A)
   elseif sum(B)<sum(A)
   Cost(:,end+1)=zeros(1,size(A,2))
       B(end+1)=sum(A)-sum(B)  
   end
end
ICost=Cost
X=zeros(size(Cost))   % Initialize allocation
[m,n]=size(Cost)      % Finding No. of rows and columns
BFS=m+n-1             % Total No. of BFS
%% Finding the cell(with minimum cost) for the allocations
for i=1:size(Cost,1)
    for j=1:size(Cost,2)
hh=min(Cost(:))   % Finding minimum cost value
[Row_index, Col_index]=find(hh==Cost)  % Finding position of minimum cost cell
x11=min(A(Row_index),B(Col_index))
[Value,index]=max(x11)            % Find maximum allocation
ii=Row_index(index)       % Identify Row Position
jj=Col_index(index)        % Identify Column Position
y11=min(A(ii),B(jj))        % Find the value
X(ii,jj)=y11
A(ii)=A(ii)-y11
B(jj)=B(jj)-y11
Cost(ii,jj)=inf
    end
end
%% Print the initial BFS
fprintf('Initial BFS =\n')
IBFS=array2table(X)
disp(IBFS)
%% Check for Degenerate and Non Degenerate
TotalBFS=length(nonzeros(X))
if TotalBFS==BFS
    fprintf('Initial BFS is Non-Degenerate \n')
else
    fprintf('Initial BFS is Degenerate \n')
end
%% Compute the Initial Transportation cost
InitialCost=sum(sum(ICost.*X))
fprintf('Initial BFS Cost is = %d \n',InitialCost)



%Steepest descent method
clc;
syms x1 x2;
f = x1^2 -x1*x2 +x2^2;
f1=inline(f);
g=gradient(f);
g1=inline(g);
hes=hessian(f);
hes1=inline(hes);
x0=[1,1/2];
x=[ ];
miter=2;
iter=0;
tol=0.05;
while (f1(x0(1),x0(2)))>tol && iter<miter
    x=[x;x0];
    s=-g1(x0(1),x0(2));
    h = hes1(x0)
    lamb=s'*s./(s'*h*s);
    xnew=x0+lamb.*s';
    x0=xnew;
    iter=iter+1;
end
fprintf('optimal sol %f,%f\n',x0(1),x0(2));
fprintf('optial ol value %f',f1(x0(1),x0(2)));


