%%%%%%%%%%%%%%% PVO的改进版本，PVO+二维直方图+动态阈值double-mapping
tic
clear all;
clc

% imgfile = ['J:\Matlab_Program\MATLAB\Kodak-gray\'];
imgfile = ['J:\Matlab_Program\MATLAB\My_Histogram\2D-PVO\TEST\'];
imgdir = dir([imgfile,'\*.bmp']);
fid=fopen('fileName.txt','wt');
performance = zeros(length(imgdir)*2,100);


for i_img = 7:7
    %%%%% input
    i_img
    img = 2*(i_img-1)+1;
    I = double(imread([imgfile,'\',imgdir(i_img).name]));
    [A B] = size(I);
    
    
    %%%%% 定义嵌入容量和嵌入失真矩阵
    ECmat = zeros(4,4);
    EDmat = 2*ones(4,4);
    ECmat(1,1) = 1;
    ECmat(2,1) = 1;
    ECmat(2,2) = log2(3);
    ECmat(3,3) = 1;
    for i = 3:4
        ECmat(2,i) = 1;
        ECmat(i,2) = 1;
    end
    EDmat(1,1) = 0.5;
    EDmat(2,1) = 0.5;
    EDmat(2,2) = 2/3;
    EDmat(3,3) = 1;
    for i = 3:4
        EDmat(2,i) = 1.5;
        EDmat(i,2) = 1.5;
        EDmat(1,i) = 1;
        EDmat(i,1) = 1;
    end
    EDmat(1,2) = 1;
    
    %%%%% 主循环
    pIndex = 1;
    for PS = 10000:1000:10000
        R  = zeros(6,1024*15);
        RR = zeros(6,1024^2*15);
        
        index = 0;
        indexbis = 0;
        for a = 2:5
            for b = 2:5
%                 if a*b == 6 || a*b == 9 || a*b == 16 || a*b == 25
%                 else
%                     continue
%                 end
                if a*b > 4
%                     [a b]
                    [LM,bin_LM_len,I] = LocationMap(I,a,b);
                    %%%%% 计算复杂度
                    Capacity = PS + bin_LM_len+34;
                    nLM = 0;
                    NL = zeros(floor((A-2)/a),floor((B-2)/b));
                    for i = 1:floor((A-2)/a)
                        for j = 1:floor((B-2)/b)
                            nLM = nLM + 1;
                            if LM(nLM) == 1
                                continue
                            end
                            for ii = 1:a+1
                                for jj = 1:b+2
                                    if ii == a+1 || jj == b+1 || jj == b+2
                                        NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii+1,b*(j-1)+jj));
                                    end
                                end
                            end
                            for ii = 1:a+2
                                for jj = 1:b+1
                                    if ii == a+1 || ii == a+2 || jj == b+1
                                        NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii,b*(j-1)+jj+1));
                                    end
                                end
                            end
                        end
                    end
                    %%%%% 模拟嵌入 tri -- Peng
                    Hmax = zeros(4,4,1024);
                    Hmin = zeros(4,4,1024);
                    Hbismax = zeros(2,1024);
                    Hbismin = zeros(2,1024);
                    Htrimax = zeros(4,1024);
                    Htrimin = zeros(4,1024);
                    for i = 1:floor((A-2)/a)
                        for j = 1:floor((B-2)/b)
                            X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
                            X = X(:);
                            [Y In] = sort(X);
                            u = min(In(a*b),In(a*b-1));
                            v = max(In(a*b),In(a*b-1));
                            d1 = X(u)-Y(a*b-2);
                            d2 = X(v)-Y(a*b-2);
                            if d1 > 0 && d2 > 0
                                if d1-d2 == 0 || d1-d2 == 1
                                    Hbismax(1,NL(i,j)+1:1024) = Hbismax(1,NL(i,j)+1:1024)+1;
                                else
                                    Hbismax(2,NL(i,j)+1:1024) = Hbismax(2,NL(i,j)+1:1024)+1;
                                end
                            else
                                if d1-d2 == 0
                                    Htrimax(1,NL(i,j)+1:1024) = Htrimax(1,NL(i,j)+1:1024)+1;
                                else
                                    if d1-d2 == 1
                                        Htrimax(2,NL(i,j)+1:1024) = Htrimax(2,NL(i,j)+1:1024)+1;
                                    else
                                        if d1-d2 == 2
                                            Htrimax(3,NL(i,j)+1:1024) = Htrimax(3,NL(i,j)+1:1024)+1;
                                        else
                                            Htrimax(4,NL(i,j)+1:1024) = Htrimax(4,NL(i,j)+1:1024)+1;
                                        end
                                    end
                                end
                            end
                            d1 = min(d1,3);
                            d2 = min(d2,3);
                            Hmax(d1+1,d2+1,NL(i,j)+1:1024) = Hmax(d1+1,d2+1,NL(i,j)+1:1024)+1;
                            u = max(In(2),In(1));
                            v = min(In(2),In(1));
                            d1 = Y(3)-X(u);
                            d2 = Y(3)-X(v);
                            if d1 > 0 && d2 > 0
                                if d1-d2 == 0 || d1-d2 == 1
                                    Hbismin(1,NL(i,j)+1:1024) = Hbismin(1,NL(i,j)+1:1024)+1;
                                else
                                    Hbismin(2,NL(i,j)+1:1024) = Hbismin(2,NL(i,j)+1:1024)+1;
                                end
                            else
                                if d1-d2 == 0
                                    Htrimin(1,NL(i,j)+1:1024) = Htrimin(1,NL(i,j)+1:1024)+1;
                                else
                                    if d1-d2 == 1
                                        Htrimin(2,NL(i,j)+1:1024) = Htrimin(2,NL(i,j)+1:1024)+1;
                                    else
                                        if d1-d2 == 2
                                            Htrimin(3,NL(i,j)+1:1024) = Htrimin(3,NL(i,j)+1:1024)+1;
                                        else
                                            Htrimin(4,NL(i,j)+1:1024) = Htrimin(4,NL(i,j)+1:1024)+1;
                                        end
                                    end
                                end
                            end
                            d1 = min(d1,3);
                            d2 = min(d2,3);
                            Hmin(d1+1,d2+1,NL(i,j)+1:1024) = Hmin(d1+1,d2+1,NL(i,j)+1:1024)+1;
                        end
                    end
%                     Hmax(:,:,1024)
%                     Hmin(:,:,1024)
%                     Htrimax(:,1024)'
%                     Htrimin(:,1024)'
                    %%%%% 计算嵌入结果
                    H = Hmax+Hmin;
                    Hbis = Hbismax+Hbismin;
                    Htri = Htrimax+Htrimin;
                    for Scale = 1:1024-1
                        index = index+1;
                        EC = 0;
                        ED = 0;
                        for i = 1:4
                            for j = 1:4
                                EC = EC+ECmat(i,j)*H(i,j,Scale);
                                ED = ED+EDmat(i,j)*H(i,j,Scale);
                            end
                        end
                        R(1,index) = floor(EC)-4-12-18;
                        R(2,index) = ED;
                        if H(3,1,Scale) > H(1,1,Scale)
                            d = H(3,1,Scale) - H(1,1,Scale);
                            EC = EC+d;
                            ED = ED-d/2;
                            R(1,index) = floor(EC)-4-12-18;
                            R(2,index) = ED;
                        end
                        for Scalebis = Scale+1:1024
                            indexbis = indexbis+1;
                            ECbis = Hbis(1,Scalebis)-Hbis(1,Scale);
                            EDbis = Hbis(2,Scalebis)-Hbis(2,Scale)+ECbis/2;
                            ECtri = Htri(2,Scalebis)-Htri(2,Scale)+max(Htri(1,Scalebis)-Htri(1,Scale),Htri(3,Scalebis)-Htri(3,Scale));
                            EDtri = sum(Htri(:,Scalebis))-sum(Htri(:,Scale))-ECtri/2;
                            RR(1,indexbis) = R(1,index)+ECbis+ECtri;
                            RR(2,indexbis) = R(2,index)+EDbis+EDtri;
                            RR(2,indexbis) = 10*log10(A*B*255^2/RR(2,indexbis));
                            RR(3,indexbis) = Scale;
                            RR(4,indexbis) = Scalebis;
                            RR(5,indexbis) = a;
                            RR(6,indexbis) = b;
                        end
                        R(2,index) = 10*log10(A*B*255^2/ED);
                        R(3,indexbis) = Scale;
                        R(4,indexbis) = Scale;
                        R(5,indexbis) = a;
                        R(6,indexbis) = b;
                        if EC >= Capacity
                            T = Scale;
                            break
                        end
                    end
                end
            end
        end
        %%%%% 实验结果
        RRR = zeros(6,indexbis);
        k = 0;
        for i = 1:indexbis
            if RR(1,i) > 10000 && RR(2,i) > 0
                k = k+1;
                RRR(:,k) = RR(:,i);
            end
        end
        RR = RRR(:,1:k);
        RRR = zeros(6,indexbis);
        k = 0;
        for i = 1:1024*15
            if R(1,i) > 10000 && R(2,i) > 0
                k = k+1;
                RRR(:,k) = R(:,i);
            end
        end
        R = RRR(:,1:k);
        %     plot(RR(1,:),RR(2,:),'r.')
        %     plot(R(1,:),R(2,:),'k.')
        [~,ind] = sort(RR(2,:),'descend');
        RR(1,:) = RR(1,ind);
        RR(2,:) = RR(2,ind);
        RR(3,:) = RR(3,ind);
        RR(4,:) = RR(4,ind);
        RR(5,:) = RR(5,ind);
        RR(6,:) = RR(6,ind);

        ind = find(RR(1,:)>=Capacity,1,'first');
        if isempty(ind)
            break
        end
        [RR(2,ind) RR(3,ind) RR(4,ind) RR(5,ind) RR(6,ind)]
        
        performance(img,pIndex) = RR(1,ind) - (bin_LM_len+34);
        performance(img+1,pIndex) = RR(2,ind);
        pIndex = pIndex + 1;
        
    end
    
    
    
    
end
% hold on
% save ProD_45.mat performance

toc