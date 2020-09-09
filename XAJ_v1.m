

%Modified by ZIYU DING in OU ; 2020/1/14

 %FR 产流面积比重； S自由水蓄水量（mm）

%  SM KG KI   CG CI  与时段尺度有关   


%1)  K%蒸发系数实际上上是会随着季节而改变的【0.5 2】
%2  SM表层土自由永容量【0 100】重要参数
%3  KG%地下水从自由水库的出流系数【0 0.7】
%4   KI攘中流的出流系数]【0 0.7】
%KG+KI=0.7左右
%5   CG地下水库的消退系数[0.98 0.998]相当于汇流时间为50-500天
%6  CI 深层攘中流的消退系数[0 0.9]
%7  CS河网蓄水消退系数 大概为单个时段长的   根据程序推测取值大概为【0  1】
%8   WUM[0 40] 9  WLM[40 150]   10  WDM[40 150]
%11   IMP不透水面积的比例。在天然流域此值很小,约为[0.01 0.1]
%12  B张力水蓄水容量曲线的方次在0.4左右  对于几百平方千米0.2-0.3
%13  C深层蒸散发系数【0 0.3】
%14   EX表层自由水蓄水容量曲线的方次[1 2]
%15   L滞后时间【1 10】自己估计

function F=XAJ_v1(varargin)  
global DATA  AREA WU WL WD FR S  Spin

% XAJ是新安江的运行程序,用于单纯形和遗传算法调用,也用于新安江模型的预报
% XX是调用的优化参数
% fit 返回目标函数的适值
% dc返回有效性系数.
% result是一个数组,返回格式为[时间,雨量,实测流量,计算流量]; %只能返回最后一组参数的值
% 
%clear  original
% XX=variable(weizhi,:);
%XX=variable;
% 输入起始值 ,WU,WL,WD,S FR  【起始值可以用过第一次参数率定重新确定】
%WU=20;WL=50;WD=10;FR=0.7; S=2; AREA=13988;

%WU =oripara(weizhi-1,1);WL=oripara(weizhi-1,2);WD=oripara(weizhi-1,3);FR=oripara(weizhi-1,4);S=oripara(weizhi-1,5);


W=WU+WL+WD; 
%AREA=17910; %JS的面积
U=AREA/3.6/24; 
%输入雨量E,蒸散发能力P,实测流量QS

%TIME=DATA(:,1); 
P=DATA(:,1); 
EM=DATA(:,2); 
QS=DATA(:,3); 
 
QI0=0.3.*QS(1);      %初始壤中流 
QG0=0.4.*QS(1); %初始地下径流
%  参数处理
x= varargin{1};
%[num,numvars]=size(x); 
% 所优化参数



%K SM KG KI CG CI CS WUM WLM WDM IMP B C EX L 
% K%蒸发系数   SM;%自由水蓄水库容量 KG%地下水从自由水库的出流系数
%A_K=XX(:,1); A_SM=XX(:,2); A_KG=XX(:,3); A_KI=XX(:,4); 
%A_CG=XX(:,5); A_CI=XX(:,6); A_CS=XX(:,7); A_WUM=XX(:,8); 
%A_WLM=XX(:,9); A_WDM=XX(:,10); A_IMP=XX(:,11); A_B=XX(:,12); 
%A_C=XX(:,13); A_EX=XX(:,14); A_L=XX(:,15); 

A_K=x(1); A_SM=x(2); A_KG=x(3); A_KI=x(4); 
A_CG=x(5); A_CI=x(6); A_CS=x(7); A_WUM=x(8); 
A_WLM=x(9); A_WDM=x(10); A_IMP=x(11); A_B=x(12); 
A_C=x(13); A_EX=x(14); A_L=x(15); 


A_WM=A_WUM+A_WLM+A_WDM; 
 
for I=1:1      %%%% %%% 对每组参数计算
 
    K=A_K(I);    % 原来是KC=A_KC(I);

    SM=A_SM(I); 
    KG=A_KG(I); 
    KI=A_KI(I); 
    CG=A_CG(I); 
    CI=A_CI(I); 
    CS=A_CS(I); 
    WUM=A_WUM(I); 
    WLM=A_WLM(I); 
    WDM=A_WDM(I); 
    WM=WUM+WLM+WDM; 
    IMP=A_IMP(I); 
    B=A_B(I); 
    C=A_C(I); 
    EX=A_EX(I); 
    L=A_L(I); 
    L=round(L);   
    
    
    WMM=(1+B).*WM/(1-IMP); 
     M=size(P,1);  
    PE=P-K.*EM;  %修改KC为K
    clear QJ
    for T=1:M                %%  T以时段为单位计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  
  
          %以下为产流计算 % PE 净雨【  已检查修正】
        if PE(T)<0 
            R=0; 
        else 
            if  W>=WM   %Ｗ　土壤含水量
                A=WMM; % WMM为流域最大点蓄水容量
            else 
                A=WMM*(1-(1-W/WM).^(1/(1+B))); 
            end 
 
            if A+PE(T)>0 
                if A+PE(T)<WMM 
                    R=PE(T)-WM+W+WM.*(1-(PE(T)+A)./WMM).^(1+B); 
                else 
                    R=PE(T)+W-WM; 
                end 
            else 
                R=0; 
            end 
        end 
% 以下为蒸发计算zhengfa　
        if PE(T)<0 
            if WU+PE(T)>0 
                EU=K*EM(T); % 上层蒸发量
                ED=0; 
                EL=0; 
                WU=WU+PE(T); 
            else 
                EU=WU+P(T); 
                WU=0; 
                if WL>C*WLM 
                    EL=(K.*EM(T)-EU).*WL/WLM; 
                    WL=WL-EL; 
                    ED=0; 
                else 
                    if WL>C.*(K.*EM(T)-EU) 
                        EL=C.*(K.*EM(T)-EU); 
                        WL=WL-EL; 
                        ED=0; 
                    else 
                        EL=WL; 
                        WL=0; 
                        ED=C.*(K*EM(T)-EU)-EL; 
                        WD=WD-ED; 
                    end 
                end 
            end 
        else 
            EU=K.*EM(T); 
            ED=0; 
            EL=0; 
            if WU+PE(T)-R<WUM 
                WU=WU+PE(T)-R; 
            else 
                if WU+WL+PE(T)-WUM>WLM 
                    WU=WUM; 
                    WL=WLM; 
                    WD=W+PE(T)-R-WU-WL; 
                else 
                    WU=WUM; 
                    WL=WU+WL+PE(T)-R-WUM; 
                end 
            end 
        end 
        E=EU+EL+ED; 
        W=WU+WL+WD; 
        
       
%分水源计算%%%%%%%%%%%%%%%%%%%%%% 地面，壤中流，地下径流  
        SMM=(1+EX).*SM; 
        if (PE(T)<=0)|(R<=0) 
            RS=0; 
            RG=S.*KG.*FR;   %  S 总的自由水
            RSS=RG.*KI./KG; 
        else 
            X=FR; 
            FR=(R-PE(T).*IMP)./PE(T); %本时段产流面积比例
            S=X.*S./FR; 
            SS=S; 
            Q=R./FR;   %本时段透水面积上的净雨
           % G=fix(Q./5)+1; 
            %Q=Q./G; 
            %KID=KI.^(1/G); 
            %KGD=KID.*KG./KI; 
            RS=0; 
            RG=0; 
            %RI=0; 
            RSS=0;% 自己加的
           % for J=1:G    % 分区计算  感觉有问题  删除修改
                if S>=SM 
                    AU=SMM; 
                else 
                    AU=SMM.*(1-(1-S./SM).^(1./(1+EX))); 
                end 
 
                if AU+Q<SMM 
                    RS=(Q-SM+S+SM.*(1-(Q+AU)./SMM).^(1+EX)).*FR; 
                else 
                    RS=(Q+S-SM).*FR; 
                end 
                S=Q-RS./FR+S; 
                RG=S.*KG.*FR; 
                RSS=S.*KI.*FR; %修改RSS
                S=Q+SS-(RS+RSS+RG)./FR; %修改RSS
            
        end 
        OUT(T,:)=[RS,RSS,RG]; 
    end     %
%一次模型运算结束%%%%%%%%%%%%%%%%%%%%%%%%%% 
%汇流计算
    RS=OUT(:,1); RSS=OUT(:,2);RG=OUT(:,3); 

       TRS(1)=RS(1).*U;
    TRSS(1)=QI0 ;
    TRG(1)=QG0 ;
    TR(1)=TRS(1)+TRSS(1)+TRG(1);
    for T=2:M
        TRS(T)=RS(T).*U;
        TRSS(T)=TRSS(T-1).*CI+RSS(T).*(1-CI).*U;
        TRG(T)=TRG(T-1).*CG+RG(T).*(1-CG).*U;
        TR(T)=TRS(T)+TRSS(T)+TRG(T);
    end
    QJ=TR;
    if L<0  L=0;end
    for T=L+2:M
        QJ(T)=CS.*QJ(T-1)+(1-CS).*TR(T-L);
    end
%%以下为目标函数计算%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    alf=0.6; 
    y1=0;y2=0; 
    n1=1;n2=1; 
    for T=1:M 
        if QJ(T)>800 
            y1=(QJ(T)-QS(T)).^2+y1; 
            n1=n1+1; 
        else 
            y2=(QJ(T)-QS(T)).^2+y2; 
            n2=n2+1; 
        end 
    end 
    q0=mean(QS); 
    q1=mean(QJ); 
    y=(y1*alf/n1+y2*(1-alf)/n2)*(1+abs(q0-q1)/q0); 
fit(I)=y; 
%?以下为(有效性系数)确定性系数计算%%%%%%%%%%%%%%%%%% 
   QS1=QS(Spin+1:end);QJ1=QJ(Spin+1:end);
    f1=sum( (QS1-QJ1').^2);
    f2=sum((QS1-mean(QS1).*ones(M-Spin,1)).^2);

    dq=1-f1/f2; 
    dc(I)=dq; 
   % result =[TIME,P,QS,QJ']; 
   result =[P,QS,QJ']; 
   
  oripara(I,:)=[WU WL WD FR S];
end  %?I %一组参数计算结束I
 
fit=-fit'; %%遗传算法为了求最大值,在此加负号. 
F=-dc'; 
 

 