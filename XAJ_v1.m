

%Modified by ZIYU DING in OU ; 2020/1/14

 %FR ����������أ� S����ˮ��ˮ����mm��

%  SM KG KI   CG CI  ��ʱ�γ߶��й�   


%1)  K%����ϵ��ʵ�������ǻ����ż��ڶ��ı�ġ�0.5 2��
%2  SM�����������������0 100����Ҫ����
%3  KG%����ˮ������ˮ��ĳ���ϵ����0 0.7��
%4   KI�������ĳ���ϵ��]��0 0.7��
%KG+KI=0.7����
%5   CG����ˮ�������ϵ��[0.98 0.998]�൱�ڻ���ʱ��Ϊ50-500��
%6  CI ���������������ϵ��[0 0.9]
%7  CS������ˮ����ϵ�� ���Ϊ����ʱ�γ���   ���ݳ����Ʋ�ȡֵ���Ϊ��0  1��
%8   WUM[0 40] 9  WLM[40 150]   10  WDM[40 150]
%11   IMP��͸ˮ����ı���������Ȼ�����ֵ��С,ԼΪ[0.01 0.1]
%12  B����ˮ��ˮ�������ߵķ�����0.4����  ���ڼ���ƽ��ǧ��0.2-0.3
%13  C�����ɢ��ϵ����0 0.3��
%14   EX�������ˮ��ˮ�������ߵķ���[1 2]
%15   L�ͺ�ʱ�䡾1 10���Լ�����

function F=XAJ_v1(varargin)  
global DATA  AREA WU WL WD FR S  Spin

% XAJ���°��������г���,���ڵ����κ��Ŵ��㷨����,Ҳ�����°���ģ�͵�Ԥ��
% XX�ǵ��õ��Ż�����
% fit ����Ŀ�꺯������ֵ
% dc������Ч��ϵ��.
% result��һ������,���ظ�ʽΪ[ʱ��,����,ʵ������,��������]; %ֻ�ܷ������һ�������ֵ
% 
%clear  original
% XX=variable(weizhi,:);
%XX=variable;
% ������ʼֵ ,WU,WL,WD,S FR  ����ʼֵ�����ù���һ�β����ʶ�����ȷ����
%WU=20;WL=50;WD=10;FR=0.7; S=2; AREA=13988;

%WU =oripara(weizhi-1,1);WL=oripara(weizhi-1,2);WD=oripara(weizhi-1,3);FR=oripara(weizhi-1,4);S=oripara(weizhi-1,5);


W=WU+WL+WD; 
%AREA=17910; %JS�����
U=AREA/3.6/24; 
%��������E,��ɢ������P,ʵ������QS

%TIME=DATA(:,1); 
P=DATA(:,1); 
EM=DATA(:,2); 
QS=DATA(:,3); 
 
QI0=0.3.*QS(1);      %��ʼ������ 
QG0=0.4.*QS(1); %��ʼ���¾���
%  ��������
x= varargin{1};
%[num,numvars]=size(x); 
% ���Ż�����



%K SM KG KI CG CI CS WUM WLM WDM IMP B C EX L 
% K%����ϵ��   SM;%����ˮ��ˮ������ KG%����ˮ������ˮ��ĳ���ϵ��
%A_K=XX(:,1); A_SM=XX(:,2); A_KG=XX(:,3); A_KI=XX(:,4); 
%A_CG=XX(:,5); A_CI=XX(:,6); A_CS=XX(:,7); A_WUM=XX(:,8); 
%A_WLM=XX(:,9); A_WDM=XX(:,10); A_IMP=XX(:,11); A_B=XX(:,12); 
%A_C=XX(:,13); A_EX=XX(:,14); A_L=XX(:,15); 

A_K=x(1); A_SM=x(2); A_KG=x(3); A_KI=x(4); 
A_CG=x(5); A_CI=x(6); A_CS=x(7); A_WUM=x(8); 
A_WLM=x(9); A_WDM=x(10); A_IMP=x(11); A_B=x(12); 
A_C=x(13); A_EX=x(14); A_L=x(15); 


A_WM=A_WUM+A_WLM+A_WDM; 
 
for I=1:1      %%%% %%% ��ÿ���������
 
    K=A_K(I);    % ԭ����KC=A_KC(I);

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
    PE=P-K.*EM;  %�޸�KCΪK
    clear QJ
    for T=1:M                %%  T��ʱ��Ϊ��λ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  
  
          %����Ϊ�������� % PE ���꡾  �Ѽ��������
        if PE(T)<0 
            R=0; 
        else 
            if  W>=WM   %�ס�������ˮ��
                A=WMM; % WMMΪ����������ˮ����
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
% ����Ϊ��������zhengfa��
        if PE(T)<0 
            if WU+PE(T)>0 
                EU=K*EM(T); % �ϲ�������
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
        
       
%��ˮԴ����%%%%%%%%%%%%%%%%%%%%%% ���棬�����������¾���  
        SMM=(1+EX).*SM; 
        if (PE(T)<=0)|(R<=0) 
            RS=0; 
            RG=S.*KG.*FR;   %  S �ܵ�����ˮ
            RSS=RG.*KI./KG; 
        else 
            X=FR; 
            FR=(R-PE(T).*IMP)./PE(T); %��ʱ�β����������
            S=X.*S./FR; 
            SS=S; 
            Q=R./FR;   %��ʱ��͸ˮ����ϵľ���
           % G=fix(Q./5)+1; 
            %Q=Q./G; 
            %KID=KI.^(1/G); 
            %KGD=KID.*KG./KI; 
            RS=0; 
            RG=0; 
            %RI=0; 
            RSS=0;% �Լ��ӵ�
           % for J=1:G    % ��������  �о�������  ɾ���޸�
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
                RSS=S.*KI.*FR; %�޸�RSS
                S=Q+SS-(RS+RSS+RG)./FR; %�޸�RSS
            
        end 
        OUT(T,:)=[RS,RSS,RG]; 
    end     %
%һ��ģ���������%%%%%%%%%%%%%%%%%%%%%%%%%% 
%��������
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
%%����ΪĿ�꺯������%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
%?����Ϊ(��Ч��ϵ��)ȷ����ϵ������%%%%%%%%%%%%%%%%%% 
   QS1=QS(Spin+1:end);QJ1=QJ(Spin+1:end);
    f1=sum( (QS1-QJ1').^2);
    f2=sum((QS1-mean(QS1).*ones(M-Spin,1)).^2);

    dq=1-f1/f2; 
    dc(I)=dq; 
   % result =[TIME,P,QS,QJ']; 
   result =[P,QS,QJ']; 
   
  oripara(I,:)=[WU WL WD FR S];
end  %?I %һ������������I
 
fit=-fit'; %%�Ŵ��㷨Ϊ�������ֵ,�ڴ˼Ӹ���. 
F=-dc'; 
 

 