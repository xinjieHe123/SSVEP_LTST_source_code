%% transfer_cross_subject--叠加性质
clear;
clc;
%% 获取数据
choose = 2; % choose = 2 beta 数据集；choose = 1;benchmark 数据集
method = 1; % CSSFT方法 -- method = 1 所提方法和FBCCA；method = 2 为CSSFT方法；method = 3 为TT-CCA方法；
% for choose = 1:2
if choose == 1
% 获取处理之后的数据
data = importdata("eeg_data.mat"); 
% 数据形式为：通道 * 采样点 * 子代 * 目标数 * block * subject
channel = size(data,1);
samplePoint =  size(data,2);
% subband = size(data,3);
subband = 1;
FB_coef0=[1:subband].^(-1.25)+0.25; % for filter bank analysis

target = size(data,4);
block = size(data,5);
FB_coef=FB_coef0'*ones(1,length(target));
subject = size(data,6);
source_subject = 1:subject;

%% 相关参数
pha_val=[0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 ...
        0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5]*pi;
sti_f=[8.0:1:15.0, 8.2:1:15.2,8.4:1:15.4,8.6:1:15.6,8.8:1:15.8];
tw_1 = 0.6 : 0.1 : 1.5;
Fs = 250;
simulation_target = 40;
Nh = 5;
% target_templete = sincostemplate(sti_f,pha_val,tw,Fs,simulation_target,Nh);
for TW = 1:length(tw_1)
    tw = tw_1(TW);
    target_templete = sincostemplate(sti_f,pha_val,tw,Fs,simulation_target,Nh);
%% 划分源域和目标域，留一交叉验证
for cross = 1:subject
    source_subject = 1:subject;
    target_subject = cross;
    source_subject(cross)= [];
    source_data = data(:,:,:,:,:,source_subject);
    % source_data = source_data_1(:,1:Fs * tw,:,:,:,:);
    if (method == 2)
        % 进行CSSFT方法复现
        % 被试选择（此方法不可行，忽略），利用所有源域被试数据
        % 获取FBCCA识别准确率最高的被试作为源域被试
        addpath('D:\个人信息\个人论文\cross_subject_transfer');
        result_benchmark = xlsread("result_benchmark.xlsx");
        if (TW > 5)
            % 获取1.1s到1.6s的数据结果
            FBCCA_result = result_benchmark(40:74,5*(TW-6)+2);
        else
            FBCCA_result = result_benchmark(1:35,5*(TW-1)+2);
        end

        % 选择被试最好的源域结果
        FBCCA_result(target_subject) = [];
        [~,index] =  max(FBCCA_result);
        % 获取index号源域被试数据
        souece_choose_data = squeeze(source_data(:,:,:,:,:,index));
        % 获取模板，进行大平均
        % average_template = squeeze(mean(source_data,6));
        all_average_template = squeeze(mean(souece_choose_data,5));
        % 获取每个子代和每个目标下的迁移空间滤波器
        % 当前时间窗,训练空间滤波器
        for sub = 1:subband 
            for tg = 1:target
                % 获取模板数据
                template_data = all_average_template(:,1:Fs * tw,sub,tg);
                template_reference = target_templete(10 * (tg-1)+1:10 * (tg-1)+10,:);
                % 模板之间进行典型相关分析，获取空间滤波器参数
                [u1,v1,r] = canoncorr(template_data',template_reference');
                % 存储
                U(sub,tg,:,:) = u1(:,1);
                V(sub,tg,:,:) = v1(:,1);
            end
        end
        % 测试阶段
        test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
        test_data = test_data_1(:,1:Fs*tw,:,:,:,:);

        count_1 = 0;
        for bk = 1:block
            for tg = 1:target
            % count = count + 1;
            tic
            test_test = squeeze(test_data(:,:,:,tg,bk));
                for sub = 1:subband
                    for tt = 1:target
                        sub_test = test_data(:,:,sub,tg,bk);
                        % CCA
                        template_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                        template = template_1(:,1:tw * Fs);
                        [~,~,r] = canoncorr(sub_test',template');
                        % CCA(sub,tt) = r(1);
                        r_2 = r(1);
                        % 迁移
                        % 获取空间滤波器
                        u_filter = squeeze(U(sub,tt,:,:))';
                        v_filter = squeeze(V(sub,tt,:,:))';
                        template_reference_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                        r_1 = corrcoef(u_filter * sub_test,v_filter * template_reference_1);
                        % r_1(sub,tt) = abs(r_11(1,2));
                        R(sub,tt) = sign(r_2) * r_2^2 + sign(r_1(1,2)) *  r_1(1,2)^2;
                    end
                end
                % r_21 = sum((r_2).*FB_coef,1);
                % r_11 = sum((r_1).*FB_coef,1);
                % R1 = r_21./norm(r_21) + r_11./norm(r_11);
                R1 = sum((R).*FB_coef,1);
                [~,idx] =  max(R1);
                if idx == tg
                    count_1 = count_1 + 1;
                end
            end
        end
        accurary_2(TW,cross) = count_1 / (bk * tg);
        a =1;
    elseif(method == 3)
        % tt-CCA + adaptive update
        accurary_TT(TW,cross) = tt_CCA(source_data,target_templete,Fs,tw,data,target_subject,FB_coef);
    else
    % 对于源域被试，每个被试模板平均化，去除相关脑电噪声
    ss_data = mean(source_data,5);
    % 叠加原理，卷积分解，获得脉冲序列信号
    for ss = 1:length(source_subject)
        for sub = 1:subband
            for tg = 1:target
                eeg_data_average = ss_data(:,:,sub,tg,1,ss);
                dataLength = size(eeg_data_average,2);
                % 获取每个子代下的周期信号
                fs_0=sti_f(tg);
                ph_0=pha_val(tg);
                frequency_period=1.05*1/sti_f(tg);
                % 生成周期信号成分H
                [H0,h0]=my_conv_H(fs_0,ph_0,Fs,dataLength/Fs,60,frequency_period);
                
                h_len=size(H0,1);
                % simultaneously
                ssvep0 = eeg_data_average;
                w0_old=randn(1,size(eeg_data_average,1));
                x_hat_old=w0_old*ssvep0*H0'*inv(H0*H0');
                e_old=norm(w0_old*ssvep0-x_hat_old*H0);
                iter_err=100;
                iter=1;
                while (iter_err(iter)>0.0001 && iter<200)
                    w0_new=x_hat_old*H0*ssvep0'*inv(ssvep0*ssvep0');
                    x_hat_new=w0_new*ssvep0*H0'*inv(H0*H0');
                    e_new=norm(w0_new*ssvep0-x_hat_new*H0);
                    iter=iter+1;
                    iter_err(iter)=abs(e_old-e_new);
                    w0_old=w0_new;
                    w0_old=w0_old/std(w0_old);
                    x_hat_old=x_hat_new;
                    x_hat_old=x_hat_old/std(x_hat_old);
                    e_old=e_new;
                end
                x_hat_=x_hat_new;
                x_hat=x_hat_(1:h_len);

                % Reconstructed SSVEP
                y_re=x_hat*H0;
                
                
                y_=w0_new*ssvep0;
                y_re(:,1:length(find(y_re==0)))=0.8*y_re(:,Fs+1:Fs+length(find(y_re==0)));
                 
                % figure;
                % subplot(2,1,1)
                % 时域信号
                % plot(1:length(y_re),y_re,'b-',LineWidth=1);
                % hold on
                % plot(1:length(y_),y_,'k-',LineWidth=1);
                % legend('Transfer superposition tempalte','Source average tempalte');
                % xlabel('Time (s)' )
                % ylabel('Amplitude');
                % ylim([-5,7])
                % title('corrcoef == 0.8454')
                % 频域信号
                % subplot(2,1,2)
                [fft_signal_1,f]=fft_5F(y_re,Fs,5);
                % plot(f,fft_signal_1,'b-',LineWidth=1);
                % hold on 
                [fft_signal,f]=fft_5F(y_,Fs,5);
                % plot(f,fft_signal,'k-',LineWidth=1);
                % hold on
                % xline(8,'-.m',LineWidth=1);
                % hold on
                % xline(16,'-.m',LineWidth=1);
                % hold on
                % xline(24,'-.m',LineWidth=1);
                % legend('Transfer superposition tempalte','Source average tempalte','target frequency');
                % xlabel('Time (s)' )
                % ylabel('Amplitude');
                % title('corrcoef == 0.9083')
                r_f = corrcoef(fft_signal,fft_signal_1);
                r_f_1(ss,sub,tg) = abs(r_f(1,2));
                r = corrcoef(y_,y_re);     % the similarity between the reconstructed and original ssvep
                r_1(ss,sub,tg) = abs(r(1,2));
                % ycor(sn,stim_no)=abs(r(1,2));
                y_storage(ss,sub,tg,:,:) = y_;
                y_re_storage(ss,sub,tg,:,:) = y_re;
                
                source(tg).x_hat(:,:,sub,ss) = x_hat;
                source(tg).w_new(:,:,sub,ss)= w0_new;
            end
            %% 计算两两之间的皮尔逊相关性系数
                % for ii = 1:size(y_storage,1)
                    % for jj = 1:size(y_re_storage,1)
                        % 时域
                        % C0RR_1 = corrcoef(squeeze(y_storage(ii,:,:)),squeeze(y_re_storage(jj,:,:)));
                        % CORR(ii,jj) = C0RR_1(1,2);
                        % 频域
                        % C0RR_1_F = corrcoef(squeeze(fft_signal_1(ii,:,:)),squeeze(fft_signal(jj,:,:)));
                        % CORR_F(ii,jj) = C0RR_1_F(1,2);
                       

                    % end
                % end
        end
    end

    %% 平均化所有源域被试模板，能否增强还是抑制信号？
    average_template_ = mean(y_storage,1);
    average_template_re = mean(y_re_storage,1);
    
    %% 单独模板

    
    % 获取不同时间窗下的测试数据
    test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
    test_data = test_data_1(:,1:Fs*tw,:,:,:,:);

    % 决策
   
    % 绘图观察
    % figure;
    % y_re_1 = squeeze(average_template_(1,sub,tg,:,:));
    % y_1 = squeeze(average_template_re(1,sub,tg,:,:));
    % [fft_signal_1,f]=fft_5F(y_re_1,Fs,5);
    % plot(f,fft_signal_1,'r-');
    % hold on
    % [fft_signal_1,f]=fft_5F(y_1,Fs,5);
    % plot(f,fft_signal_1,'b:');

    % 获取目标被试数据,所有源域(不做源域选择的情况)参与决策过程
    % test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
    % test_data = test_data_1(:,1:Fs*tw,:,:,:,:);
    % 策略1，采用平均模板+个体迁移滤波器进行融合决策
    % count = 0 ;
    % count_1 = 0;
    count_2 = 0;
    % count_3 = 0;
    count_4 = 0;
    count_5 = 0;
    count_6 = 0;
    count_7 = 0;
    for bk = 1:block
        for tg = 1:target
            % count = count + 1;
            tic
            test_test = squeeze(test_data(:,:,:,tg,bk));
            [subject_index,K] = subject_choose(test_test,target,subband,source,tw,Fs,FB_coef,source_subject,y_re_storage,target_templete);
            for tt = 1:target
                for sub = 1:subband
                    sub_test = test_data(:,:,sub,tg,bk);
                    %% 源域被试选择策略
                    % 基于FBCCA决策以获得预测标签值以及参考模板值
                    % [subject_index,K] = subject_choose(test_test,target,subband,source,tw,Fs,FB_coef,source_subject,y_re_storage,target_templete);
                    for ii = 1:length(subject_index)
                        source_transfer_filter_11 = source(tt).w_new(:,:,sub,subject_index(ii));
                        test_variable_11 = source_transfer_filter_11 * sub_test;
                        ss_tempalte_11 = squeeze(y_re_storage(subject_index(ii),sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        rr1_ss_choose = corrcoef(test_variable_11,ss_tempalte_11);
                        r_ss_choose(ii) = abs(rr1_ss_choose(1,2));
                    end
                    % 求和
                    SS_CCRS_choose = sum(r_ss_choose);
                    SS_CCRS1_choose(sub,tt) = SS_CCRS_choose;
                    %% 所有源域被试融合测试
                    % for ss = 1:length(source_subject)
                        % 获取源迁移空间滤波器
                        % source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        % test_variable = source_transfer_filter * sub_test;
                        % [fft_signal_1,f]=fft_5F(test_variable,Fs,5);
                        % plot(f,fft_signal_1,'r-');
                        % 获取模板
                        % template_sub_1 = squeeze(average_template_re(1,sub,tt,:,:))';
                        % template_sub = template_sub_1(1,1:Fs * tw);
                        % [fft_signal_1,f]=fft_5F(template_sub,Fs,5);
                        % hold on
                        % plot(f,fft_signal_1,'g-');
                        % r1a = corrcoef(test_variable,template_sub);
                        % [A,B,R_1] = canoncorr(sub_test',template_sub');
                        % R(ss) = abs(R_1);
                        % TSCORR_ss(ss) = abs(r1a(1,2));
                    % end
                    % TSCORR(sub,tt) = sum(TSCORR_ss);
                    % CRS(sub,tt) = sum(R);

                    %% 原型空间滤波器获取迁移空间滤波器共同成分+正余弦参考模板
                    % [CCRS(sub,tt),u1] = commen_filter(source_subject,source,target_templete,sub_test,tt,sub,tw,Fs);
                    %% 单独被试 + 共迁移空间滤波器
                    % for ss = 1:length(source_subject)
                        % source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        % test_variable = u1' * sub_test;
                        % ss_tempalte = squeeze(y_re_storage(ss,sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        % rr1_ss = corrcoef(test_variable,ss_tempalte);
                        % wr_ss(ss) = abs(rr1_ss(1,2));
                    
                    % end
                    % 求和
                    % WS_CCRS = sum(wr_ss);
                    % WS_CCRS1(sub,tt) = WS_CCRS;
                    %% 单独被试 + 单独迁移空间滤波器
                    for ss = 1:length(source_subject)
                        source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        test_variable = source_transfer_filter * sub_test;
                        ss_tempalte = squeeze(y_re_storage(ss,sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        rr1_ss = corrcoef(test_variable,ss_tempalte);
                        r_ss(ss) = abs(rr1_ss(1,2));
                    
                    end
                    % 求和
                    SS_CCRS = sum(r_ss);
                    SS_CCRS1(sub,tt) = SS_CCRS;
                    %% CCA 
                    template_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                    template = template_1(:,1:tw * Fs);
                    [~,~,r] = canoncorr(sub_test',template');
                    CCA(sub,tt) = r(1);

                end
            end
            % TSCORR1=sum((TSCORR).*FB_coef,1);
            % CRS1 = sum((CRS).*FB_coef,1);
            CCA1 = sum((CCA).*FB_coef,1);
            % CCRS1 = sum((CCRS).*FB_coef,1);
            SS_CCRS1 = sum((SS_CCRS1).*FB_coef,1);
            SS_CCRS1_choose1 = sum((SS_CCRS1_choose).*FB_coef,1);
            fusion_choose = CCA1/norm(CCA1) + SS_CCRS1_choose1/norm(SS_CCRS1_choose1);
            fusion_CCA_SS1 = CCA1/norm(CCA1) + SS_CCRS1/norm(SS_CCRS1);
            norm_CCA1 = CCA1/norm(CCA1);
            fusion_feature(bk,tg,:,:) = fusion_choose;
            norm_feature(bk,tg,:,:) = norm_CCA1;
            % WS_CCRS1 = sum((WS_CCRS1).*FB_coef,1);
            % 所提方法的结果
            [~,idx] =  max(CCA1);
            if idx == tg
                count_2 = count_2 + 1;
            end
            % 所提方法的结果
            % [~,idx] =  max(CRS1);
            % if idx == tg
                % count_1 = count_1 + 1;
            % end
            % 所提方法的结果
            % [~,idx] =  max(TSCORR1);
            % if idx == tg
                % count = count + 1;
            % end
            % 所提方法的结果
            % [~,idx] =  max(CCRS1);
            % if idx == tg
                % count_3 = count_3 + 1;
            % end

            [~,idx] =  max(SS_CCRS1);
            if idx == tg
                count_4 = count_4 + 1;
            end
            % 融合决策方法
            [~,idx] =  max(fusion_CCA_SS1);
            if idx == tg
                count_5 = count_5 + 1;
            end

            [~,idx] =  max(SS_CCRS1_choose1);
            if idx == tg
                count_6 = count_6 + 1;
            end

            [~,idx] =  max(fusion_choose);
            if idx == tg
                count_7 = count_7 + 1;
            end
            toc
            
        end
    end
    % accurary(cross) = count / (bk * tg)
    % accurary_1(cross) = count_1 / (bk * tg)
    accurary_2(TW,cross) = count_2 / (bk * tg)
    % accurary_3(cross) = count_3 / (bk * tg)
    accurary_4(TW,cross) = count_4 / (bk * tg)
    accurary_5(choose,TW,cross) = count_5 / (bk * tg)
    accurary_6(TW,cross) = count_6 / (bk * tg)
    accurary_7(choose,TW,cross) = count_7 / (bk * tg)
    end
end
end
else
    % 获取处理之后的数据
data = importdata("beta_eeg_data.mat");
% 数据形式为：通道 * 采样点 * 子代 * 目标数 * block * subject
channel = size(data,1);
samplePoint =  size(data,2);
% subband = size(data,3);
subband = 1;
FB_coef0=[1:subband].^(-1.25)+0.25; % for filter bank analysis

target = size(data,4);
block = size(data,5);
FB_coef=FB_coef0'*ones(1,length(target));
subject = size(data,6);
source_subject = 1:subject;

%% 相关参数
pha_val=[0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 ...
        0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5]*pi;
sti_f=[8.6:0.2:15.8,8.0 8.2 8.4];
tw_1 = 0.6 : 0.1 : 1.5;
Fs = 250;
simulation_target = 40;
Nh = 5;
% target_templete = sincostemplate(sti_f,pha_val,tw,Fs,simulation_target,Nh);
for TW = 1:length(tw_1)
    tw = tw_1(TW);
    target_templete = sincostemplate(sti_f,pha_val,tw,Fs,simulation_target,Nh);
%% 划分源域和目标域，留一交叉验证
for cross = 1:subject
    source_subject = 1:subject;
    target_subject = cross;
    source_subject(cross)= [];
    source_data = data(:,:,:,:,:,source_subject);
    % source_data = source_data_1(:,1:Fs * tw,:,:,:,:);
    if (method == 2)
        % 进行CSSFT方法复现
        % 被试选择（此方法不可行，忽略），利用所有源域被试数据
        % 获取FBCCA识别准确率最高的被试作为源域被试
        addpath('D:\个人信息\个人论文\cross_subject_transfer');
        result_benchmark = xlsread("result_beta.xlsx");
        if (TW > 5)
            % 获取1.1s到1.6s的数据结果
            FBCCA_result = result_benchmark(75:144,5*(TW-6)+2);
        else
            FBCCA_result = result_benchmark(1:70,5*(TW-1)+2);
        end

        % 选择被试最好的源域结果
        FBCCA_result(target_subject) = [];
        [~,index] =  max(FBCCA_result);
        % 获取index号源域被试数据
        souece_choose_data = squeeze(source_data(:,:,:,:,:,index));
        % 获取模板，进行大平均
        % average_template = squeeze(mean(source_data,6));
        all_average_template = squeeze(mean(souece_choose_data,5));
        % 获取每个子代和每个目标下的迁移空间滤波器
        % 当前时间窗,训练空间滤波器
        for sub = 1:subband 
            for tg = 1:target
                % 获取模板数据
                template_data = all_average_template(:,1:Fs * tw,sub,tg);
                template_reference = target_templete(10 * (tg-1)+1:10 * (tg-1)+10,:);
                % 模板之间进行典型相关分析，获取空间滤波器参数
                [u1,v1,r] = canoncorr(template_data',template_reference');
                % 存储
                U(sub,tg,:,:) = u1(:,1);
                V(sub,tg,:,:) = v1(:,1);
            end
        end
        % 测试阶段
        test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
        test_data = test_data_1(:,1:Fs*tw,:,:,:,:);

        count_1 = 0;
        for bk = 1:block
            for tg = 1:target
            % count = count + 1;
            tic
            test_test = squeeze(test_data(:,:,:,tg,bk));
                for sub = 1:subband
                    for tt = 1:target
                        sub_test = test_data(:,:,sub,tg,bk);
                        % CCA
                        template_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                        template = template_1(:,1:tw * Fs);
                        [~,~,r] = canoncorr(sub_test',template');
                        % CCA(sub,tt) = r(1);
                        r_2 = r(1);
                        % 迁移
                        % 获取空间滤波器
                        u_filter = squeeze(U(sub,tt,:,:))';
                        v_filter = squeeze(V(sub,tt,:,:))';
                        template_reference_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                        r_1 = corrcoef(u_filter * sub_test,v_filter * template_reference_1);
                        % r_1(sub,tt) = abs(r_11(1,2));
                        R(sub,tt) = sign(r_2) * r_2^2 + sign(r_1(1,2)) *  r_1(1,2)^2;
                    end
                end
                % r_21 = sum((r_2).*FB_coef,1);
                % r_11 = sum((r_1).*FB_coef,1);
                % R1 = r_21./norm(r_21) + r_11./norm(r_11);
                R1 = sum((R).*FB_coef,1);
                [~,idx] =  max(R1);
                if idx == tg
                    count_1 = count_1 + 1;
                end
            end
        end
        accurary_2(TW,cross) = count_1 / (bk * tg);
        a =1;
    elseif(method == 3)
        % tt-CCA + adaptive update
        accurary_TT(TW,cross) = tt_CCA(source_data,target_templete,Fs,tw,data,target_subject,FB_coef);
    else
    % 对于源域被试，每个被试模板平均化，去除相关脑电噪声
    ss_data = mean(source_data,5);
    % 叠加原理，卷积分解，获得脉冲序列信号
    for ss = 1:length(source_subject)
        for sub = 1:subband
            for tg = 1:target
                eeg_data_average = ss_data(:,:,sub,tg,1,ss);
                dataLength = size(eeg_data_average,2);
                % 获取每个子代下的周期信号
                fs_0=sti_f(tg);
                ph_0=pha_val(tg);
                frequency_period=1.05*1/sti_f(tg);
                % 生成周期信号成分H
                [H0,h0]=my_conv_H(fs_0,ph_0,Fs,dataLength/Fs,60,frequency_period);
                
                h_len=size(H0,1);
                % simultaneously
                ssvep0 = eeg_data_average;
                w0_old=randn(1,size(eeg_data_average,1));
                x_hat_old=w0_old*ssvep0*H0'*inv(H0*H0');
                e_old=norm(w0_old*ssvep0-x_hat_old*H0);
                iter_err=100;
                iter=1;
                while (iter_err(iter)>0.0001 && iter<200)
                    w0_new=x_hat_old*H0*ssvep0'*inv(ssvep0*ssvep0');
                    x_hat_new=w0_new*ssvep0*H0'*inv(H0*H0');
                    e_new=norm(w0_new*ssvep0-x_hat_new*H0);
                    iter=iter+1;
                    iter_err(iter)=abs(e_old-e_new);
                    w0_old=w0_new;
                    w0_old=w0_old/std(w0_old);
                    x_hat_old=x_hat_new;
                    x_hat_old=x_hat_old/std(x_hat_old);
                    e_old=e_new;
                end
                x_hat_=x_hat_new;
                x_hat=x_hat_(1:h_len);

                % Reconstructed SSVEP
                y_re=x_hat*H0;
                
                
                y_=w0_new*ssvep0;
                y_re(:,1:length(find(y_re==0)))=0.8*y_re(:,Fs+1:Fs+length(find(y_re==0)));
                 
                % figure;
                % 时域信号
                % plot(1:length(y_re),y_re);
                % hold on
                % plot(1:length(y_),y_);
                
                % 频域信号
                % figure;
                % [fft_signal_1,f]=fft_5F(y_re,Fs,5);
                % plot(f,fft_signal_1);
                % hold on 
                % [fft_signal,f]=fft_5F(y_,Fs,5);
                % plot(f,fft_signal);
                % r_f = corrcoef(fft_signal,fft_signal_1);
                r=corrcoef(y_,y_re);     % the similarity between the reconstructed and original ssvep
                % ycor(sn,stim_no)=abs(r(1,2));
                y_storage(ss,sub,tg,:,:) = y_;
                y_re_storage(ss,sub,tg,:,:) = y_re;
                
                source(tg).x_hat(:,:,sub,ss) = x_hat;
                source(tg).w_new(:,:,sub,ss)= w0_new;
            end
            %% 计算两两之间的皮尔逊相关性系数
                % for ii = 1:size(y_storage,1)
                    % for jj = 1:size(y_re_storage,1)
                        % 时域
                        % C0RR_1 = corrcoef(squeeze(y_storage(ii,:,:)),squeeze(y_re_storage(jj,:,:)));
                        % CORR(ii,jj) = C0RR_1(1,2);
                        % 频域
                        % C0RR_1_F = corrcoef(squeeze(fft_signal_1(ii,:,:)),squeeze(fft_signal(jj,:,:)));
                        % CORR_F(ii,jj) = C0RR_1_F(1,2);
                       

                    % end
                % end
        end
    end

    %% 平均化所有源域被试模板，能否增强还是抑制信号？
    average_template_ = mean(y_storage,1);
    average_template_re = mean(y_re_storage,1);
    
    %% 单独模板

    
    % 获取不同时间窗下的测试数据
    test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
    test_data = test_data_1(:,1:Fs*tw,:,:,:,:);

    % 决策
   
    % 绘图观察
    % figure;
    % y_re_1 = squeeze(average_template_(1,sub,tg,:,:));
    % y_1 = squeeze(average_template_re(1,sub,tg,:,:));
    % [fft_signal_1,f]=fft_5F(y_re_1,Fs,5);
    % plot(f,fft_signal_1,'r-');
    % hold on
    % [fft_signal_1,f]=fft_5F(y_1,Fs,5);
    % plot(f,fft_signal_1,'b:');

    % 获取目标被试数据,所有源域(不做源域选择的情况)参与决策过程
    % test_data_1 = squeeze(data(:,:,:,:,:,target_subject));
    % test_data = test_data_1(:,1:Fs*tw,:,:,:,:);
    % 策略1，采用平均模板+个体迁移滤波器进行融合决策
    % count = 0 ;
    % count_1 = 0;
    count_2 = 0;
    % count_3 = 0;
    count_4 = 0;
    count_5 = 0;
    count_6 = 0;
    count_7 = 0;
    for bk = 1:block
        for tg = 1:target
            % count = count + 1;
            tic
            test_test = squeeze(test_data(:,:,:,tg,bk));
            [subject_index,K] = subject_choose(test_test,target,subband,source,tw,Fs,FB_coef,source_subject,y_re_storage,target_templete);
            for tt = 1:target
                for sub = 1:subband
                    sub_test = test_data(:,:,sub,tg,bk);
                    %% 源域被试选择策略
                    % 基于FBCCA决策以获得预测标签值以及参考模板值
                    % [subject_index,K] = subject_choose(test_test,target,subband,source,tw,Fs,FB_coef,source_subject,y_re_storage,target_templete);
                    for ii = 1:length(subject_index)
                        source_transfer_filter_11 = source(tt).w_new(:,:,sub,subject_index(ii));
                        test_variable_11 = source_transfer_filter_11 * sub_test;
                        ss_tempalte_11 = squeeze(y_re_storage(subject_index(ii),sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        rr1_ss_choose = corrcoef(test_variable_11,ss_tempalte_11);
                        r_ss_choose(ii) = abs(rr1_ss_choose(1,2));
                    end
                    % 求和
                    SS_CCRS_choose = sum(r_ss_choose);
                    SS_CCRS1_choose(sub,tt) = SS_CCRS_choose;
                    %% 所有源域被试融合测试
                    % for ss = 1:length(source_subject)
                        % 获取源迁移空间滤波器
                        % source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        % test_variable = source_transfer_filter * sub_test;
                        % [fft_signal_1,f]=fft_5F(test_variable,Fs,5);
                        % plot(f,fft_signal_1,'r-');
                        % 获取模板
                        % template_sub_1 = squeeze(average_template_re(1,sub,tt,:,:))';
                        % template_sub = template_sub_1(1,1:Fs * tw);
                        % [fft_signal_1,f]=fft_5F(template_sub,Fs,5);
                        % hold on
                        % plot(f,fft_signal_1,'g-');
                        % r1a = corrcoef(test_variable,template_sub);
                        % [A,B,R_1] = canoncorr(sub_test',template_sub');
                        % R(ss) = abs(R_1);
                        % TSCORR_ss(ss) = abs(r1a(1,2));
                    % end
                    % TSCORR(sub,tt) = sum(TSCORR_ss);
                    % CRS(sub,tt) = sum(R);

                    %% 原型空间滤波器获取迁移空间滤波器共同成分+正余弦参考模板
                    % [CCRS(sub,tt),u1] = commen_filter(source_subject,source,target_templete,sub_test,tt,sub,tw,Fs);
                    %% 单独被试 + 共迁移空间滤波器
                    % for ss = 1:length(source_subject)
                        % source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        % test_variable = u1' * sub_test;
                        % ss_tempalte = squeeze(y_re_storage(ss,sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        % rr1_ss = corrcoef(test_variable,ss_tempalte);
                        % wr_ss(ss) = abs(rr1_ss(1,2));
                    
                    % end
                    % 求和
                    % WS_CCRS = sum(wr_ss);
                    % WS_CCRS1(sub,tt) = WS_CCRS;
                    %% 单独被试 + 单独迁移空间滤波器
                    for ss = 1:length(source_subject)
                        source_transfer_filter = source(tt).w_new(:,:,sub,ss);
                        test_variable = source_transfer_filter * sub_test;
                        ss_tempalte = squeeze(y_re_storage(ss,sub,tt,:,1:tw * Fs))';
                        % 计算被试中皮尔逊相关系数
                        rr1_ss = corrcoef(test_variable,ss_tempalte);
                        r_ss(ss) = abs(rr1_ss(1,2));
                    
                    end
                    % 求和
                    SS_CCRS = sum(r_ss);
                    SS_CCRS1(sub,tt) = SS_CCRS;
                    %% CCA 
                    template_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
                    template = template_1(:,1:tw * Fs);
                    [~,~,r] = canoncorr(sub_test',template');
                    CCA(sub,tt) = r(1);

                end
            end
            % TSCORR1=sum((TSCORR).*FB_coef,1);
            % CRS1 = sum((CRS).*FB_coef,1);
            CCA1 = sum((CCA).*FB_coef,1);
            norm_CCA1 = CCA1/norm(CCA1);
            % CCRS1 = sum((CCRS).*FB_coef,1);
            SS_CCRS1 = sum((SS_CCRS1).*FB_coef,1);
            SS_CCRS1_choose1 = sum((SS_CCRS1_choose).*FB_coef,1);
            fusion_choose = CCA1/norm(CCA1) + SS_CCRS1_choose1/norm(SS_CCRS1_choose1);
            fusion_CCA_SS1 = CCA1/norm(CCA1) + SS_CCRS1/norm(SS_CCRS1);
            fusion_feature_1(bk,tg,:,:) = fusion_choose;
            norm_feature_1(bk,tg,:,:) = norm_CCA1;
            % WS_CCRS1 = sum((WS_CCRS1).*FB_coef,1);
            % 所提方法的结果
            [~,idx] =  max(CCA1);
            if idx == tg
                count_2 = count_2 + 1;
            end
            % 所提方法的结果
            % [~,idx] =  max(CRS1);
            % if idx == tg
                % count_1 = count_1 + 1;
            % end
            % 所提方法的结果
            % [~,idx] =  max(TSCORR1);
            % if idx == tg
                % count = count + 1;
            % end
            % 所提方法的结果
            % [~,idx] =  max(CCRS1);
            % if idx == tg
                % count_3 = count_3 + 1;
            % end

            [~,idx] =  max(SS_CCRS1);
            if idx == tg
                count_4 = count_4 + 1;
            end
            % 融合决策方法
            [~,idx] =  max(fusion_CCA_SS1);
            if idx == tg
                count_5 = count_5 + 1;
            end

            [~,idx] =  max(SS_CCRS1_choose1);
            if idx == tg
                count_6 = count_6 + 1;
            end

            [~,idx] =  max(fusion_choose);
            if idx == tg
                count_7 = count_7 + 1;
            end
            toc
            
        end
    end
    % accurary(cross) = count / (bk * tg)
    % accurary_1(cross) = count_1 / (bk * tg)
    accurary_2(TW,cross) = count_2 / (bk * tg)
    % accurary_3(cross) = count_3 / (bk * tg)
    accurary_4(TW,cross) = count_4 / (bk * tg)
    accurary_5(TW,cross) = count_5 / (bk * tg)
    accurary_6(TW,cross) = count_6 / (bk * tg)
    accurary_7(TW,cross) = count_7 / (bk * tg)
    end
    toc
end

end

end
% end
accurary_2 = accurary_2';
accurary_4 = accurary_4';
accurary_5 = accurary_5';
accurary_6 = accurary_6';
accurary_7 = accurary_7';

mean_accurary_2 = mean(accurary_2,1);
mean_accurary_4 = mean(accurary_4,1);
mean_accurary_5 = mean(accurary_5,1);
mean_accurary_6 = mean(accurary_6,1);
mean_accurary_7 = mean(accurary_7,1);

%% 正余弦参考模板生成
function target_templete = sincostemplate(simulation_frequency,simulation_phase,datalength,fN,simulation_target,Nh)
target_templete = [];
for i = 1:simulation_target
      for k = 1:Nh
            for j = 1:(datalength * fN)
                  target_sample(2*k-1,j) = sin(2 * pi * k * simulation_frequency(1,i) * (1/fN * j) + simulation_phase(i));
                  target_sample(2*k,j) = cos(2 * pi * k * simulation_frequency(1,i) * (1/fN * j) + simulation_phase(i));
            end
      end
      target_templete = [target_templete;target_sample];
end
end
%% 5倍频傅里叶变换
function [fft_signal,f]=fft_5F(signal,F,num)
NFFT=num*F;
Y=fft(signal,NFFT);

P2=abs(Y/length(Y));

P1=P2(1:length(Y)/2+1);
P1(2:end-1) = P1(2:end-1);
P1(1) = 0;
fft_signal=P1;
f=F/2*linspace(0,1,NFFT/2+1); 
end

%% 共空间迁移滤波器
function [CCRS,u1] = commen_filter(source_subject,source,target_templete,sub_test,tt,sub,tw,Fs)
for ss = 1:length(source_subject)
    % 获取源迁移空间滤波器
    source_fliter_1 = source(tt).w_new(:,:,sub,ss);
    % 计算自协方差矩阵
    source_fliter = source_fliter_1/norm(source_fliter_1);
    Cxx(ss,:,:) = source_fliter' * source_fliter;
end
% 求和
C_xx_sum = squeeze(sum(Cxx,1));
[V,D] = eig(C_xx_sum);
[~, loc] = max(diag(D));
u1=V(:,loc);
W_commen_filter(sub,tt,:,:) = u1;
% 跨被试共同模板获取

% 计算皮尔逊相关性系数值
% template_sub_1 = squeeze(average_template_re(1,sub,tt,:,:))';
% template_sub = template_sub_1(1,1:Fs * tw);
% CCRS(sub,tt) = corrcoef(u1' * sub_test,template_sub);
template_1 = target_templete(10 * (tt-1)+1:10 * (tt-1)+10,:);
template = template_1(:,1:tw * Fs);
[~,~,r_CCRS] = canoncorr((u1' * sub_test)',template');
CCRS = r_CCRS;
end
