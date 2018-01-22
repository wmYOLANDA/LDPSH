close all;
clear; clc;

addpath('./utils/');
bitsNum = [16 32 64 128];
hash_methods = {'ITQ','ITQ-Ours'};

% construct data
load 'E:\feat_samples\sift.mat';
Xtrain = double(sift_learn');
Xtest = double(sift_query');
Xdata = double(sift_base');
gt = uint32(sift_groundtruth'); % each collum of Gt is a ground truth vector
clear sift_groundtruth sift_base sift_learn sift_query;

Xtrain = Xdata(1:size(Xtrain, 1) / 5,:);
meanVal = sum(Xtrain, 1) / size(Xtrain, 1);
Xtrain = bsxfun(@minus, Xtrain, meanVal);
Xtest = bsxfun(@minus, Xtest, meanVal);
Xdata = bsxfun(@minus, Xdata, meanVal);

% training and testing
nSel = size(Xdata, 1);
runtimes = 1;
iterNum = 5;
nNeighborNum = 10;
numepoch1 = 50;
numepoch2 = 10;
for i = 1:runtimes
    for j = 1:length(bitsNum)
        % data prepare
        [pcaMat, ~] = eigs(cov(Xtrain), bitsNum(j));
        Xtrainj = Xtrain * pcaMat;
        Xdataj = Xdata * pcaMat;
        Xtestj = Xtest * pcaMat;
        
        k = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Referred method %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        addpath('./ITQ/');
        
        [pcaMat, ~] = eigs(cov(Xtrainj), bitsNum(j)); % return the corresponding collum of the nbits max eigenvalue
        ITQparam.nbits = bitsNum(j);
        ITQparam.pcaW = pcaMat;
        ITQparam = trainITQ(Xtrainj, ITQparam);
        
        [Btrain_refer, Utrain_refer] = compressITQ(Xtrainj, ITQparam);
%         save('itq.mat', 'Btrain', 'Utrain', 'ITQparam');
        
        [Btrain, ~] = compressITQ(Xdataj, ITQparam);
        [Btest, ~] = compressITQ(Xtestj, ITQparam);
        
        [recall{i}{j, k}, precision{i}{j, k}, rec{i}{j, k}, pre{i}{j, k}, mAP1{i}{j, k}] = RecallAndPrecision_v2(Btrain', Btest', size(Btrain, 1), gt', size(gt, 2));
        
        k = 2;
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Our LDPSH method %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        addpath('./TSH/');
%         load('itq.mat');
        
        [Btrain, Btest] = trainLR_v1(Xtrainj, Btrain_refer, Utrain_refer, Xdataj, Xtestj, iterNum, nNeighborNum, numepoch1, numepoch2, bitsNum(j));
        [recall{i}{j, k}, precision{i}{j, k}, rec{i}{j, k}, pre{i}{j, k}, mAP1{i}{j, k}] = RecallAndPrecision_v2(Btrain', Btest', size(Btrain, 1), gt', size(gt, 2));
    end
end

mAP = mAP1;
for j = 1:2
    for i =1: length(bitsNum)
        tmp = zeros(size(mAP{1, 1}{i, j}));
        for k =1:runtimes
            tmp = tmp+mAP{1, k}{i, j};
        end
        MAP{i, j} = tmp/runtimes;
    end
    clear tmp;
end

% save result
% result_name = ['evaluations_sift1M_result' '.mat'];
% save(result_name, 'precision', 'recall', 'rec', 'pre', 'MAP', 'mAP', ...
%     'hash_methods', 'nhashmethods', 'bitsNum');

% plot attribution
line_width = 2;
marker_size = 8;
xy_font_size = 14;
legend_font_size = 12;
linewidth = 1.6;
title_font_size = xy_font_size;

%% show precision vs. recall , i is the selection of which bits.
figure('Color', [1 1 1]); hold on;
choose_times = 1;
choose_bits = 2;
for j = 1:2
    p = plot(rec{choose_times}{choose_bits, j}, pre{choose_times}{choose_bits, j});
    color = gen_color(j);
    marker = gen_marker(j);
    set(p,'Color', color)
    set(p,'Marker', marker);
    set(p,'LineWidth', line_width);
    if j < 8
        set(p, 'LineStyle', '--');
    else
        set(p, 'LineStyle', '-');
    end
    set(p,'MarkerSize', marker_size);
end

str_nbits =  num2str(bitsNum(choose_bits));
h1 = xlabel(['Recall @ ', str_nbits, ' bits']);
h2 = ylabel('Precision');
set(h1, 'FontSize', xy_font_size);
set(h2, 'FontSize', xy_font_size);
axis square;
hleg = legend(hash_methods);
set(hleg, 'FontSize', legend_font_size);
set(hleg,'Location', 'best');
set(gca, 'linewidth', linewidth);
box on; grid on; hold off;

%% show recall vs. the number of retrieved sample.
figure('Color', [1 1 1]); hold on;
pos = [1:100:1001 1200];
for j = 1: 2
    recc = recall{choose_times}{choose_bits, j};
    %p = plot(pos(1,1:posEnd), recc(1,1:posEnd));
    p = plot(pos, recc(pos, 1));
    color = gen_color(j);
    marker = gen_marker(j);
    set(p,'Color', color)
    set(p,'Marker', marker);
    if j < 8
        set(p, 'LineStyle', '--');
    else
        set(p, 'LineStyle', '-');
    end
    set(p,'LineWidth', line_width);
    set(p,'MarkerSize', marker_size);
end

str_nbits =  num2str(bitsNum(choose_bits));
set(gca, 'linewidth', linewidth);
h1 = xlabel('The number of retrieved samples');
h2 = ylabel(['Recall @ ', str_nbits, ' bits']);
set(h1, 'FontSize', xy_font_size);
set(h2, 'FontSize', xy_font_size);
axis square;
hleg = legend(hash_methods);
set(hleg, 'FontSize', legend_font_size);
set(hleg,'Location', 'best');
box on; grid on; hold off;

%% show precision vs. the number of retrieved sample.
figure('Color', [1 1 1]); hold on;
for j = 1: 2
    prec = precision{choose_times}{choose_bits, j};
    %p = plot(pos(1,1:posEnd), recc(1,1:posEnd));
    p = plot(pos, prec(pos, 1));
    color = gen_color(j);
    marker = gen_marker(j);
    set(p,'Color', color)
    set(p,'Marker', marker);
    if j < 8
        set(p, 'LineStyle', '--');
    else
        set(p, 'LineStyle', '-');
    end
    set(p,'LineWidth', line_width);
    set(p,'MarkerSize', marker_size);
end

str_nbits =  num2str(bitsNum(choose_bits));
set(gca, 'linewidth', linewidth);
h1 = xlabel('The number of retrieved samples');
h2 = ylabel(['Precision @ ', str_nbits, ' bits']);
set(h1, 'FontSize', xy_font_size);
set(h2, 'FontSize', xy_font_size);
axis square;
hleg = legend(hash_methods);
set(hleg, 'FontSize', legend_font_size);
set(hleg,'Location', 'best');
box on; grid on; hold off;

%% show mAP. This mAP function is provided by Yunchao Gong
figure('Color', [1 1 1]); hold on;
for j = 1: 2
    map = [];
    for i = 1: length(bitsNum)
        map = [map, MAP{i, j}];
    end
    p = plot(log2(bitsNum), map);
    color=gen_color(j);
    marker=gen_marker(j);
    set(p,'Color', color);
    set(p,'Marker', marker);
    if j <= 7
        set(p, 'LineStyle', '--');
    else
        set(p, 'LineStyle', '-');
    end
    set(p,'LineWidth', line_width);
    set(p,'MarkerSize', marker_size);
end

h1 = xlabel('Number of bits');
h2 = ylabel('mean Average Precision (mAP)');
set(h1, 'FontSize', xy_font_size);
set(h2, 'FontSize', xy_font_size);
axis square;
set(gca, 'xtick', log2(bitsNum));
set(gca, 'XtickLabel', {'8', '16', '32', '64', '128'});
set(gca, 'linewidth', linewidth);
hleg = legend(hash_methods);
set(hleg, 'FontSize', legend_font_size);
set(hleg, 'Location', 'best');
box on; grid on; hold off;