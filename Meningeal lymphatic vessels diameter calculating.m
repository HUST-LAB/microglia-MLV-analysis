clc,clear,close all

%% Parameters adjustment

pathname = 'N:\My article\PS-MLVs\20230811\Original data';  % The path of image.     
filename = '2.tif';  % The nameof image.
CropName = '2';

eps = 1e-3;                       % Minimum precision of the algorithm (used as a compensation when the denominator is zero).
rate = 0.8;                      % Threshold Segmentation Ratio (using the OTSU algorithm as a baseline).
R1 = 4;                           % Radius of image-closing computation circular formwork.
S = 2500;                          % Remove small connected domains of the image.
R2 = 1;                           % Radius of Morphological Corrosion circular Formwork.
Boundary = 0;                     % Boundary = 0 for preserving boundaries at image edges,Boundary = 1 for excluding boundaries at image edges.
% mode = 'up_down';                 % Lymphatic vessel wall extraction mode, 'up_down' for up-down direction, 'left_right' for left-right direction.
mode = 'left_right';
Slope = 5;                       % Median filter range of slope.
L = 30;                           % Radius of search for contralateral points.
gap = 1;                          % Approximate lymphatic vessel wall segment display interval, gap = 0 for not displaying, gap = 1 for displaying.
saveData = 1;                         % saveData = 1 for saving, saveData = 0 for not saving.

%% File reading

img = imread(fullfile(pathname, filename));
[m, n, k] = size(img);
if k > 3
    img(:, :, 4:k) = [];
    Igray = rgb2gray(img);
elseif k == 3
    Igray = rgb2gray(img);
elseif k ==1
    Igray = img;
end
%% Image pre-processing

T = graythresh(Igray);
Ibw = imbinarize(Igray, rate*T);

h1 = strel('disk', R1);   %%
Ibw_close = imclose(Ibw, h1);

Ibw_fill = imfill(Ibw_close, 'holes');

Ibw_fill = bwareaopen(Ibw_fill, S);    
% h2 = strel('line', 11, 135);

if Boundary == 1
    Ibw_edge = bwperim(Ibw_fill);
else
    h2 = strel('disk', R2);   
    Ibw_edge = Ibw_fill - imerode(Ibw_fill, h2);
end

%% Image display

figure, imshow(Igray)
figure, imshow(Ibw)
figure, imshow(Ibw_close)
figure, imshow(Ibw_fill)
figure, imshow(Ibw_edge)
figure, imshow(Ibw_edge)

%% Extraction of the lymphatic vessel wall, four types of extraction.

P_up = zeros(n, 2);
P_down = zeros(n, 2);
P_left = zeros(m, 2);
P_right = zeros(m, 2);

for j = 1:n

    for i = 1:m
        if Ibw_edge(i, j) == 1
            P_up(j, :) = [i, j];
            break;
        end
    end
end

for j = 1:n

    for i = m:-1:1
        if Ibw_edge(i, j) == 1
            P_down(j, :) = [i, j];
            break;
        end
    end
end


for i = 1:m
    for j = 1:n
        if Ibw_edge(i, j) == 1
            P_left(i, :) = [i, j];
            break;
        end
    end
end
            

for i = 1:m

    for j = n:-1:1
        if Ibw_edge(i, j) == 1
            P_right(i, :) = [i, j];
            break;
        end
    end
end

if strcmp(mode, 'up_down')
    P1 = P_up;
    P2 = P_down;
else
    P1 = P_left;
    P2 = P_right;
end
    
%% Calculation of the diameter

len = length(P1(:, 1));
D1 = zeros(len-2, 1);
D2 = zeros(len-2, 1);


K1 = P1(2:end, :) - P1(1:end-1, :);
K1 = K1(:, 1) ./ (K1(:, 2)+eps);
K1 = medfilt1(K1, Slope);

K2 = P2(2:end, :) - P2(1:end-1, :);
K2 = K2(:, 1) ./ (K2(:, 2)+eps);
K2 = medfilt1(K2, Slope);
%% First calculation

for i = 1:len-2
    x1 = P1(i+1, 1);
    y1 = P1(i+1, 2);
    if x1 == 0 || y1 == 0
        continue;
    end
    kt1 = K1(i+1);
    if strcmp(mode, 'up_down')
        Lmin = round(y1 - L  / sqrt(1 + kt1*kt1));
    else
        Lmin = round(x1 - L  / sqrt(1 + kt1*kt1));
    end
    if Lmin < 2
        Lmin = 2;
    end
    
    if strcmp(mode, 'up_down')
        Lmax = round(y1 + L  / sqrt(1 + kt1*kt1));
    else
        Lmax = round(x1 + L  / sqrt(1 + kt1*kt1));
    end
    if Lmax > len-1
        Lmax = len-1;
    end
    d = ones(Lmax - Lmin + 1, 1)*10000;
    f = ones(Lmax - Lmin + 1, 1)*10000;
    Line_d = zeros(Lmax - Lmin + 1, 2);
    for j = Lmin:1:Lmax
        
        
        
        x4 = P2(j, 1);
        y4 = P2(j, 2);
        
%         if ~isempty(find(P1(:, 1) == x4) & find(P1(:, 2) == y4))
%             continue;
%         end
        
        if x4 == 0 || y4 == 0
            continue;
        end
        Line_d(j - Lmin + 1, :) = [x4, y4];
        kn = (y4 - y1) / (x4 - x1 + eps);

        kt2 = K2(j);
        
        f(j - Lmin + 1) = (abs(kt1*kn + 1) + abs(kt2*kn + 1)) / 2;
        d(j - Lmin + 1) = sqrt((x4 - x1)^2 + (y4 - y1)^2);

    end
    
    if ~flag
        continue;
    end
    num = find(f == min(f));
    if ~(length(num) == 1)
        dtemp = d(num);
        num = find(d == min(dtemp));
        num = num(1);
    end
    D1(i) = d(num);
    
    if gap ~= 0
        if mod(i, gap) == 0
            line([y1, Line_d(num, 2)], [x1, Line_d(num, 1)], 'color', 'r');
        end
    end
    
end

%% Second calculation

for i = 1:len-2
    x1 = P2(i+1, 1);
    y1 = P2(i+1, 2);
    if x1 == 0 || y1 == 0
        continue;
    end
    kt1 = K2(i+1);
    if strcmp(mode, 'up_down')
        Lmin = round(y1 - L  / sqrt(1 + kt1*kt1));
    else
        Lmin = round(x1 - L  / sqrt(1 + kt1*kt1));
    end
    if Lmin < 2
        Lmin = 2;
    end
    
    if strcmp(mode, 'up_down')
        Lmax = round(y1 + L  / sqrt(1 + kt1*kt1));
    else
        Lmax = round(x1 + L  / sqrt(1 + kt1*kt1));
    end
    
    if Lmax > len-1
        Lmax = len-1;
    end
    d = ones(Lmax - Lmin + 1, 1)*10000;
    f = ones(Lmax - Lmin + 1, 1)*10000;
    Line_d = zeros(Lmax - Lmin + 1, 2);
    for j = Lmin:1:Lmax
        x4 = P1(j, 1);
        y4 = P1(j, 2);
        
%         if ~isempty(find(P2(:, 1) == x4) & find(P2(:, 2) == y4))
%             continue;
%         end
        
        
        if x4 == 0 || y4 == 0
            continue;
        end
        
        Line_d(j - Lmin + 1, :) = [x4, y4];
        kn = (y4 - y1 + eps) / (x4 - x1 + eps);
        
        kt2 = K1(j);
        
        f(j - Lmin + 1) = (abs(kt1*kn + 1) + abs(kt2*kn + 1)) / 2;
        d(j - Lmin + 1) = sqrt((x4 - x1)^2 + (y4 - y1)^2);
        
    end
    
    num = find(f == min(f));
    if ~(length(num) == 1)
        dtemp = d(num);
        num = find(d == min(dtemp));
        num = num(1);
    end
    D2(i) = d(num);
    if gap ~= 0
        if mod(i, gap) == 0
            line([y1, Line_d(num, 2)], [x1, Line_d(num, 1)], 'color', 'g');
        end
    end
    
end

x = 1;
%% Calculation of average diameter

D = (D1 + D2) / 2;
Ds = smooth(D, 30);

%% Display of the diameter change curves

% figure, plot(D)
position = [1 : length(Ds)];
% position = [1 : 1022] * 0.42;    % 0.42 indicates the pixel size()
figure, plot(position, Ds * 0.42', 'b-', 'linewidth', 1.5)
set(gca, 'FontName', 'Arial', 'FontSize' , 24);
% xlabel('position/\mum', 'FontName', 'Arial', 'FontSize' , 24);
% ylabel('diameter/\mum', 'FontName', 'Arial', 'FontSize' , 24);
% axis([position(1), position(2)])
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
,
if saveData == 1
 %% Data storage

ExcelName1 = [CropName, '_', 'D.xlsx'];
ExcelName2 = [CropName, '_', 'Ds.xlsx'];
MatName1 = [CropName, '_', 'D.mat'];
MatName2 = [CropName, '_', 'Ds.mat'];
txtname = [CropName, '_', 'config.txt'];
    
xlswrite(fullfile(pathname, ExcelName1), D);
xlswrite(fullfile(pathname, ExcelName2), Ds);
save(fullfile(pathname, MatName1), 'D');
save(fullfile(pathname, MatName2), 'Ds');

end