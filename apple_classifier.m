clc;
clear;
close all;

%% Display info
disp('Apple Classifier Program is started');

%% Loading the Dataset

cd sibhaye_bimar
a1=dir;
for i=3:length(a1)
    b1=a1(i).name;
    cd(b1);
    a2=dir;
    for j1=3:length(a2)
        b2=a2(j1).name;
        cd(b2);
        a3=dir;
        for k=3:length(a3)
            b3=a3(k).name;
            sibhaye_bimar{i-2,j1-2,k-2}=imresize(imread(b3),0.01);
        end
        cd ..
    end
    cd ..
end
cd ..

cd sibhaye_salem
a4=dir;
for i=3:length(a4)
    b4=a4(i).name;
    cd(b4);
    a5=dir;
    for j2=3:length(a5)
        b5=a5(j2).name;
        cd(b5);
        a6=dir;
        for k=3:length(a6)
            b6=a6(k).name;
            if k==3 && j2==3 && i==3
                X=imread(b6);
            end
            sibhaye_salem{i-2,j2-2,k-2}=imresize(imread(b6),0.01);
        end
        cd ..
    end
    cd ..
end
cd ..

gal_va_kerme_sibe_zard1=sibhaye_bimar(1,:,:);
gal_va_senak_sibe_ghermez1=sibhaye_bimar(2,:,:);
gal_va_senake_sibe_zard1=sibhaye_bimar(3,:,:);
gal_va_zangare_sibe_ghermez1=sibhaye_bimar(4,:,:);
gal_va_zangare_sibe_zard1=sibhaye_bimar(5,:,:);
gale_ghermez1=sibhaye_bimar(6,:,:);
gale_zard1=sibhaye_bimar(7,:,:);
kerm_va_senake_sibe_ghermez1=sibhaye_bimar(8,:,:);
kerm_va_senake_sibe_zard1=sibhaye_bimar(9,:,:);
kerm_va_zangare_sibe_zard1=sibhaye_bimar(10,:,:);
kerme_ghermez1=sibhaye_bimar(11,:,:);
kerme_zard1=sibhaye_bimar(12,:,:);
senak_va_zangare_sibe_zard1=sibhaye_bimar(13,:,:);
senake_ghermez1=sibhaye_bimar(14,:,:);
senake_zard1=sibhaye_bimar(15,:,:);
zangare_zard1=sibhaye_bimar(16,:,:);


for i1=1:j1-1
    if isempty(gal_va_kerme_sibe_zard1{i1})==0
        gal_va_kerme_sibe_zard{i1}=gal_va_kerme_sibe_zard1{i1};
    end
    if isempty(gal_va_senak_sibe_ghermez1{i1})==0
        gal_va_senak_sibe_ghermez{i1}=gal_va_senak_sibe_ghermez1{i1};
    end
    if isempty(gal_va_senake_sibe_zard1{i1})==0
        gal_va_senake_sibe_zard{i1}=gal_va_senake_sibe_zard1{i1};
    end
    if isempty(gal_va_zangare_sibe_ghermez1{i1})==0
        gal_va_zangare_sibe_ghermez{i1}=gal_va_zangare_sibe_ghermez1{i1};
    end
    if isempty(gal_va_zangare_sibe_zard1{i1})==0
        gal_va_zangare_sibe_zard{i1}=gal_va_zangare_sibe_zard1{i1};
    end
    if isempty(gale_ghermez1{i1})==0
        gale_ghermez{i1}=gale_ghermez1{i1};
    end
    if isempty(gale_zard1{i1})==0
        gale_zard{i1}=gale_zard1{i1};
    end
    if isempty(kerm_va_senake_sibe_ghermez1{i1})==0
        kerm_va_senake_sibe_ghermez{i1}=kerm_va_senake_sibe_ghermez1{i1};
    end
    if isempty(kerm_va_senake_sibe_zard1{i1})==0
        kerm_va_senake_sibe_zard{i1}=kerm_va_senake_sibe_zard1{i1};
    end
    if isempty(kerm_va_zangare_sibe_zard1{i1})==0
        kerm_va_zangare_sibe_zard{i1}=kerm_va_zangare_sibe_zard1{i1};
    end
    if isempty(kerme_ghermez1{i1})==0
        kerme_ghermez{i1}=kerme_ghermez1{i1};
    end
    if isempty(kerme_zard1{i1})==0
        kerme_zard{i1}=kerme_zard1{i1};
    end
    if isempty(senak_va_zangare_sibe_zard1{i1})==0
        senak_va_zangare_sibe_zard{i1}=senak_va_zangare_sibe_zard1{i1};
    end
    if isempty(senake_ghermez1{i1})==0
        senake_ghermez{i1}=senake_ghermez1{i1};
    end
    if isempty(senake_zard1{i1})==0
        senake_zard{i1}=senake_zard1{i1};
    end
    if isempty(zangare_zard1{i1})==0
        zangare_zard{i1}=zangare_zard1{i1};
    end
end

sibeghermez1=sibhaye_salem(1,:,:);
sibezard1=sibhaye_salem(2,:,:);
for i2=1:j2-1
    if isempty(sibezard1{i2})==0
        sibezard{i2}=sibezard1{i2};
    end
    if isempty(sibeghermez1{i2})==0
        sibeghermez{i2}=sibeghermez1{i2};
    end
end

%% Creating an Sparse image to form the initial dictionary and weights



A = rand(64);
num_trials=10000;
batch_size=100;

BUFF=4;
[L M]=size(A);
sz=sqrt(L);

eta = 1.0;
noise_var= 0.01;
beta= 2.2;
sigma=0.316;
tol=.01;

VAR_GOAL=0.1;
S_var=VAR_GOAL*ones(M,1);
var_eta=.001;
alpha=.02;
gain=sqrt(sum(A.*A))';

%% Loading a sample image to form the initial dictionary and weights

D=X;

if (exist('disp_handle','var'))
  update_network(A,S_var,disp_handle);
else
  disp_handle=display_network(A,S_var);
end


D=D(:,:,1);
D=im2double(D);

%% Display info

N = size(X,1);
tekst = char(' ', ['Apple Classifier ']);
if exist('L', 'var') && (L < size(X,2))
    I = randperm(size(X,2));
    X = im2double(X(:,I(1:L)));  
    tekst = char(tekst, ['Use ',int2str(L),' of ',int2str(size(X,2)), ...
        ' randomly selected vectors from X.']);
end
tekst = char(tekst, 'w1 and w2 are mean values for 1 or 2-norms ');
tekst = char(tekst, 'of X (image), R (error) or W (Dictionary).');
tekst = char(tekst, 'Summary: ');
K = size(D,2);
L = size(X,2);

format1 = ' %30s :   time     w1     w2   accuracy';
format2 = ' %30s :  %5.2f   %5.2f  %5.2f   %5.2f';
 
tekst = char(tekst, sprintf(format1, 'Method'));
    
%% Reloading sparseapprox function to optimize considered dictionary and weights


[W, res] = sparseapprox(X, D, 'Graph Regulized', 'tnz',4, 'v',2, 'doOP',0);
t = sprintf(format2, 'Graph Regulized', ...
    res.time, res.snr, mean(res.norm0W), mean(res.norm1W) );
tekst = char(tekst, t);

disp(tekst);

%% Test of an image for Classification

fprintf('\n\n\n');
disp('Please enter the path and name of the given image, between two quotes');
disp('Example: C:\Users\...\Desktop\Program\sibhaye_salem\sibe_ghermez\1\DSC00729.png');
X_test=input('Enter here:');
if exist('X')
    X_test=imread(X_test);
    W= sparseapprox(X_test, D, 'Test', 'tnz',4, 'v',2, 'doOP',0);
end
fprintf('The input image belongs to %s\n',W);