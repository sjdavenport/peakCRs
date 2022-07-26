contour_image_loc = '/vols/Scratch/ukbiobank/nichols/ContourInf/MNI/';

MNImask = imgload('MNImask');
cope_files = filesindir(contour_image_loc, '_cope5_MNI');

% Obtain the bounded mask
bounds = mask_bounds( MNImask );
bounded_mask = MNImask(bounds{:});

nsubj = 200;
imgs = loadsubs( 1:nsubj, '/vols/Scratch/ukbiobank/nichols/ContourInf/MNI/', 0, MNImask, 1, cope_files );

%%
global PIloc
for nsubj = 100:25:150
    imgs_2D = squeeze(imgs(:,31,:,1:nsubj));
    mask_2D = logical(squeeze(bounded_mask(:,31,:)));
    
    lat_data = Field(imgs_2D, mask_2D);
    lat_data_mean = mean(lat_data).field;
    
    FWHM = 2;
    
    smooth_data = convfield(lat_data, FWHM);
    smooth_mean = mean(smooth_data);
    smooth_std = std(smooth_data);
    smooth_CD = smooth_mean./smooth_std;
    smooth_CD_lat = smooth_CD.field;

    init_peak_locs = lmindices(smooth_CD.field, 2, mask_2D);
    
    save([PIloc, 'UKBanalysis/fMRI_application/smooth_CD_nsubj_', num2str(nsubj), '.mat'], 'smooth_CD_lat')

    meanfn = @(x) 0;
    
    out = convCR_t(lat_data, FWHM, meanfn, {init_peak_locs(:,1),init_peak_locs(:,2)} );
    save([PIloc, 'UKBanalysis/fMRI_application/CD_CIout_nsubj_', num2str(nsubj), '.mat'], 'out')

    if nsubj == 125
        smooth_mean = convfield(lat_data_mean, FWHM);
        init_peak_locs = lmindices(smooth_mean.field, 2, mask_2D);
        save([PIloc, 'UKBanalysis/fMRI_application/lat_mean_nsubj_', num2str(nsubj), '.mat'], 'lat_data_mean')

        out = convCR(lat_data, FWHM, meanfn, {init_peak_locs(:,1),init_peak_locs(:,2)} );
        save([PIloc, 'UKBanalysis/fMRI_application/CIout_nsubj_', num2str(nsubj), '.mat'], 'out')
    end
end

%%
MNImask = imgload('MNImask');
bounds = mask_bounds( MNImask );
bounded_mask = MNImask(bounds{:});
mask_2D = logical(squeeze(bounded_mask(:,31,:)));

color = zeros([size(mask_2D), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-mask_2D);

nsubj = 150;
load([PIloc, 'UKBanalysis/fMRI_application/CIout_nsubj_',num2str(nsubj),'.mat'], 'out')
index = 1;
specify_covmate = 0;

temp_ells = zeros(1, length(out.MFTD{index}));
for MFTD_iter = 1:length(out.MFTD{index})
    if specify_covmate
        covmateinv = nsubj*inv(cov(out.MFTD{index}));
       	temp_ells(MFTD_iter) = inellipse(zeros(2,1), covmateinv, out.MFTD{index}(MFTD_iter,:)');
    else
        temp_ells(MFTD_iter) = inellipse(zeros(2,1), nsubj*inv(out.cltSigmas{index}), out.MFTD{index}(MFTD_iter,:)');
    end
end
quant = 0.975
chi2quant_MFTD = prctile(temp_ells, 100*quant)
chi2quant_asym = chi2inv(quant, 2)

%%
cov(out.MFTD{index})
out.cltSigmas{index}/sqrt(nsubj)

%% Graphics
FWHM = 2;
load([PIloc, 'UKBanalysis/fMRI_application/lat_mean_nsubj_',num2str(nsubj),'.mat'], 'lat_data_mean')
load([PIloc, 'UKBanalysis/fMRI_application/CIout_nsubj_',num2str(nsubj),'.mat'], 'out')
resadd = 11;
params = ConvFieldParams([FWHM,FWHM], resadd, 0);
lat_data_mean_field = Field(lat_data_mean, logical(mask_2D));
smooth_mean = convfield(lat_data_mean_field, params);
imagesc(smooth_mean)

colored_slice = cell(1, 2);
slice = cell(1,2);
centre_point = cell(1,2);
quant = 0.975; % 0.95 for individual, 0.975 for joint (at the 0.05 level!)

useMFTDcov = 0;

useMFTD = 0;

for index = [1,2]
    temp_ells = zeros(1, length(out.MFTD{index}));
    if useMFTDcov
        covmateinv = inv(cov(out.MFTD{index}));
    else
        covmateinv = nsubj*inv(out.cltSigmas{index});
    end
    if useMFTD
        for MFTD_iter = 1:length(out.MFTD{index})
        %         temp_ells(MFTD_iter) = inellipse(zeros(2,1), nsubj*inv(out.cltSigmas{index}), out.MFTD{index}(MFTD_iter,:)');
        temp_ells(MFTD_iter) = inellipse(zeros(2,1), covmateinv, out.MFTD{index}(MFTD_iter,:)');
        end
        chi2quant = prctile(temp_ells, 100*quant);
    else
        chi2quant = chi2inv(quant, 2);
    end
%     Sigma = out.cltSigmas{index};
    slice{index} = zeros(smooth_mean.fieldsize);
    centre_point{index} = zeros(smooth_mean.fieldsize);
    for I = 1:length(smooth_mean.xvals{1})
        I
        for J = 1:length(smooth_mean.xvals{2})
            point = [smooth_mean.xvals{1}(I), smooth_mean.xvals{2}(J)]';
            [~, slice{index}(I, J)] = inellipse(point, covmateinv, out.max_locs(:,index), chi2quant);
            [~, centre_point{index}(I, J)] = inellipse(point, eye(2), out.max_locs(:,index), 0.01);
        end
    end
    colored_slice{index} = zeros([size(fliplr(slice{index})'), 3]);
    colored_slice{index}(:,:,1) = fliplr(slice{index})';
end

imagesc(fliplr(smooth_mean.field)')
hold on
for index = [1,2]
    im2 = imagesc(colored_slice{index});
    set(im2,'AlphaData',fliplr(slice{index})');
    hold on
end
color = zeros([size(fliplr(smooth_mean.mask)'), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-fliplr(smooth_mean.mask)');
axis off

global PIloc
save_filename = ['wholebrain_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(FWHM), '_resadd_', num2str(resadd)];
if useMFTDcov
    save_filename = [save_filename, '_MFTDcov'];
else
    save_filename = [save_filename, '_asymcov'];
end
if useMFTD
    save_filename = [save_filename, '_usedMTFD'];
else
    save_filename = [save_filename, '_usedasym'];
end
    
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

% First Peak
flipped_mean = fliplr(smooth_mean.field)';
surround = 5;
resadd = smooth_mean.resadd;
lm_locs = lmindices(flipped_mean, 2);
subset_indices = {lm_locs(1,1)-(surround+1)*resadd:lm_locs(1,1)+(surround-1)*resadd, ...
                                    lm_locs(2,1)-surround*resadd:lm_locs(2,1)+(surround)*resadd};
clf
imagesc(flipped_mean(subset_indices{:}))
subset_indices{3} = ':';
hold on
colored_slice_section = colored_slice{1}(subset_indices{:});
color_section = color(subset_indices{:});
slice_section = fliplr(slice{1})';
slice_section = slice_section(subset_indices{:});
im2 = imagesc(colored_slice_section);
set(im2,'AlphaData',slice_section);

point_section = fliplr(centre_point{1})';
point_section = point_section(subset_indices{:});
im3 = imagesc(color_section);
set(im3,'AlphaData',point_section);

xticks([])
yticks([])

save_filename = ['firstpeak_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(FWHM), '_resadd_', num2str(resadd)];
if useMFTDcov
    save_filename = [save_filename, '_MFTDcov'];
else
    save_filename = [save_filename, '_asymcov'];
end
if useMFTD
    save_filename = [save_filename, '_usedMTFD'];
else
    save_filename = [save_filename, '_usedasym'];
end
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

% Second Peak
flipped_mean = fliplr(smooth_mean.field)';
surround = 5;
resadd = smooth_mean.resadd;
lm_locs = lmindices(flipped_mean, 2);
subset_indices = {lm_locs(1,2)-(surround+1)*resadd:lm_locs(1,2)+(surround-1)*resadd, ...
                                    lm_locs(2,2)-surround*resadd:lm_locs(2,2)+(surround)*resadd};
clf
imagesc(flipped_mean(subset_indices{:}))
subset_indices{3} = ':';
hold on
colored_slice_section = colored_slice{2}(subset_indices{:});
slice_section = fliplr(slice{2})';
slice_section = slice_section(subset_indices{:});
im2 = imagesc(colored_slice_section);
set(im2,'AlphaData',slice_section);

point_section = fliplr(centre_point{2})';
point_section = point_section(subset_indices{:});
im3 = imagesc(color_section);
set(im3,'AlphaData',point_section);

xticks([])
yticks([])

save_filename = ['secondpeak_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(FWHM), '_resadd_', num2str(resadd)];
if useMFTDcov
    save_filename = [save_filename, '_MFTDcov'];
else
    save_filename = [save_filename, '_asymcov'];
end
if useMFTD
    save_filename = [save_filename, '_usedMTFD'];
else
    save_filename = [save_filename, '_usedasym'];
end
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

%%
flipped_mean = fliplr(lat_data_mean)';
surround = 5;
imagesc(flipped_mean(52-surround+1:52+surround-1, 17-surround:17+surround-1))
%%
surround = 3
imagesc(flipped_mean(52-surround:52+surround, 57-surround:57+surround-1))



%% Save tstat image
tstat = mvtstat(imgs);

tstat2save = zeros(91,109,91);
tstat2save(bounds{:}) = tstat;
global PIloc
imgsave( tstat2save, 'tstat_100', [PIloc, 'UKBanalysis/fMRI_application/'] )

%%
global PIloc
tstat = imgload([PIloc, 'UKBanalysis/fMRI_application/tstat_100']);

viewbrain(tstat, [0,45,0])

%%
lmindices(lat_data_mean, 2)

%%
MNImask = imgload('MNImask');
bounds = mask_bounds( MNImask );
bounded_mask = MNImask(bounds{:});
bounded_tstat = tstat(bounds{:});
viewbrain(bounded_tstat, [0,31,0], bounded_mask, 0)

%%
imagesc(squeeze(bounded_tstat(:,31,:)))

%%

bounded_tstat_extended = pad_vals(bounded_tstat, 4 );
bounded_mask_extended = pad_vals(bounded_mask, 4 );
% tstat2plot = zeros(size(bounded_tstat) + 4);
% bounded_mask_extended = zeros(size(bounded_tstat) + 4);
% tstat2plot(3:end-2,3:end-2,3:end-2) = bounded_tstat;
% bounded_mask_extended =
viewbrain(bounded_tstat_extended, [0,35,0], bounded_mask_extended)

%%
viewbrain(tstat, [0,35,0], bounded_mask)

% %%
% index = 1;
% temp_ells = zeros(1, length(MFdist));
% % MFdist_bounded = MFdist(abs(MFdist(:,1)) < 100, 2);
% for MFTD_iter = 1:length(MFdist)
%     temp_ells(MFTD_iter) = inellipse(zeros(2,1), nsubj*inv(out.cltSigmas{index}), MFdist(MFTD_iter,:)');
% end
% quant = 0.975;
% chi2quant = prctile(temp_ells, 100*quant)
% chi2quant = chi2inv(quant, 2)

%% Analyse MFTD distbn
index = 1; D = 2;
Delta = zeros(D, D*(D+1)/2);
covmate = [out.Lambda{index}, Delta; Delta', out.Omega{index}];
[MFdist, ~, DH, DH_nobounds] = MFTD( reshape(out.peakderiv2{index}, [D,D]), covmate, 100 );

% h1 = histogram(MFdist(abs(MFdist(:,1)) < 10, 1));
% h1 = histogram(MFdist(:,1));
histogram(DH_nobounds)

%% Generate from the empirical distbn
niters = 100000;
nulldist = mvnrnd(zeros(1,D), out.Lambda{index}, niters)/sqrt(nsubj);
denom = inv(reshape(out.peakderiv2{index}, [D,D]));
for I = 1:length(nulldist)
    nulldist(I,:) = (denom*( nulldist(I,:)'))';
end

%%
h1 = histogram(MFdist(abs(MFdist(:,1)) < 10, 1));
hold on
h2 = histogram(nulldist(:,1));
h2.BinWidth = h1.BinWidth

%%
histogram(nulldist(:,1))