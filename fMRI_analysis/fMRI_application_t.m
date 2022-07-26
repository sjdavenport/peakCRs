%%
MNImask = imgload('MNImask');
bounds = mask_bounds( MNImask );
bounded_mask = MNImask(bounds{:});
mask_2D = logical(squeeze(bounded_mask(:,31,:)));

color = zeros([size(mask_2D), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-mask_2D);

nsubj = 125;
load([PIloc, 'UKBanalysis/fMRI_application/CD_CIout_nsubj_',num2str(nsubj),'.mat'], 'out')
index = 1;

%%
cov(out.MFTD{index})
out.cltSigmas{index}/sqrt(nsubj)

%% Graphics
FWHM = 1.3;
load([PIloc, 'UKBanalysis/fMRI_application/lat_mean_nsubj_',num2str(nsubj),'.mat'], 'lat_data_mean')
load([PIloc, 'UKBanalysis/fMRI_application/smooth_CD_nsubj_',num2str(nsubj),'.mat'], 'smooth_CD_lat')

resadd = 11;
params = ConvFieldParams([FWHM,FWHM], resadd, 0);
smooth_CD_field = Field(smooth_CD_lat, logical(mask_2D));
smooth_CD = convfield(smooth_CD_field, params);

imagesc(smooth_CD)
figure
imagesc(smooth_CD_lat)
%%
colored_slice = cell(1, 2);
slice = cell(1,2);
centre_point = cell(1,2);
quant = 0.975; % 0.95 for individual, 0.975 for joint (at the 0.05 level!)

useMFTDcov = 0;

useMFTD = 0;

for index = [1,2]
    covmateinv = nsubj*inv(out.cltSigmas{index});
    chi2quant = chi2inv(quant, 2);

    slice{index} = zeros(smooth_CD.fieldsize);
    centre_point{index} = zeros(smooth_CD.fieldsize);
    for I = 1:length(smooth_CD.xvals{1})
        I
        for J = 1:length(smooth_CD.xvals{2})
            point = [smooth_CD.xvals{1}(I), smooth_CD.xvals{2}(J)]';
            [~, slice{index}(I, J)] = inellipse(point, covmateinv, out.max_locs(:,index), chi2quant);
            [~, centre_point{index}(I, J)] = inellipse(point, eye(2), out.max_locs(:,index), 0.01);
        end
    end
    colored_slice{index} = zeros([size(fliplr(slice{index})'), 3]);
    colored_slice{index}(:,:,1) = fliplr(slice{index})';
end

imagesc(fliplr(smooth_CD.field)')
hold on
for index = [1,2]
    im2 = imagesc(colored_slice{index});
    set(im2,'AlphaData',fliplr(slice{index})');
    hold on
end
color = zeros([size(fliplr(smooth_CD.mask)'), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-fliplr(smooth_CD.mask)');
axis off

global PIloc
save_filename = ['CD_wholebrain_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(10*FWHM), '_resadd_', num2str(resadd)];
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

%% First Peak
flipped_mean = fliplr(smooth_CD.field)';
surround = 5;
resadd = smooth_CD.resadd;
lm_locs = lmindices(nan2zero(flipped_mean), 2);
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

save_filename = ['CD_firstpeak_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(10*FWHM)];
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

%% Second Peak
flipped_mean = fliplr(smooth_CD.field)';
surround = 5;
resadd = smooth_CD.resadd;
lm_locs = lmindices(nan2zero(flipped_mean), 2);
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

save_filename = ['CD_secondpeak_nsubj_', num2str(nsubj), '_extraFWHM_', num2str(10*FWHM)];
export_fig([PIloc, 'Figures/fMRI/', save_filename, '.tif'], '-transparent')

%% Plotting Cohen's d
load([PIloc, 'UKBanalysis/fMRI_application/smooth_CD_nsubj_125.mat'], 'smooth_CD_lat')

