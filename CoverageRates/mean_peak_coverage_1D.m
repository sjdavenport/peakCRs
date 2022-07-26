run = str2num(getenv('SGE_TASK_ID')); % Need 3*7*2 = 42 different jobs!

%%
peakparams = {[2,2], [2,1.5], [1.5, 3]};

if mod(run,2) == 0
    sim_type = 'N';
else
    sim_type = 'tnoise';
end

modulo3 = mod(run,3) + 1;
params2use = peakparams(modulo3);
npeaks = 3;

global PIloc
mainsaveloc = [PIloc, 'CoverageRates/1Dmeanmultpeaks/ServerRuns/'];
mkdir(mainsaveloc)
modulo3 = num2str(modulo3);
addon = ['server_npeaks_', num2str(npeaks),'_paramsetting_', modulo3];

FWHM_vec = 3:9;
rc.FWHM_vec = FWHM_vec(mod(run, 7) + 1);
addon = [addon, '_FWHM_', num2str(rc.FWHM_vec)];
addon = [addon, '_sim_type_', sim_type];
saveloc = [mainsaveloc, addon];

peak_width = 7;
peak_half_width = peak_width/2;
peakspec = cell(1,npeaks);
peak_centres = 5:10:(10*(npeaks-1)+5);
for peak = 1:npeaks
    peakspec{peak} = [peak_centres(peak) - peak_half_width, peak_centres(peak) + peak_half_width];
end

increm = 0.001;
smo = 0.01; %I.e. the peaks are not smoothed before doing stuff! 
npeaks = length(peakspec);
xvals = 1:increm:(10*npeaks); 
[sigstore, meanfn] = peakgen1D( xvals, peakspec, params2use, 1, smo);

%% Saves the coverage for a parabolic signal for different FWHM and numbers of subjects
rc.nsubj_vec = 20:20:200;
rc.coverage_norm = zeros(length(rc.FWHM_vec), length(rc.nsubj_vec));
rc.niters = 5000;
rc.ellvals = cell(length(rc.FWHM_vec), length(rc.nsubj_vec));
rc.coverage = cell(length(rc.FWHM_vec), length(rc.nsubj_vec));

data_info.xvals_vecs = {xvals}; data_info.resadd = 1/(data_info.xvals_vecs{1}(2) - data_info.xvals_vecs{1}(1)) - 1;
data_info.resadd = round(data_info.resadd); %As CFparams doesn't like non integers
data_info.lat_signal = sigstore;

data_info.npeaks = npeaks;
data_info.maxlocs = xvals(lmindices(sigstore, npeaks));

data_info.lat_std = 1;
data_info.dim = max(data_info.xvals_vecs{1});
data_info.D = 1;
data_info.meanfn = meanfn;

rc.data_info = data_info;

for I = 1:length(rc.FWHM_vec)
    data_info.FWHM = rc.FWHM_vec(I);
    rc.FWHM_vec(I)
    
    for J = 1:length(rc.nsubj_vec)
        data_info.nsubj = rc.nsubj_vec(J);
        rc.nsubj_vec(J)
        
        [rc.ellvals{I,J}, rc.maxhist{I,J}, rc.coverage{I,J}, rc.bonf_coverage{I,J}] = convCIcoverage(data_info, rc.niters, 'm', sim_type);
    
        save(saveloc, 'rc')
    end
end

