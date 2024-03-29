%% Obtaining the coverage over the different settings for peaks of the mean in 2D
% Initialize the run for parallelization (this should be set at a different
% integer between 1 and 210 for each parallel run)
run = str2num(getenv('SGE_TASK_ID')); % Need 210 jobs

%% Initialize parameters for the model
modulo2 = mod(run,2) + 1;
modulo3 = mod(run,3) + 1;
modulo5 = mod(run,5) + 1;
modulo7 = mod(run,7) + 1;

nsubj_vec = 20:20:100;
nsubj_init = nsubj_vec(modulo5);

FWHM_vec = 3:9;
FWHM = FWHM_vec(modulo7);

smo_settings = [10, 15, 19];
smo = smo_settings(modulo3);

% Specify the directory to save the results
global PIloc
mainsaveloc = [PIloc, 'CoverageRates/2Dmean/ServerRunsReal/'];
mkdir(mainsaveloc)

%% Main run
for nsubj = [nsubj_init, 220-nsubj_init]
    
    addon = ['smo_', num2str(smo), '_FWHM_', num2str(FWHM), '_nsubj_', num2str(nsubj)];
    saveloc = [mainsaveloc, addon];
    
    Mag = 5;
    Rad = [10,13,20]*(3/4)/10;
    Smo = smo*[1,1,1];
    Dim = [50,50];
    centre_locs = {[20,40]'/2,[60,25]'/2,[75,75]'/2};
    
    data_info = init_data_info2D( FWHM, nsubj, Rad, Smo, Dim, centre_locs, Mag );
    
    % imagesc(data_info.meanonfinelat)
    niters = 5000;
    
    [rc.ellvals, rc.maxhist, rc.coverage, rc.bonf_coverage] = convCIcoverage(data_info, niters);
    
    save(saveloc, 'rc')
end