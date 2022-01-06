%% Initialize things (for the given setting)
run = str2num(getenv('SGE_TASK_ID')); % Need 210 jobs

modulo2 = mod(run,2);
modulo3 = mod(run,3) + 1;
modulo5 = mod(run,5) + 1;
modulo7 = mod(run,11) + 1;

if modulo2 == 0
    nsubj_vec = 20:20:100;
else
    nsubj_vec = 120:20:200;
end
nsubj = nsubj_vec(modulo5);

FWHM_vec = 10:20;
FWHM = FWHM_vec(modulo7);

smo_settings = [10, 15, 19];
smo = smo_settings(modulo3);
% smo = 5;

global PIloc
mainsaveloc = [PIloc, 'CoverageRates/2Dtstat/ServerRuns/'];
addon = ['smo_', num2str(smo), '_FWHM_', num2str(FWHM), '_nsubj_', num2str(nsubj)];
saveloc = [mainsaveloc, addon];
mkdir(mainsaveloc)

Mag = 5;
Rad = [10,13,20]*(3/4)/10;
Smo = smo*[1,1,1];
Dim = [50,50];
centre_locs = {[20,40]'/2,[60,25]'/2,[75,75]'/2};

data_info = init_data_info2D( FWHM, nsubj, Rad, Smo, Dim, centre_locs, Mag );

% imagesc(data_info.meanonfinelat)
niters = 5000;

[rc.ellvals, rc.maxhist, rc.coverage, rc.bonf_coverage] = convCIcoverage(data_info, niters, 't');

save(saveloc, 'rc')

% %% test
% FWHM = 7; smo = 19; Dim = [50,50]; Smo = smo*[1,1,1]; Mag = 5; nsubj = 10; 
% Rad = [10,13,20]*(3/4)/10;
% centre_locs = {[20,40]'/2,[60,25]'/2,[75,75]'/2};
% data_info = init_data_info2D( FWHM, nsubj, Rad, Smo, Dim, centre_locs, Mag );
