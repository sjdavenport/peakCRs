run = 1

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

%% Example Run
Mag = 5;
Rad = [10,13,20]*(3/4)/10;
Smo = smo*[1,1,1];
Dim = [50,50];
centre_locs = {[20,40]'/2,[60,25]'/2,[75,75]'/2};

data_info = init_data_info2D( FWHM, nsubj, Rad, Smo, Dim, centre_locs, Mag );

nsubj = 50;
out = convCR(noise, FWHM)