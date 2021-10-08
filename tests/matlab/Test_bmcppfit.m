%% test bmcppfit
% kai.herz@tuebingen.mpg.de
% 2021

%% where are we ?
script_path = [];
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    script_path = fileparts(matlab.desktop.editor.getActiveFilename);
else
    script_path = fileparts(which(mfilename));
end

%% add pulseq-cest to the path
% this is assuming that you chose the folder names according to the readme
pulseq_cest_path = fullfile(script_path, '..', '..', 'build',  '_deps', 'pulseq_cest-src');
addpath(genpath(pulseq_cest_path));
install_pulseqcest;

%% make results dir
results_path = fullfile(script_path, 'results');
mkdir(results_path);

%% seq file for simulation
seq_fn = 'APTw_3T_003_2uT_8block_DC95_834ms_braintumor.seq';
seq_fp = fullfile(script_path,seq_fn);

%% simulation settings file for pulseq-cest sim
yaml_fn = fullfile(script_path,'GM_3T_001_bmcppfit.yaml');

%% run simulation to get z-spec 
M_z = simulate_pulseqcest(seq_fp, yaml_fn);
figure;
hold on;
seq = mr.Sequence();
seq.read(seq_fp);
offsets = seq.definitions('offsets_ppm');
plot(offsets,M_z);
set ( gca, 'xdir', 'reverse');

%% add fit data to yaml fit file
y_init = yaml.ReadYaml(yaml_fn);
y_fit = y_init;
y_fit.fit_data = M_z;
y_fit.pulseq_file = seq_fn;

%% lets weight the regions of the cest pool a bit more
y_fit.weights = ones(size(M_z));
[~, id] = min(abs(offsets-y_fit.cest_pool.amide.dw));
y_fit.weights(id-3:id+3) = 2;
[~, id] = min(abs(offsets-y_fit.cest_pool.NOE_1.dw));
y_fit.weights(id-3:id+3) = 2;

%% openmp threads
y_fit.threads = uint8(4);

%% add fit parameter here.
% water t2
y_fit.fit_params.water_t2.start = y_fit.water_pool.t2 * 0.9;
y_fit.fit_params.water_t2.upper = 0.1;
y_fit.fit_params.water_t2.lower = 0.01;

% amide exchange rate
y_fit.fit_params.cest_1_k.start = y_fit.cest_pool.amide.k * 2;
y_fit.fit_params.cest_1_k.upper = 1000;
y_fit.fit_params.cest_1_k.lower = 5;

% noe concentration
y_fit.fit_params.cest_2_f.start = y_fit.cest_pool.NOE_1.f * 0.8;
y_fit.fit_params.cest_2_f.upper = y_fit.cest_pool.NOE_1.f * 5;
y_fit.fit_params.cest_2_f.lower = y_fit.cest_pool.NOE_1.f * 0.1;

%% simulate with start paramers
y_init.water_pool.t2 = y_fit.fit_params.water_t2.start;
y_init.cest_pool.amide.k = y_fit.fit_params.cest_1_k.start;
y_init.cest_pool.NOE_1.f = y_fit.fit_params.cest_2_f.start;
pre_fit_yaml_fn = fullfile(results_path,'GM_3T_001_bmcppfit_before_fit.yaml');
yaml.WriteYaml(pre_fit_yaml_fn,y_init);
M_z = simulate_pulseqcest(seq_fp, pre_fit_yaml_fn);
plot(offsets,M_z);

%% make yaml file for sim
yaml_fit_fn = fullfile(script_path,'yaml_fit.yaml');
yaml.WriteYaml(yaml_fit_fn, y_fit);

%% results file
yaml_res = fullfile(results_path,'fit_results.yaml');

%% run bmcppfit
bin_dir =  fullfile(script_path, '..', '..', 'install', 'bin');
bin_path = fullfile(bin_dir, 'bmcppfit');
system([bin_path ' -p=' yaml_fit_fn ' -o=' yaml_res]);

%% get results 
y_res = yaml.ReadYaml(yaml_res);
disp(['Simulated water t2 was ' num2str(y_fit.water_pool.t2) ', fitted k is ' num2str(y_res.water_t2)]);
disp(['Simulated amide k was ' num2str(y_fit.cest_pool.amide.k) ', fitted k is ' num2str(y_res.cest_1_k)]);
disp(['Simulated NOE f was ' num2str(y_fit.cest_pool.NOE_1.f) ', fitted f is ' num2str(y_res.cest_2_f)]);
y_init.water_pool.t2 = y_res.water_t2;
y_init.cest_pool.amide.k = y_res.cest_1_k;
y_init.cest_pool.NOE_1.f = y_res.cest_2_f;

%% run new sim for comparison
new_yaml_fn = fullfile(results_path,'GM_3T_001_bmcppfit_fitted.yaml');
yaml.WriteYaml(new_yaml_fn,y_init);
M_z = simulate_pulseqcest(seq_fp, new_yaml_fn);
plot(offsets,M_z);
%% add legend to plot
legend('Original Z-spectrum','Z-spectrum with Initial Parameters', 'Z-spectrum with Fitted Parameters');


