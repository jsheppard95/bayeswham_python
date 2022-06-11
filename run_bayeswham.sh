dim=1;

periodicity=[1];

T=298;

harmonicBiasesFile='./diala_phi_EXAMPLE/bias/harmonic_biases.txt';

histBinEdgesFile='./diala_phi_EXAMPLE/hist/hist_binEdges.txt';

histDir='./diala_phi_EXAMPLE/hist';

tol_WHAM=1E-15;

maxIter_WHAM=1E6;

steps_MH=1E6;

saveMod_MH=1E3;

printMod_MH=1E3;

maxDpStep_MH=5E-4;

seed_MH=200184;

prior='Dirichlet';

alpha=2;

python BayesWHAM.py $dim $periodicity $T $harmonicBiasesFile $histBinEdgesFile $histDir $tol_WHAM $maxIter_WHAM $steps_MH $saveMod_MH $printMod_MH $maxDpStep_MH $seed_MH $prior $alpha

