## Geodetic modelling software in Matlab


### 1. Installation

#### 1.1 Download

```
git clone https://github.com/falkamelung/GeodMod.git
```

#### 1.2 Setup environment path

Add the following to your _~/.cshrc_:

```
setenv GEODMOD_HOME       ~/development/matlab/GeodMod
setenv GEODMOD_TESTDATA   ~/insarlab/test/geodmod
setenv GEODMOD_TESTBENCH  ~/insarlab/test/geodmod_testbench
```

#### 1.3 Setup matlab path

Add the following to your _startup.m_ file for Matlab to initiate the configuration:

```matlab
%% Setting for GeodMod
disp('Setting paths for geodmod...')
run( [ getenv('GEODMOD_HOME') filesep 'addpath_geodmod'] )
```

### 2. Test Run

```matlab
geodmod Darwin.min
```
