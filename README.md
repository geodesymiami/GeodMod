## Geodetic modelling software in Matlab


### 1. Installation

#### 1.1 Setup environment path

Add the following to your _~/.bashrc_:

```
alias matlab=/Applications/MATLAB_R2019a.app/bin/matlab
export MATLAB_USE_USERPATH=1

##--------- GEODMOD ------------------##
export GEODMOD_HOME=~/Documents/test/test1/geodmod
export GEODMOD_TESTDATA=~/Documents/testdata/geodmod
export GEODMOD_TESTBENCH=~/Downloads/TESTBENCH
```

#### 1.2 Download

```
git clone https://github.com/geodesymiami/GeodMod.git $GEODMOD_HOME/geodmod
```

#### 1.3 Setup matlab path

Add the following to your _startup.m_ file in _USER/matlab_ for Matlab to initiate the configuration:

```matlab
%% Setting for GeodMod
disp('Setting paths for geodmod...')
run( [ getenv('GEODMOD_HOME') filesep 'addpath_geodmod'] )
```

### 2. Test Run
Download test data

```
mkdir -p $GEODMOD_TESTDATA
cd $GEODMOD_TESTDATA
wget https://zenodo.org/record/3348961/files/geodmod_testdata.xz
tar -xvJf  geodmod_testdata.xz

Start matlab program from terminal to initiate the environment variables defined above.

```
matlab
>>geodmod Darwin.min
```
