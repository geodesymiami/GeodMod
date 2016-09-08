function [dataset]=InitializeDataset(dataset,inverseopt);
%InitializeDataset   - initializes the dataset structure by adding the design matrix G_phaseramp  for faster inversions
%
%usage: [dataset]=InitializeDataset(dataset,inverseopt)
%
%Note: This can't be done in Qt2dataset because it depends on the inversion options
%
%FA, June 2007 
[d.d,d.coord,d.normalization,d.radarlook,d.datind,d.hgt,d.G_phaseramp,d.D_1,d.D_2,d.D_3,d.D_4,d.D_5,d.D_6,d.D_7,d.D_8]=datasetstructure2data(dataset,inverseopt);
dataset(1).fulldata=d; 
