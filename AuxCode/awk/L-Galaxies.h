struct LGalaxy {
   long long GalID;
   long long HaloID;
   long long FirstProgGal;
   long long NextProgGal;
   long long LastProgGal;
   long long FOFCentralGal;
   long long FileTreeNr;
   long long DescendantGal;
   long long MainLeafId;
   long long TreeRootId;
   long long SubID;
   long long MMSubID;
   int PeanoKey;
   float Redshift;
   int Type;
   int SnapNum;
   float LookBackTimeToSnap;
   float CentralMvir;
   float CentralRvir;
   float DistanceToCentralGal[3];
   float Pos[3];
   float Vel[3];
   int Len;
   float Mvir;
   float Rvir;
   float Vvir;
   float Vmax;
   float GasSpin[3];
   float StellarSpin[3];
   float InfallVmax;
   float InfallVmaxPeak;
   int InfallSnap;
   float InfallHotGas;
   float HotRadius;
   float OriMergTime;
   float MergTime;
   float ColdGas;
   float StellarMass;
   float BulgeMass;
   float DiskMass;
   float HotGas;
   float EjectedMass;
   float BlackHoleMass;
   float ICM;
   float MetalsColdGas;
   float MetalsStellarMass;
   float MetalsBulgeMass;
   float MetalsDiskMass;
   float MetalsHotGas;
   float MetalsEjectedMass;
   float MetalsICM;
   float PrimordialAccretionRate;
   float CoolingRadius;
   float CoolingRate;
   float CoolingRate_beforeAGN;
   float QuasarAccretionRate;
   float RadioAccretionRate;
   float Sfr;
   float SfrBulge;
   float XrayLum;
   float BulgeSize;
   float StellarDiskRadius;
   float GasDiskRadius;
   float CosInclination;
   int DisruptOn;
   int MergeOn;
   float ObsMagDust[20];
   float ObsMag[20];
   float ObsMagBulge[20];
   float MassWeightAge;
   float rbandWeightAge;
   int sfh_ibin;
   int sfh_numbins;
   float sfh_DiskMass[20];
   float sfh_BulgeMass[20];
   float sfh_ICM[20];
   float sfh_MetalsDiskMass[20];
   float sfh_MetalsBulgeMass[20];
   float sfh_MetalsICM[20];
};
