;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_zerofloats, LGs 
; test whether floats are NaN or too small for SQLServer
; if so, set offending values to 0
; assumes the existence of a function testFloat accepting an array of floats
 sel = testFloat(LGs.Redshift)
 if(sel(0) gt -1) then begin
     LGs[sel].Redshift = 0
 endif
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     LGs[sel].LookBackTimeToSnap = 0
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralMvir = 0
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralRvir = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(0) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(1) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(2) = 0
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(0) = 0
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(1) = 0
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(2) = 0
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(0) = 0
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(1) = 0
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(2) = 0
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Mvir = 0
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Rvir = 0
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Vvir = 0
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     LGs[sel].Vmax = 0
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(0) = 0
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(1) = 0
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(2) = 0
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(0) = 0
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(1) = 0
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(2) = 0
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmax = 0
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmaxPeak = 0
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallHotGas = 0
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].HotRadius = 0
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].OriMergTime = 0
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].MergTime = 0
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGas = 0
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarMass = 0
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeMass = 0
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskMass = 0
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].HotGas = 0
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].EjectedMass = 0
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BlackHoleMass = 0
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     LGs[sel].ICM = 0
 endif
 sel = testFloat(LGs.MetalsColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsColdGas = 0
 endif
 sel = testFloat(LGs.MetalsStellarMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsStellarMass = 0
 endif
 sel = testFloat(LGs.MetalsBulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsBulgeMass = 0
 endif
 sel = testFloat(LGs.MetalsDiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsDiskMass = 0
 endif
 sel = testFloat(LGs.MetalsHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsHotGas = 0
 endif
 sel = testFloat(LGs.MetalsEjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsEjectedMass = 0
 endif
 sel = testFloat(LGs.MetalsICM)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsICM = 0
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].PrimordialAccretionRate = 0
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRadius = 0
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate = 0
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate_beforeAGN = 0
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].QuasarAccretionRate = 0
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].RadioAccretionRate = 0
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     LGs[sel].Sfr = 0
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     LGs[sel].SfrBulge = 0
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     LGs[sel].XrayLum = 0
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeSize = 0
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarDiskRadius = 0
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].GasDiskRadius = 0
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     LGs[sel].CosInclination = 0
 endif
 sel = testFloat(LGs.ObsMagDust(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(0) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(1) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(2) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(3) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(4) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(5))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(5) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(6))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(6) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(7))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(7) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(8))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(8) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(9))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(9) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(10))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(10) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(11))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(11) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(12))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(12) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(13))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(13) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(14))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(14) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(15))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(15) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(16))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(16) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(17))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(17) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(18))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(18) = 0
 endif
 sel = testFloat(LGs.ObsMagDust(19))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagDust(19) = 0
 endif
 sel = testFloat(LGs.ObsMag(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(0) = 0
 endif
 sel = testFloat(LGs.ObsMag(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(1) = 0
 endif
 sel = testFloat(LGs.ObsMag(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(2) = 0
 endif
 sel = testFloat(LGs.ObsMag(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(3) = 0
 endif
 sel = testFloat(LGs.ObsMag(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(4) = 0
 endif
 sel = testFloat(LGs.ObsMag(5))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(5) = 0
 endif
 sel = testFloat(LGs.ObsMag(6))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(6) = 0
 endif
 sel = testFloat(LGs.ObsMag(7))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(7) = 0
 endif
 sel = testFloat(LGs.ObsMag(8))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(8) = 0
 endif
 sel = testFloat(LGs.ObsMag(9))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(9) = 0
 endif
 sel = testFloat(LGs.ObsMag(10))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(10) = 0
 endif
 sel = testFloat(LGs.ObsMag(11))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(11) = 0
 endif
 sel = testFloat(LGs.ObsMag(12))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(12) = 0
 endif
 sel = testFloat(LGs.ObsMag(13))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(13) = 0
 endif
 sel = testFloat(LGs.ObsMag(14))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(14) = 0
 endif
 sel = testFloat(LGs.ObsMag(15))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(15) = 0
 endif
 sel = testFloat(LGs.ObsMag(16))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(16) = 0
 endif
 sel = testFloat(LGs.ObsMag(17))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(17) = 0
 endif
 sel = testFloat(LGs.ObsMag(18))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(18) = 0
 endif
 sel = testFloat(LGs.ObsMag(19))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMag(19) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(0) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(1) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(2) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(3))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(3) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(4))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(4) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(5))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(5) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(6))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(6) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(7))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(7) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(8))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(8) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(9))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(9) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(10))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(10) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(11))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(11) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(12))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(12) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(13))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(13) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(14))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(14) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(15))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(15) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(16))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(16) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(17))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(17) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(18))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(18) = 0
 endif
 sel = testFloat(LGs.ObsMagBulge(19))
 if(sel(0) gt -1) then begin
     LGs[sel].ObsMagBulge(19) = 0
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].MassWeightAge = 0
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].rbandWeightAge = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(3) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(4) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(5) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(6) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(7) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(8) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(9) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(10) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(11) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(12) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(13) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(14) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(15) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(16) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(17) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(18) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(19) = 0
 endif
end
