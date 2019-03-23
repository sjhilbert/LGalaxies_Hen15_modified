;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LGalaxy_testfloats, LGs, nstart 
; test whether floats are NaN or too small for SQLServer
; assumes the existence of a function testFloat accepting an array of floats
 badranges = []
 bad = 0
 sel = testFloat(LGs.Redshift)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Redshift --- ', nstart+sel
     print, 'Redshift --- ', LGs[sel].Redshift
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'LookBackTimeToSnap --- ', nstart+sel
     print, 'LookBackTimeToSnap --- ', LGs[sel].LookBackTimeToSnap
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralMvir --- ', nstart+sel
     print, 'CentralMvir --- ', LGs[sel].CentralMvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralRvir --- ', nstart+sel
     print, 'CentralRvir --- ', LGs[sel].CentralRvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[0] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[1] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[2] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[0] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[1] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[2] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[0] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[1] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[2] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mvir --- ', nstart+sel
     print, 'Mvir --- ', LGs[sel].Mvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Rvir --- ', nstart+sel
     print, 'Rvir --- ', LGs[sel].Rvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vvir --- ', nstart+sel
     print, 'Vvir --- ', LGs[sel].Vvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vmax --- ', nstart+sel
     print, 'Vmax --- ', LGs[sel].Vmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[0] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[1] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[2] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[0] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[1] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[2] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmax --- ', nstart+sel
     print, 'InfallVmax --- ', LGs[sel].InfallVmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmaxPeak --- ', nstart+sel
     print, 'InfallVmaxPeak --- ', LGs[sel].InfallVmaxPeak
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallHotGas --- ', nstart+sel
     print, 'InfallHotGas --- ', LGs[sel].InfallHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotRadius --- ', nstart+sel
     print, 'HotRadius --- ', LGs[sel].HotRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'OriMergTime --- ', nstart+sel
     print, 'OriMergTime --- ', LGs[sel].OriMergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MergTime --- ', nstart+sel
     print, 'MergTime --- ', LGs[sel].MergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGas --- ', nstart+sel
     print, 'ColdGas --- ', LGs[sel].ColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarMass --- ', nstart+sel
     print, 'StellarMass --- ', LGs[sel].StellarMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeMass --- ', nstart+sel
     print, 'BulgeMass --- ', LGs[sel].BulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskMass --- ', nstart+sel
     print, 'DiskMass --- ', LGs[sel].DiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotGas --- ', nstart+sel
     print, 'HotGas --- ', LGs[sel].HotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'EjectedMass --- ', nstart+sel
     print, 'EjectedMass --- ', LGs[sel].EjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BlackHoleMass --- ', nstart+sel
     print, 'BlackHoleMass --- ', LGs[sel].BlackHoleMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ICM --- ', nstart+sel
     print, 'ICM --- ', LGs[sel].ICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsColdGas --- ', nstart+sel
     print, 'MetalsColdGas --- ', LGs[sel].MetalsColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsStellarMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsStellarMass --- ', nstart+sel
     print, 'MetalsStellarMass --- ', LGs[sel].MetalsStellarMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsBulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsBulgeMass --- ', nstart+sel
     print, 'MetalsBulgeMass --- ', LGs[sel].MetalsBulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsDiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsDiskMass --- ', nstart+sel
     print, 'MetalsDiskMass --- ', LGs[sel].MetalsDiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsHotGas --- ', nstart+sel
     print, 'MetalsHotGas --- ', LGs[sel].MetalsHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsEjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsEjectedMass --- ', nstart+sel
     print, 'MetalsEjectedMass --- ', LGs[sel].MetalsEjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsICM --- ', nstart+sel
     print, 'MetalsICM --- ', LGs[sel].MetalsICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'PrimordialAccretionRate --- ', nstart+sel
     print, 'PrimordialAccretionRate --- ', LGs[sel].PrimordialAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRadius --- ', nstart+sel
     print, 'CoolingRadius --- ', LGs[sel].CoolingRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate --- ', nstart+sel
     print, 'CoolingRate --- ', LGs[sel].CoolingRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate_beforeAGN --- ', nstart+sel
     print, 'CoolingRate_beforeAGN --- ', LGs[sel].CoolingRate_beforeAGN
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'QuasarAccretionRate --- ', nstart+sel
     print, 'QuasarAccretionRate --- ', LGs[sel].QuasarAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'RadioAccretionRate --- ', nstart+sel
     print, 'RadioAccretionRate --- ', LGs[sel].RadioAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Sfr --- ', nstart+sel
     print, 'Sfr --- ', LGs[sel].Sfr
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'SfrBulge --- ', nstart+sel
     print, 'SfrBulge --- ', LGs[sel].SfrBulge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'XrayLum --- ', nstart+sel
     print, 'XrayLum --- ', LGs[sel].XrayLum
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeSize --- ', nstart+sel
     print, 'BulgeSize --- ', LGs[sel].BulgeSize
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarDiskRadius --- ', nstart+sel
     print, 'StellarDiskRadius --- ', LGs[sel].StellarDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasDiskRadius --- ', nstart+sel
     print, 'GasDiskRadius --- ', LGs[sel].GasDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CosInclination --- ', nstart+sel
     print, 'CosInclination --- ', LGs[sel].CosInclination
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[0] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[1] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[2] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[3] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[4] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[5] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[6] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[7] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[8] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[9] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[10] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[11] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[12] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[13] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[14] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[15] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[16] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[17] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[18] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagDust(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagDust[19] --- ', nstart+sel
     print, 'ObsMagDust --- ', LGs[sel].ObsMagDust(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[0] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[1] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[2] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[3] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[4] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[5] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[6] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[7] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[8] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[9] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[10] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[11] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[12] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[13] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[14] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[15] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[16] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[17] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[18] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMag(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMag[19] --- ', nstart+sel
     print, 'ObsMag --- ', LGs[sel].ObsMag(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[0] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[1] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[2] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[3] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[4] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[5] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[6] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[7] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[8] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[9] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[10] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[11] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[12] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[13] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[14] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[15] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[16] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[17] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[18] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ObsMagBulge(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ObsMagBulge[19] --- ', nstart+sel
     print, 'ObsMagBulge --- ', LGs[sel].ObsMagBulge(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassWeightAge --- ', nstart+sel
     print, 'MassWeightAge --- ', LGs[sel].MassWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'rbandWeightAge --- ', nstart+sel
     print, 'rbandWeightAge --- ', LGs[sel].rbandWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[0] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[1] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[2] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[3] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[4] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[5] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[6] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[7] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[8] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[9] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[10] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[11] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[12] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[13] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[14] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[15] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[16] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[17] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[18] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_DiskMass(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_DiskMass[19] --- ', nstart+sel
     print, 'sfh_DiskMass --- ', LGs[sel].sfh_DiskMass(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[0] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[1] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[2] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[3] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[4] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[5] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[6] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[7] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[8] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[9] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[10] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[11] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[12] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[13] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[14] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[15] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[16] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[17] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[18] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_BulgeMass(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_BulgeMass[19] --- ', nstart+sel
     print, 'sfh_BulgeMass --- ', LGs[sel].sfh_BulgeMass(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[0] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[1] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[2] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[3] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[4] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[5] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[6] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[7] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[8] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[9] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[10] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[11] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[12] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[13] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[14] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[15] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[16] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[17] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[18] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_ICM(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_ICM[19] --- ', nstart+sel
     print, 'sfh_ICM --- ', LGs[sel].sfh_ICM(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[0] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[1] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[2] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[3] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[4] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[5] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[6] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[7] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[8] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[9] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[10] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[11] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[12] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[13] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[14] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[15] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[16] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[17] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[18] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsDiskMass[19] --- ', nstart+sel
     print, 'sfh_MetalsDiskMass --- ', LGs[sel].sfh_MetalsDiskMass(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[0] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[1] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[2] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[3] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[4] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[5] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[6] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[7] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[8] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[9] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[10] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[11] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[12] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[13] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[14] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[15] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[16] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[17] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[18] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsBulgeMass[19] --- ', nstart+sel
     print, 'sfh_MetalsBulgeMass --- ', LGs[sel].sfh_MetalsBulgeMass(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[0] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[1] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[2] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[3] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[4] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[5] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[6] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[7] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[8] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[9] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[10] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[11] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[12] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[13] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[14] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[15] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[16] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[17] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[18] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.sfh_MetalsICM(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'sfh_MetalsICM[19] --- ', nstart+sel
     print, 'sfh_MetalsICM --- ', LGs[sel].sfh_MetalsICM(19)
     badranges=[badranges,sel]
 endif
if(bad) then begin 
     print, 'badranges found: ',badranges
endif
return, badranges
end
