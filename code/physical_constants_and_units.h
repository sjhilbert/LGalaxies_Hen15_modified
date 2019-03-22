/*  Copyright (C) <2016+>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/** @file   physical_constants_and_units.h
 *
 *  @author Stefan Hilbert
 *  @date   2018
 *
 *  @brief  Sets up some variables used to convert from physical to internal
 *          units obtained from
 *          UNITLENGTH_IN_CM (cm to Mpc), 
 *          UNITMASS_IN_G (g to 1e10Msun), and
 *          UNITVELOCITY_IN_CM_PER_S (cm/s to km/s).
 *          
 *          \f$\rm{UnitLength}_{\rm{cm}}= 3.08568\times 10^{24}\rm{cm}\f$,
 *          converts from cm into Mpc, and
 *          \f$\rm{UnitVelocity}_{\rm{cm/s}}=10000\rm{cm/s}\f$
 *          converts from cm/s to Km/s (cm to km). 
 *          \f$\rm{UnitTime}_{\rm{s}}\f$
 *          is derived from these two quantities:
 *          
 *          \f$\frac{\rm{UnitLength}_{\rm{cm}}}{\rm{UnitVelocity}_{\rm{cm/s}}}
 *          =3.08568\times 10^{19}\rm{Mpc}~\rm{Km}^{-1}\rm{s}^{-1}\f$,
 *          
 *          through the code \f$t_{\rm{dyn}}\f$ has internal units and its never
 *          converted (note that \f$t_{\rm{dyn}}\f$ has an h factor, as the code internal
 *          units which despite not being included is \f$\rm{UnitTime}_{\rm{s}}\f$ is
 *          included in the output of time_to_present() - so it is consistent).
 *          
 *          \f$ \rm{UnitDensity}_{\rm{cgs}} =
 *          \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}}^3}=6.769898\times 10^{-31}\f$,
 *          converts density in \f$\rm{g}~\rm{cm}^{-3}\f$ into internal units
 *          \f$(10^{10}M_{\odot}\rm{Mpc}^{-3})\f$
 *          
 *          \f$ \rm{UnitPressure}_{\rm{cgs}} =
 *          \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}} \times \rm{UnitTime}_{\rm{s}}^2}
 *          =6.769898\times 10^{-21}\f$, converts pressure in
 *          \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-2}\f$ into internal units
 *          \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s) \f$
 *          
 *          \f$ \rm{UnitCoolingRate}_{\rm{cgs}} =
 *          \frac{\rm{UnitPressure}_{\rm{cgs}}}{\rm{UnitTime}_{\rm{s}}}=2.193973\times 10^{-40}\f$,
 *          converts the cooling rate in \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-3}\f$ into
 *          internal units \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s)^{-3}) \f$
 *          
 *          \f$ \rm{UnitEnergy}_{\rm{cgs}} =
 *          \frac{\rm{UnitMass}_{\rm{g}} \times \rm{UnitLength}_{\rm{cm}}^2}{\rm{UnitTime}_{\rm{s}}^2}
 *          =1.989000\times 10^{53}\f$, converts energy in
 *          \f$\rm{g}~\rm{cm}^2\rm{s}^{-2}\f$ into internal units
 *          \f$(10^{10}M_{\odot}~\rm{Mpc}^{2}(Mpc/Mk/s)^{-2})\f$
 *          
 *          \f$ \rm{Hubble} = \rm{HUBBLE} \times \rm{UnitTime}_{\rm{s}}=100.0001\f$, where
 *          \f$\rm{HUBBLE}=3.2407789\times 10^{-18} h~\rm{s}^{-1}\f$, is the hubble
 *          constante in \f$(h~\rm{Km}~\rm{s}^{-1}\rm{Mpc}^{-1})\f$.
 *
 *
 * @note  the values for physical constants and basic internal units are updated, and differ
 *        slightly from those in the public Henriques 2015 version of the code.
 *        for using the old public Henriques 2015 values, comment the new values and uncomment
 *        the old values below.
 **/
#ifndef PHYSICAL_CONSTANTS_AND_UNITS_H
#define PHYSICAL_CONSTANTS_AND_UNITS_H

// /* physical constants: in external (CGS or astron.) units: */
// #define GRAVITY                           6.67408e-8        /* in cgs */
// #define SOLAR_MASS                        1.98855e33        /* in cgs */
// #define SOLAR_LUM                         3.826e33          /* in cgs */
// #define RAD_CONST                         7.565e-15         /* in cgs */
// #define AVOGADRO                          6.02214086e23     /* in cgs */
// #define BOLTZMANN                         1.38064852e-16    /* in cgs */
// #define GAS_CONST                         8.3142598e7       /* in cgs */
// #define SPEED_OF_LIGHT                    2.99792458e10     /* in cgs */
// #define PLANCK                            6.62607004e-27    /* in cgs */
// #define PROTONMASS                        1.67262190e-24    /* in cgs */
// #define HUBBLE                            3.24077893e-18    /* in h/sec */
// #define D_HUBBLE                          2997.92458        /* in Mpc/h */
// 
// 
// /* standard unit conversions: */
// #define SEC_PER_YEAR                      3.1536000e7
// #define SEC_PER_MEGAYEAR                  3.1536000e13
// #define SEC_PER_GIGAYEAR                  3.1536000e16
// #define LENGTH_10_PC_IN_CM                3.08567758e19
// 
// 
// /* basic internal units: in terms of external units (ignoring Hubble_h for masses and lengths): */
// #define UNITLENGTH_IN_CM                  3.08567758e+24    /* Mpc - WATCH OUT, distances in the code are in Mpc/h */
// #define UNITMASS_IN_G                     1.98855e+43	      /* 10^10Msun - WATCH OUT, masses in the code are in 10^10Msun/h */
// #define UNITVELOCITY_IN_CM_PER_S          1.e5              /* Km/s - WATCH OUT, this are the correct units in the code km/s */



/* for reference: 
 * these are the old values used by the public Henriques 2015 version
 * if you prefer these, uncomment them and comment the new values above
 */
#define GRAVITY                           6.672e-8
#define SOLAR_MASS                        1.989e33
#define SOLAR_LUM                         3.826e33
#define RAD_CONST                         7.565e-15
#define AVOGADRO                          6.0222e23
#define BOLTZMANN                         1.3806e-16
#define GAS_CONST                         8.31425e7
#define SPEED_OF_LIGHT                    2.9979e10
#define PLANCK                            6.6262e-27
#define PROTONMASS                        1.6726e-24
#define HUBBLE                            3.2407789e-18
#define D_HUBBLE                          2997.92458 


#define SEC_PER_YEAR                      3.155e7
#define SEC_PER_MEGAYEAR                  3.155e13
#define SEC_PER_GIGAYEAR                  3.155e16
#define LENGTH_10_PC_IN_CM                3.08567758e19


#define UNITLENGTH_IN_CM                   3.08568e+24
#define UNITMASS_IN_G                      1.989e+43
#define UNITVELOCITY_IN_CM_PER_S           100000


#ifdef DERIVED_UNITS_AS_STATIC_CONSTANTS

/* derived internal untits: in terms of external units (ignoring Hubble_h for masses and lengths): */
// time: 
static const double UnitTime_in_s                    = (UNITLENGTH_IN_CM / UNITVELOCITY_IN_CM_PER_S);
static const double UnitTime_in_years                = (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_YEAR));
static const double UnitTime_in_Megayears            = (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_MEGAYEAR));
static const double UnitTime_in_Gigayears            = (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_GIGAYEAR));

// converts g.cm^-3 into internal units (1e10Msun Mpc^-3)
// static const double UnitDensity_in_cgs               = UNITMASS_IN_G / pow3(UNITLENGTH_IN_CM);//6.769898e-31
static const double UnitDensity_in_cgs               = (UNITMASS_IN_G / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)); //6.769898e-31

// converts g.cm^-1s^-2 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-2) \f$
// static const double UnitPressure_in_cgs              = UNITMASS_IN_G / UNITLENGTH_IN_CM / pow2(UnitTime_in_s);//6.769898e-21
static const double UnitPressure_in_cgs              = (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)); //6.769898e-21

// converts g.cm^-1.s^-3 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-3)
// static const double UnitCoolingRate_in_cgs           = UnitPressure_in_cgs / UnitTime_in_s;//2.193973e-40
static const double UnitCoolingRate_in_cgs           = (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)); //2.193973e-40

// converts g.cm^2.s^-2 into internal units (10^10Msun.Mpc^2(Mpc/Km/s)^-2)
// static const double UnitEnergy_in_cgs                = UNITMASS_IN_G * pow2(UNITLENGTH_IN_CM) / pow2(UnitTime_in_s);//1.989000e+53
static const double UnitEnergy_in_cgs                = (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S); //1.989000e+53


/* phyical constants: in internal units (ignoring Hubble_h for masses and lengths): */
// gravity in internal units
// static const double Gravity                          = GRAVITY / pow3(UNITLENGTH_IN_CM) * UNITMASS_IN_G * pow2(UnitTime_in_s);//43.00708
static const double Gravity                          = (GRAVITY * UNITMASS_IN_G / (UNITLENGTH_IN_CM * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S)); //43.00708

// speed of light in internal units
// static const double SpeedOfLight                     = SPEED_OF_LIGHT / UNITVELOCITY_IN_CM_PER_S;
static const double SpeedOfLight                     = (SPEED_OF_LIGHT / UNITVELOCITY_IN_CM_PER_S); // 2.9979e5 (km/s)

//solar mass in internal units
static const double SolarMass                        = (SOLAR_MASS / UNITMASS_IN_G); // 1e-10 (10^10 Msolar)

// converts the Hubble constant from h.s^-1 into h.Km.s^-1.Mpc-1:
// static const double Hubble                           = HUBBLE * UnitTime_in_s;//100.000
static const double Hubble                           = (HUBBLE * UNITLENGTH_IN_CM / UNITVELOCITY_IN_CM_PER_S); //100.000

// critical density of the universe in internal units
// static const double RhoCrit                          = 3 * Hubble * Hubble / (8 * M_PI * Gravity);//27.75505 (h^2.10^10Msun.Mpc^-3)
static const double RhoCrit                          = (3 * HUBBLE * HUBBLE * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM / (8 * M_PI * GRAVITY * UNITMASS_IN_G)); //27.75505 (h^2.10^10Msun.Mpc^-3)

#else  /* not defined DERIVED_UNITS_AS_STATIC_CONSTANTS */
/* use macros for derived units */

/* derived internal untits: in terms of external units: (ignoring Hubble_h for masses and lengths)*/
// time: 
#define UnitTime_in_s                      (UNITLENGTH_IN_CM / UNITVELOCITY_IN_CM_PER_S)
#define UnitTime_in_years                  (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_YEAR))
#define UnitTime_in_Megayears              (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_MEGAYEAR))
#define UnitTime_in_Gigayears              (UNITLENGTH_IN_CM / (UNITVELOCITY_IN_CM_PER_S * SEC_PER_GIGAYEAR))

//converts g.cm^-3 into internal units (1e10Msun Mpc^-3)
#define UnitDensity_in_cgs                 (UNITMASS_IN_G / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)) //6.769898e-31

//converts g.cm^-1s^-2 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-2) \f$
// #define UnitPressure_in_cgs                (UNITMASS_IN_G / UNITLENGTH_IN_CM / ((UnitTime_in_s) * (UnitTime_in_s))) //6.769898e-21
#define UnitPressure_in_cgs                (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)) //6.769898e-21

//converts g.cm^-1.s^-3 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-3)
// #define UnitCoolingRate_in_cgs             ((UnitPressure_in_cgs) / (UnitTime_in_s)) //2.193973e-40
#define UnitCoolingRate_in_cgs             (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM)) //2.193973e-40

//converts g.cm^2.s^-2 into internal units (10^10Msun.Mpc^2(Mpc/Km/s)^-2)
// #define UnitEnergy_in_cgs                  (UNITMASS_IN_G * (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM) / ((UnitTime_in_s) * (UnitTime_in_s))) //1.989000e+53
#define UnitEnergy_in_cgs                  (UNITMASS_IN_G * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S) //1.989000e+53


/* phyical constants: in internal units (without Hubble_h): */
//gravity in internal units
// #define Gravity                            (GRAVITY / (UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM) * UNITMASS_IN_G * ((UnitTime_in_s) * (UnitTime_in_s))) //43.00708
#define Gravity                            (GRAVITY * UNITMASS_IN_G / (UNITLENGTH_IN_CM * UNITVELOCITY_IN_CM_PER_S * UNITVELOCITY_IN_CM_PER_S)) //43.00708

//solar mass in internal units
#define SolarMass                          (SOLAR_MASS / UNITMASS_IN_G) // 1e-10 (10^10 Msolar)

//speed of light in internal units
#define SpeedOfLight                       (SPEED_OF_LIGHT / UNITVELOCITY_IN_CM_PER_S) // 2.9979e5 (km/s)

// converts the Hubble constant from h.s^-1 into h.Km.s^-1.Mpc-1 
// #define Hubble                             (HUBBLE * (UnitTime_in_s)) //100.000
#define Hubble                             (HUBBLE * UNITLENGTH_IN_CM / UNITVELOCITY_IN_CM_PER_S) //100.000

// critical density of the universe in internal units
// #define RhoCrit                            (3 * (Hubble) * (Hubble) / (8 * M_PI * (Gravity))) //27.75505 (h^2.10^10Msun.Mpc^-3)
#define RhoCrit                            (3 * HUBBLE * HUBBLE * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM * UNITLENGTH_IN_CM / (8 * M_PI * GRAVITY * UNITMASS_IN_G)) //27.75505 (h^2.10^10Msun.Mpc^-3)

#endif  /* not defined DERIVED_UNITS_AS_STATIC_CONSTANTS */

#endif /* header guard */
