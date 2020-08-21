// Case Notes
//
//
// so i found that for liquids, dpdt is best switched off, also it's better to use sensible internal energy 
//
//
//https://www.cfd-online.com/Forums/openfoam/119851-what-dpdt-term-chtmultiregionfoam.html
//
// https://openfoam.org/release/2-2-0/thermophysical-multiphase-energy/
//
// of course, best is to use buoyantBousinesqPimpleFoam which assumes fluids are incompressible under pressure
//
// (no pressure coupling between energy an pressure equations)
// OpenFOAM based LES of slot jet impingement heat transfer at low nozzle to plate spacing using four SGS models
//
//
// such cases deal with unstable turbulent flow, eg. LES type. 
// for temperature calculations
//
// for new case:
//
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1000;
    nCorrectors     1;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       1e5;
}

relaxationFactors
{
    nuTilda         0.7;
    "(U|h|e|k|epsilon|omega|R|rho)"               0.7;

    ".*Final"        1.0;
}
// this will ensure a standard ish verion of fvSoltuion
//
// now i must change constant/thermophysicalProperties to:
//
//
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       polynomial;
    thermo          hPolynomial;
    equationOfState icoPolynomial;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    // coefficients for DowThermA oil

    specie
    {
        molWeight       166.0;
    }
    equationOfState
    {
        rhoCoeffs<8>    ( 1.3632784512E+03 -1.3938506273E+00 1.7344173768E-03 -1.7393088213E-06 0 0 0 0 );
    // valid from 15C to 405C
    }
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpCoeffs<8>     ( 1.6665337809E+03 -6.4879359625E+00 3.4427492838E-02 -5.5134013220E-05 3.2242982175E-08 0 0 0 );

    // valid from 15C to 405C
    }
    transport
    {
        // valid from 15C to 405C

        muCoeffs<8>     ( 7.8938104591E-01 -9.4695754139E-03  4.7123657786E-05 -1.2421327939E-07 1.8270677458E-10 -1.4211773801E-13 4.5661764167E-17 0 );
        kappaCoeffs<8>  ( 1.8560400000E-01 -1.6000000000E-04 0 0 0 0 0 0 );
    }
}


dpdt            no;
// this failed, i can try turning dpdt back on...
//
//
// this failed...
//
//
// let me try then to limit pimple iterations to 3 outer corrector
//
//
// this also failed
//
//
// all right looks like i need to try buoyantBousinesqPimpleFoam
//
// because this solver seems to be able to solve things even for LES in literature
//
// for buoyantBousinesqPimpleFoam, newtonaian trasnport model and transport properties is used.
//
// constant viscosity, prandtl number, and constant beta
//
//
// i suspect that making some thermophysical properties constant may stabilise the solver...
//
//
// i made the density constant and it works!!!
//
//
// the stability increased leaps and bounds
//
// however, at 1.06s, the max temp is 406 ish and min temp is 238
//
//
// useful links
// boussinesqPimpleFoam solvers:
// Varun Report: buoyantBoussinesqPimpleFoam
//
// http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2016/VarunVenkatesh/Varun_report.pdf
//
//
// https://caefn.com/openfoam/solvers-buoyantpimplefoam
//
//
// buoyantPimpleFoam internal energy (liquid) vs enthalpy (gas) modelling
//
// https://openfoam.org/release/2-2-0/thermophysical-multiphase-energy/
//
//
// what is dpdt and how does it affect equations:
//
// https://www.cfd-online.com/Forums/openfoam/119851-what-dpdt-term-chtmultiregionfoam.html
//
// LES for buoyantBoussinesqPImpleFoam successfully used
//
// https://link-springer-com.libproxy.berkeley.edu/article/10.1007/s00231-018-2470-8?shared-article-renderer
//
// courant Number function for debugging
//
// https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html
//
// https://www.openfoam.com/documentation/user-guide/monitoring-jobs.php
