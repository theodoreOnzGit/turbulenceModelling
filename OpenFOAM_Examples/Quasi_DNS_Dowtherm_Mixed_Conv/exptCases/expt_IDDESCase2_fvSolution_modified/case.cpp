// Case Notes
//
//
//
// When I conducted tests of pimpleFoam using a spalart Allmaras RANS type model, i noticed that the mesh was odd
// 
// I also noted that we had the same problem where the timestep adjusted was getting smaller and smaller, it was a 
// classic blowup case.
//
//
// However, I also found that the case didn't blow up as fast, the timestep was going to 1e-15 and that wasn't causing 
// as bad a blow up, the Co didn't blow up as bad.
//
// My conjecture is this, I noted that i put three inner corrector loops,
// the fvSolution File was like this:
solvers
{
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "(U|k|epsilon|omega|R|nuTilda)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }

    "(U|k|epsilon|omega|R|nuTilda)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }

}
 
PIMPLE
{
    nNonOrthogonalCorrectors 1;
    nCorrectors          3;
    nOuterCorrectors    1;

    residualControl
    {
        U
        {
                tolerance  1e-3;
                relTol      0;
        }
        p
        {
                tolerance  5e-2;
                relTol      0;
        }
     }

}


relaxationFactors
{
    fields
    {
        p      1;
        pFinal   1;
    }
    equations
    {
        "(U|k|epsilon|omega|R|nuTilda)"     1;
        "(U|k|epsilon|omega|R|nuTilda)Final"   1;
    }
}
// but in this case, the fvSolution file was originally:
//
//
solvers
{
    p
    {
        solver          GAMG;
        tolerance       0;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    pFinal
    {
        $p;
        smoother        DICGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }


    p_rgh
    {
        solver           GAMG;
        tolerance        1e-7;
        relTol           0.01;

        smoother         DICGaussSeidel;

    }


    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|h|e|k|epsilon|omega|R|rho)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0.1;
    }

    "(U|h|e|k|epsilon|omega|R|rho)Final"
    {
        $U;
        relTol          0;
    }


    PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       1e5;
}


// one might conjecture that as the mesh size gets ultrasmall,
//
// the absolute tolerance for p and U and other quantities may not be good enough for that mesh size, thus introducing errors
//
// furthermore, since errors come from exploding Co, one might conjecture that the pressure equation is iffy, thus it needs more corrector loops
//
//
//
//this case is basically to adjust fvSolution to have 10 inner corrector loops and have abs tolerance 0, rel Tol 0.01
//
//
//unfortunately test blew up, and earlier at that...
//
//
//#alternative 2
//
// I think doing PIMPLE with Co = 0.5 can help stabilise the solution and dampen out the numerical errors and leave
//
// the turbulence able to do what it needs to do
//
//
//
// so this is how i set up the case:
 PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1000;
    nCorrectors     10;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       1e5;
}

relaxationFactors
{
    nuTilda         0.1;
    "(U|h|e|k|epsilon|omega|R|rho)"               0.1;

    ".*Final"        1.0;
}

// however, the case blows up pretty quickly as the residuals get large...
// At PIMPLE iteration 40, it blows up
// the initial residual velocities tend to blow up after awhile.
// Is it the relaxation factor? IDK, i can try
//
//
//
// I reduced relaxation factor to 0.01, it seems that when the x velocity residual reaches 0.14 ish, the instability continues to grow
//
// at PIMPLE iteration 47, it still blows up, one can also see continutity errors blowing up when that happens
//
// so the instabilities start growing after iteration 24
// PIMPLE: iteration 24
/*
DILUPBiCGStab:  Solving for Ux, Initial residual = 0.141764, Final residual = 5.85663e-07, No Iterations 1
DILUPBiCGStab:  Solving for Uy, Initial residual = 0.0828515, Final residual = 1.04073e-06, No Iterations 1
DILUPBiCGStab:  Solving for Uz, Initial residual = 0.000199656, Final residual = 1.50377e-09, No Iterations 1
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41812e+06
DILUPBiCGStab:  Solving for h, Initial residual = 0.283213, Final residual = 1.52996e-06, No Iterations 1
GAMG:  Solving for p_rgh, Initial residual = 0.0435579, Final residual = 0.000326054, No Iterations 19
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.12333e-06, global = -2.11262e-06, cumulative = -4.48319e-05
GAMG:  Solving for p_rgh, Initial residual = 0.000319712, Final residual = 2.3744e-06, No Iterations 19
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11274e-06, global = -2.11262e-06, cumulative = -4.69445e-05
GAMG:  Solving for p_rgh, Initial residual = 1.37398e-06, Final residual = 7.4328e-08, No Iterations 8
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -4.90571e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45283e-08, Final residual = 7.45283e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -5.11697e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45639e-08, Final residual = 7.45639e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -5.32823e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45643e-08, Final residual = 7.45643e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -5.5395e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45643e-08, Final residual = 7.45643e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -5.75076e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45643e-08, Final residual = 7.45643e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -5.96202e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45643e-08, Final residual = 7.45643e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11263e-06, global = -2.11262e-06, cumulative = -6.17328e-05
GAMG:  Solving for p_rgh, Initial residual = 7.45643e-08, Final residual = 5.51803e-10, No Iterations 16
Pressure gradient source: uncorrected Ubar = 0.22618, pressure gradient = 1.41813e+06
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 2.11262e-06, global = -2.11262e-06, cumulative = -6.38455e-05
*/

/*
PIMPLE: iteration 47
DILUPBiCGStab:  Solving for Ux, Initial residual = 0.150879, Final residual = 4.26877e-05, No Iterations 1
DILUPBiCGStab:  Solving for Uy, Initial residual = 0.104736, Final residual = 1.45675e-05, No Iterations 1
DILUPBiCGStab:  Solving for Uz, Initial residual = 0.240588, Final residual = 1.21229e-05, No Iterations 1
Pressure gradient source: uncorrected Ubar = 0.225788, pressure gradient = 2.50089e+07
DILUPBiCGStab:  Solving for h, Initial residual = 0.886997, Final residual = 5.24155e-06, No Iterations 1
GAMG:  Solving for p_rgh, Initial residual = 0.511419, Final residual = 0.00426414, No Iterations 21
Pressure gradient source: uncorrected Ubar = 0.224482, pressure gradient = 9.32367e+07
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.334255, global = -0.331302, cumulative = -1.26678
GAMG:  Solving for p_rgh, Initial residual = 0.00194401, Final residual = 1.01674e-05, No Iterations 13
Pressure gradient source: uncorrected Ubar = 0.22388, pressure gradient = 1.24653e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331359, global = -0.331302, cumulative = -1.59808
GAMG:  Solving for p_rgh, Initial residual = 2.10545e-05, Final residual = 1.93901e-07, No Iterations 12
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.25102e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331303, global = -0.331302, cumulative = -1.92939
GAMG:  Solving for p_rgh, Initial residual = 2.38575e-07, Final residual = 8.40815e-08, No Iterations 2
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -2.26069
GAMG:  Solving for p_rgh, Initial residual = 8.37245e-08, Final residual = 8.37245e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -2.59199
GAMG:  Solving for p_rgh, Initial residual = 8.37183e-08, Final residual = 8.37183e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -2.92329
GAMG:  Solving for p_rgh, Initial residual = 8.37182e-08, Final residual = 8.37182e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -3.25459
GAMG:  Solving for p_rgh, Initial residual = 8.37182e-08, Final residual = 8.37182e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -3.58589
GAMG:  Solving for p_rgh, Initial residual = 8.37182e-08, Final residual = 8.37182e-08, No Iterations 0
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.2511e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -3.9172
GAMG:  Solving for p_rgh, Initial residual = 8.37182e-08, Final residual = 8.05647e-10, No Iterations 21
Pressure gradient source: uncorrected Ubar = 0.223871, pressure gradient = 1.25104e+08
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
time step continuity errors : sum local = 0.331302, global = -0.331302, cumulative = -4.2485
*/

// it seems that h is also unable to iterate properly, the initial residual doesn't go below 0.225
