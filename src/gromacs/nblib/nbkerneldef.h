#ifndef GROMACS_NBKERNELDEF_H
#define GROMACS_NBKERNELDEF_H

//! Enum for selecting the SIMD kernel type
enum class BenchMarkKernels : int
{
    SimdAuto,
    SimdNo,
    Simd4XM,
    Simd2XMM,
    Count
};

//! Enum for selecting the combination rule
enum class BenchMarkCombRule : int
{
    RuleGeom,
    RuleLB,
    RuleNone,
    Count
};

//! Enum for selecting coulomb type
enum class BenchMarkCoulomb : int
{
    Pme,
    ReactionField,
    Count
};


#endif // GROMACS_NBKERNELDEF_H
