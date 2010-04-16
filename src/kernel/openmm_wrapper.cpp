#include <cmath>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cctype>
#include <algorithm>

using namespace std;

#include "OpenMM.h"

#include "gmx_fatal.h"
#include "typedefs.h"
#include "mdrun.h"
#include "physics.h"
#include "string2.h"
#include "gmx_gpu_utils.h"
#include "mtop_util.h"
#include "warninp.h"

#include "openmm_wrapper.h"

using namespace OpenMM;

#define MEM_ERR_MSG(str) \
    "The %s-simulation GPU memory test detected errors. As memory errors would cause incorrect " \
    "simulation results, gromacs has aborted execution.\n Make sure that your GPU's memory is not " \
    "overclocked and that the device is properly cooled.\n", (str)


/** Convert string to type T. */
template <class T>
static bool from_string(T& t, const string& s, ios_base& (*f)(ios_base&))
{
    istringstream iss(s);
    return !(iss >> f >> t).fail();
}

/**
 * Split string around a given delimiter.
 */
static vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        if (item.length() != 0)
            elems.push_back(item);
    }
    return elems;
}

/**
 * Split a string "option=value" into "option" and "value" strings.
 */
static void split_opt_val(const string &s, string &opt, string &val)
{
    size_t eqPos = s.find('=');
    if (eqPos != string::npos)
    {
        opt = s.substr(0, eqPos);
        if (eqPos != s.length())  val = s.substr(eqPos+1);
    }
}

/**
 * Compare two strings ignoring case.
 */
static bool isStringEqNCase(const string s1, const string s2)
{
    return (gmx_strncasecmp(s1.c_str(), s2.c_str(), max(s1.length(), s2.length())) == 0);
}

/**
 * Convert string to upper case.
 */
static string toUpper(const string &s)
{
    string stmp(s);
    std::transform(stmp.begin(), stmp.end(), stmp.begin(), static_cast < int(*)(int) > (toupper));
    return stmp;
}

#define SIZEOF_PLATFORMS    1
#define SIZEOF_MEMTESTS     3
#define SIZEOF_DEVICEIDS    1
#define SIZEOF_FORCE_DEV    2

/** Possible options of the mdrun -device parameter. */
static const char *devOptStrings[] = { "platform", "deviceid", "memtest", "force-device" };
enum devOpt
{
    PLATFORM     = 0,
    DEVICEID     = 1,
    MEMTEST      = 2,
    FORCE_DEVICE = 3
};

/**
 * Parse the options string provided to mdrun.
 */
class GmxOpenMMPlatformOptions
{
public:
    GmxOpenMMPlatformOptions(const char *optionString);
//    GmxOpenMMPlatformOptions(string &optionString);
    ~GmxOpenMMPlatformOptions() { options.clear(); }
    string getOptionValue(const string &optName);
    void remOption(const string &optName) { options.erase(toUpper(optName)); }
private:
    void setOption(const string &opt, const string &val);
    map<string, string> options;
    static const char * const platforms[SIZEOF_PLATFORMS];
    static const char * const memtests[SIZEOF_MEMTESTS];
    static const char * const deviceid[SIZEOF_DEVICEIDS];
    static const char * const force_dev[SIZEOF_FORCE_DEV];
};

/* possible values of the parameters that will not be passed to OMM
   first value will always be the default used in case if the option is not specified */
const char * const GmxOpenMMPlatformOptions::platforms[SIZEOF_PLATFORMS] = { "Cuda" /*,"OpenCL"*/ };
const char * const GmxOpenMMPlatformOptions::memtests[SIZEOF_MEMTESTS]   = { "15", "full", "off" }; /* also valid any positive integer >=15 */
const char * const GmxOpenMMPlatformOptions::deviceid[SIZEOF_DEVICEIDS]  = { "0" };       /* also valid any positive integer */
const char * const GmxOpenMMPlatformOptions::force_dev[SIZEOF_FORCE_DEV] = { "no", "yes" };

/*
 * TODO doc
 */
GmxOpenMMPlatformOptions::GmxOpenMMPlatformOptions(const char *optionString)
{
    // set default values
    setOption("platform",       platforms[0]);
    setOption("memtest",        memtests[0]);
    setOption("deviceid",       deviceid[0]);
    setOption("force-device",   force_dev[0]);

    string opt(optionString);

    // remove all whitespaces
    opt.erase(remove_if(opt.begin(), opt.end(), ::isspace), opt.end());
    // tokenize around ","-s
    vector<string> tokens = split(opt, ',');

    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it)
    {
        string opt = "", val = "";
        split_opt_val(*it, opt, val);

        if (isStringEqNCase(opt, "platform"))
        {
            setOption(opt, val);
            continue;
        }

        if (isStringEqNCase(opt, "memtest"))
        {
            if (!isStringEqNCase(val, "full") && !isStringEqNCase(val, "off")) /* the value has to be an integer >15 */
            {
                int secs;
                if (!from_string<int>(secs, val, std::dec))
                {
                    gmx_fatal(FARGS, "Invalid value for option memtestoption: \"%s\"!\n", val.c_str());
                }
                if (secs < 15)
                {
                    gmx_fatal(FARGS, "Incorrect value for memtest option (%d). "
                            "Memtest needs to run for at least 15s!\n", secs);
                }
            }
            setOption(opt, val);
            continue;
        }

        if (isStringEqNCase(opt, "deviceid"))
        {
            int id;
            if (!from_string<int>(id, val, std::dec) )
            {
                gmx_fatal(FARGS, "Invalid device id: \"%s\"!\n", val.c_str());
            }
            setOption(opt, val);
            continue;
        }

        if (isStringEqNCase(opt, "force-device"))
        {
            if (!isStringEqNCase(val, "yes") && !isStringEqNCase(val, "no"))
            {
                gmx_fatal(FARGS, "Invalid OpenMM force option: \"%s\"!\n", val.c_str());
            }
            setOption(opt, val);
            continue;
        }

        // if we got till here something went wrong
        gmx_fatal(FARGS, "Invalid OpenMM platform option: \"%s\"!\n", (*it).c_str());
    }

    /*    cout << "== platform = " << getOptionValue("platform") << endl
             << "== deviceID = " << getOptionValue("deviceid") << endl
             << "== memtest  = " << getOptionValue("memtest") << endl
             << "== force    = " << getOptionValue("force-device") << endl;
    */
}

/*
 * TODO doc
 */
string GmxOpenMMPlatformOptions::getOptionValue(const string &optName)
{
    if (options.find(toUpper(optName)) != options.end())
    {
        return options[toUpper(optName)];
    }
    else
    {
        return NULL;
    }
}

/*
 * TODO doc
 */
void GmxOpenMMPlatformOptions::setOption(const string &opt, const string &val)
{
    options[toUpper(opt)] = val;
}

/*
 * TODO doc
 */
class OpenMMData
{
public:
    System* system;
    Context* context;
    Integrator* integrator;
    bool removeCM;
    GmxOpenMMPlatformOptions *platformOpt;
};

/*
 * TODO doc
 */
void run_memtest(FILE* fplog, int devId, const char* pre_post, GmxOpenMMPlatformOptions *opt)
{

    char warn_buf[STRLEN];

    int which_test;
    int res;
    const char * test_type = opt->getOptionValue("memtest").c_str();
    if (!gmx_strcasecmp(test_type, "off"))
    {
        which_test = 0;
    }
    else
    {
        if (!gmx_strcasecmp(test_type, "full"))
        {
            which_test = 2;
        }
        else
        {
            from_string<int>(which_test, test_type, std::dec);
        }
    }

    switch (which_test)
    {
        case 0: /* no memtest */
            sprintf(warn_buf, "%s-simulation GPU memtest skipped. Note, that faulty memory can cause "
                "incorrect results!\n", pre_post);
            fprintf(fplog, "%s", warn_buf);
            gmx_warning(warn_buf);
            break; /* case 0 */

        case 1: /* quick memtest */
            fprintf(fplog,  "%s-simulation %s GPU memtest in progress...\n", pre_post, test_type);
            fprintf(stdout, "\n%s-simulation %s GPU memtest in progress...", pre_post, test_type);
            fflush(fplog);
            fflush(stdout);
            if (do_quick_memtest(-1) != 0)
            {
                gmx_fatal(FARGS,MEM_ERR_MSG(pre_post));
            }
            else
            {
                fprintf(fplog,  "Memory test completed without errors.\n");
                fprintf(stdout, "done, no errors detected\n");
            }
            break; /* case 1 */

        case 2: /* full memtest */
            fprintf(fplog,  "%s-simulation %s memtest in progress...\n", pre_post, test_type);
            fprintf(stdout, "\n%s-simulation %s memtest in progress...", pre_post, test_type);
            fflush(fplog);
            fflush(stdout);
            if (do_full_memtest(-1) != 0)
            {
                gmx_fatal(FARGS, MEM_ERR_MSG(pre_post) );

            }
            else
            {
                fprintf(fplog, "Memory test completed without errors.\n");
                fprintf(stdout, "done, no errors detected\n");
            }
            break; /* case 2 */

        default: /* timed memtest */
            fprintf(fplog,  "%s-simulation memtest for ~%ds in progress...\n", pre_post, which_test);
            fprintf(stdout, "\n%s-simulation memtest for ~%ds in progress...", pre_post, which_test);
            fflush(fplog);
            fflush(stdout);
            if (do_timed_memtest(-1, which_test) != 0)
            {
                gmx_fatal(FARGS, MEM_ERR_MSG(pre_post));

            }
            else
            {
                fprintf(fplog, "Memory test completed without errors.\n");
                fprintf(stdout, "done, no errors detected.\n");
            }
        }
        fflush(fplog);
        fflush(stdout);
}

/*
 * TODO doc
 */
void checkGmxOptions(t_inputrec *ir, gmx_localtop_t *top)
{

    char warn_buf[STRLEN];

    // Abort if unsupported critical options are present

    /* Integrator */
    if (ir->eI ==  eiMD)
        gmx_warning( "OpenMM does not support leap-frog, will use velocity-verlet integrator.\n");

    if ( (ir->eI !=  eiMD)   &&
            (ir->eI !=  eiVV)   &&
            (ir->eI !=  eiVVAK) &&
            (ir->eI !=  eiSD1)  &&
            (ir->eI !=  eiSD2)  &&
            (ir->eI !=  eiBD) )
    {
        gmx_fatal(FARGS, "OpenMM supports only the following integrators: md/md-vv/md-vv-avek, sd/sd1, and bd.\n");
    }

    /* Electroctstics */
    if ( (ir->coulombtype != eelPME) &&
            (ir->coulombtype != eelRF) &&
            (ir->coulombtype != eelEWALD) &&
            // no-cutoff
            ( !(ir->coulombtype == eelCUT && ir->rcoulomb == 0 &&  ir->rvdw == 0)) )
    {
        gmx_fatal(FARGS,"OpenMM supports only the following methods for electrostatics: "
                "NoCutoff (i.e. rcoulomb = rvdw = 0 ),Reaction-Field, Ewald or PME.\n");
    }

    if ( (ir->etc != etcNO) &&
            (ir->eI !=  eiSD1)  &&
            (ir->eI !=  eiSD2)  &&
            (ir->eI !=  eiBD) )
    {
        gmx_warning("OpenMM supports only Andersen thermostat with the md/md-vv/md-vv-avek integrators.\n");
    }

    if (ir->opts.ngtc > 1)
        gmx_fatal(FARGS,"OpenMM does not support multiple temperature coupling groups.\n");

    if (ir->epc != etcNO)
        gmx_fatal(FARGS,"OpenMM does not support pressure coupling.\n");

    if (ir->opts.annealing[0])
        gmx_fatal(FARGS,"OpenMM does not support simulated annealing.\n");

    if (ir->eConstrAlg != econtSHAKE)
        gmx_warning("Constraints in OpenMM are done by a combination "
                    "of SHAKE, SETTLE and CCMA. Accuracy is based on the SHAKE tolerance set "
                    "by the \"shake_tol\" option.\n");

    if (ir->nwall != 0)
        gmx_fatal(FARGS,"OpenMM does not support walls.\n");

    if (ir->ePull != epullNO)
        gmx_fatal(FARGS,"OpenMM does not support pulling.\n");

    if (top->idef.il[F_DISRES].nr > 0)
        gmx_fatal(FARGS,"OpenMM does not support distant restraints.\n");

    if (top->idef.il[F_ORIRES].nr > 0)
        gmx_fatal(FARGS,"OpenMM does not support orientation restraints.\n");

    if (top->idef.il[F_ANGRES].nr > 0)
        gmx_fatal(FARGS,"OpenMM does not support angle restraints.\n");

    if (top->idef.il[F_DIHRES].nr > 0)
        gmx_fatal(FARGS,"OpenMM does not support dihedral restraints.\n");

    if (ir->efep != efepNO)
        gmx_fatal(FARGS,"OpenMM does not support free energy calculations.\n");

    if (ir->opts.ngacc > 1)
        gmx_fatal(FARGS,"OpenMM does not support non-equilibrium MD (accelerated groups).\n");

    if (IR_ELEC_FIELD(*ir))
        gmx_fatal(FARGS,"OpenMM does not support electric fields.\n");

    if (ir->bQMMM)
        gmx_fatal(FARGS,"OpenMM does not support QMMM calculations.\n");

    if (ir->rcoulomb != ir->rvdw)
        gmx_fatal(FARGS,"OpenMM uses a single cutoff for both Coulomb "
                  "and VdW interactions. Please set rcoulomb equal to rvdw.\n");

}

/**
 * Initialize OpenMM, and return a pointer to the OpenMMData object containing 
 * the System, Integrator, and Context.
 *//*
 * TODO doc
 */
void* openmm_init(FILE *fplog, const char *platformOptStr,
                  t_commrec *cr,t_inputrec *ir,
                  gmx_mtop_t *top_global, gmx_localtop_t *top,
                  t_mdatoms *mdatoms, t_forcerec *fr,t_state *state)
{

    char warn_buf[STRLEN];

    static bool hasLoadedPlugins = false;
    string usedPluginDir;

    try
    {
        if (!hasLoadedPlugins)
        {
            vector<string> loadedPlugins;
            /*  Look for OpenMM plugins at various locations (listed in order of priority):
                - on the path in OPENMM_PLUGIN_DIR environment variable if this is specified
                - on the path in the OPENMM_PLUGIN_DIR macro that is set by the build script
                - at the default location assumed by OpenMM
            */
            /* env var */
            char *pluginDir = getenv("OPENMM_PLUGIN_DIR");
            trim(pluginDir);
            /* no env var or empty */
            if (pluginDir != NULL && *pluginDir != '\0')
            {
                loadedPlugins = Platform::loadPluginsFromDirectory(pluginDir);
                if (loadedPlugins.size() > 0)
                {
                    hasLoadedPlugins = true;
                    usedPluginDir = pluginDir;
                }
                else
                {
                    gmx_fatal(FARGS, "The directory provided in the OPENMM_PLUGIN_DIR environment variable "
                              "(%s) does not contain valid OpenMM plugins. Check your OpenMM installation!", 
                              pluginDir);
                }
            }

            /* macro set at build time  */
#ifdef OPENMM_PLUGIN_DIR
            if (!hasLoadedPlugins)
            {
                loadedPlugins = Platform::loadPluginsFromDirectory(OPENMM_PLUGIN_DIR);
                if (loadedPlugins.size() > 0)
                {
                    hasLoadedPlugins = true;
                    usedPluginDir = OPENMM_PLUGIN_DIR;
                }
            }
#endif
            /* default loocation */
            if (!hasLoadedPlugins)
            {
                loadedPlugins = Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
                if (loadedPlugins.size() > 0)
                {
                    hasLoadedPlugins = true;
                    usedPluginDir = Platform::getDefaultPluginsDirectory();
                }
            }

            /* if there are still no plugins loaded there won't be any */
            if (!hasLoadedPlugins)
            {
                gmx_fatal(FARGS, "No OpenMM plugins were found! You can provide the"
                          " plugin directory in the OPENMM_PLUGIN_DIR environment variable.", pluginDir);
            }

            fprintf(fplog, "\nPlugins loaded from directory %s:\t", usedPluginDir.c_str());
            for (int i = 0; i < loadedPlugins.size(); i++)
            {
                fprintf(fplog, "%s, ", loadedPlugins[i].c_str());
            }
            fprintf(fplog, "\n");
        }

        /* parse option string */
        GmxOpenMMPlatformOptions *opt = new GmxOpenMMPlatformOptions(platformOptStr);

        /* check wheter Gromacs options compatibility with OpenMM */
        checkGmxOptions(ir, top);

        /* check GPU support */
        char s[1000];
        is_supported_cuda_gpu(atoi(opt->getOptionValue("deviceid").c_str()), s);

        // Create the system.
        const t_idef& idef = top->idef;
        const int numAtoms = top_global->natoms;
        const int numConstraints = idef.il[F_CONSTR].nr/3;
        const int numSettle = idef.il[F_SETTLE].nr/2;
        const int numBonds = idef.il[F_BONDS].nr/3;
        const int numAngles = idef.il[F_ANGLES].nr/4;
        const int numPeriodic = idef.il[F_PDIHS].nr/5;
        const int numRB = idef.il[F_RBDIHS].nr/5;
        const int num14 = idef.il[F_LJ14].nr/3;
        System* sys = new System();
        if (ir->nstcomm > 0)
            sys->addForce(new CMMotionRemover(ir->nstcomm));

        // Set bonded force field terms.
        const int* bondAtoms = (int*) idef.il[F_BONDS].iatoms;
        HarmonicBondForce* bondForce = new HarmonicBondForce();
        sys->addForce(bondForce);
        int offset = 0;
        for (int i = 0; i < numBonds; ++i)
        {
            int type = bondAtoms[offset++];
            int atom1 = bondAtoms[offset++];
            int atom2 = bondAtoms[offset++];
            bondForce->addBond(atom1, atom2,
                               idef.iparams[type].harmonic.rA, idef.iparams[type].harmonic.krA);
        }
        const int* angleAtoms = (int*) idef.il[F_ANGLES].iatoms;
        HarmonicAngleForce* angleForce = new HarmonicAngleForce();
        sys->addForce(angleForce);
        offset = 0;
        for (int i = 0; i < numAngles; ++i)
        {
            int type = angleAtoms[offset++];
            int atom1 = angleAtoms[offset++];
            int atom2 = angleAtoms[offset++];
            int atom3 = angleAtoms[offset++];
            angleForce->addAngle(atom1, atom2, atom3, 
                    idef.iparams[type].harmonic.rA*M_PI/180.0, idef.iparams[type].harmonic.krA);
        }
        const int* periodicAtoms = (int*) idef.il[F_PDIHS].iatoms;
        PeriodicTorsionForce* periodicForce = new PeriodicTorsionForce();
        sys->addForce(periodicForce);
        offset = 0;
        for (int i = 0; i < numPeriodic; ++i)
        {
            int type = periodicAtoms[offset++];
            int atom1 = periodicAtoms[offset++];
            int atom2 = periodicAtoms[offset++];
            int atom3 = periodicAtoms[offset++];
            int atom4 = periodicAtoms[offset++];
            periodicForce->addTorsion(atom1, atom2, atom3, atom4,
                                      idef.iparams[type].pdihs.mult,
                                      idef.iparams[type].pdihs.phiA*M_PI/180.0, 
                                      idef.iparams[type].pdihs.cpA);
        }
        const int* rbAtoms = (int*) idef.il[F_RBDIHS].iatoms;
        RBTorsionForce* rbForce = new RBTorsionForce();
        sys->addForce(rbForce);
        offset = 0;
        for (int i = 0; i < numRB; ++i)
        {
            int type = rbAtoms[offset++];
            int atom1 = rbAtoms[offset++];
            int atom2 = rbAtoms[offset++];
            int atom3 = rbAtoms[offset++];
            int atom4 = rbAtoms[offset++];
            rbForce->addTorsion(atom1, atom2, atom3, atom4,
                                idef.iparams[type].rbdihs.rbcA[0], idef.iparams[type].rbdihs.rbcA[1],
                                idef.iparams[type].rbdihs.rbcA[2], idef.iparams[type].rbdihs.rbcA[3],
                                idef.iparams[type].rbdihs.rbcA[4], idef.iparams[type].rbdihs.rbcA[5]);
        }

        // Set nonbonded parameters and masses.

        int ntypes = fr->ntype;
        int* types = mdatoms->typeA;
        real* nbfp = fr->nbfp;
        real* charges = mdatoms->chargeA;
        real* masses = mdatoms->massT;
        NonbondedForce* nonbondedForce = new NonbondedForce();
        sys->addForce(nonbondedForce);

        if (ir->rcoulomb == 0)
        {
            nonbondedForce->setNonbondedMethod(NonbondedForce::NoCutoff);
        }
        else
        {
            switch (ir->coulombtype)
            {
            case eelRF: // TODO what is the correct condition?
                if (ir->ePBC == epbcXYZ)
                {
                    nonbondedForce->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
                    sys->setPeriodicBoxVectors(Vec3(state->box[0][0], 0, 0),
                                               Vec3(0, state->box[1][1], 0), Vec3(0, 0, state->box[2][2]));
                }
                else if (ir->ePBC == epbcNONE)
                    nonbondedForce->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
                else
                    gmx_fatal(FARGS,"OpenMM supports only full periodic boundary conditions "
                              "(pbc = xyz), or none (pbc = no).\n");
                nonbondedForce->setCutoffDistance(ir->rcoulomb);
                break;

            case eelEWALD:
                if (ir->ewald_geometry == eewg3DC)
                    gmx_fatal(FARGS,"OpenMM supports only Ewald 3D geometry.\n");
                if (ir->epsilon_surface != 0)
                    gmx_fatal(FARGS,"OpenMM does not support dipole correction in Ewald summation.\n");
                nonbondedForce->setNonbondedMethod(NonbondedForce::Ewald);
                nonbondedForce->setCutoffDistance(ir->rcoulomb);
                sys->setPeriodicBoxVectors(Vec3(state->box[0][0], 0, 0),
                                           Vec3(0, state->box[1][1], 0), Vec3(0, 0, state->box[2][2]));
                break;

            case eelPME:
                nonbondedForce->setNonbondedMethod(NonbondedForce::PME);
                nonbondedForce->setCutoffDistance(ir->rcoulomb);
                sys->setPeriodicBoxVectors(Vec3(state->box[0][0], 0, 0),
                                           Vec3(0, state->box[1][1], 0), Vec3(0, 0, state->box[2][2]));
                break;

            default:
                gmx_fatal(FARGS,"Internal error: you should not see this message, it that the"
                          "electrosatics option check failed. Please report this error!");
            }
        }

        for (int i = 0; i < numAtoms; ++i)
        {
            real c6 = nbfp[types[i]*2*ntypes+types[i]*2];
            real c12 = nbfp[types[i]*2*ntypes+types[i]*2+1];
            if (c12 <= 0)
                nonbondedForce->addParticle(charges[i], 1.0, 0.0);
            else
                nonbondedForce->addParticle(charges[i], pow(c12/c6, (real) (1.0/6.0)), c6*c6/(4.0*c12));
            sys->addParticle(masses[i]);
        }

        // Build a table of all exclusions.

        vector<set<int> > exclusions(numAtoms);
        for (int i = 0; i < numAtoms; i++)
        {
            int start = top->excls.index[i];
            int end = top->excls.index[i+1];
            for (int j = start; j < end; j++)
                exclusions[i].insert(top->excls.a[j]);
        }

        // Record the 1-4 interactions, and remove them from the list of exclusions.

        const int* nb14Atoms = (int*) idef.il[F_LJ14].iatoms;
        offset = 0;
        for (int i = 0; i < num14; ++i)
        {
            int type = nb14Atoms[offset++];
            int atom1 = nb14Atoms[offset++];
            int atom2 = nb14Atoms[offset++];
            real c6 = idef.iparams[type].lj14.c6A;
            real c12 = idef.iparams[type].lj14.c12A;
            real sigma, epsilon;
            if (c12 <= 0)
            {
                epsilon = (real) 0.0;
                sigma = (real) 1.0;
            }
            else
            {
                epsilon = (real) ((c6*c6)/(4.0*c12));
                sigma = (real) pow(c12/c6, (real) (1.0/6.0));
            }
            nonbondedForce->addException(atom1, atom2,
                                         fr->fudgeQQ*charges[atom1]*charges[atom2], sigma, epsilon);
            exclusions[atom1].erase(atom2);
            exclusions[atom2].erase(atom1);
        }

        // Record exclusions.

        for (int i = 0; i < numAtoms; i++)
        {
            for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            {
                if (i < *iter)
                    nonbondedForce->addException(i, *iter, 0.0, 1.0, 0.0);
            }
        }

        // Add GBSA if needed.
        t_atoms        atoms;
        atoms    = gmx_mtop_global_atoms(top_global);

        if (ir->implicit_solvent == eisGBSA)
        {
            GBSAOBCForce* gbsa = new GBSAOBCForce();
            sys->addForce(gbsa);
            gbsa->setSoluteDielectric(ir->epsilon_r);
            gbsa->setSolventDielectric(ir->gb_epsilon_solvent);
            gbsa->setCutoffDistance(nonbondedForce->getCutoffDistance());
            if (nonbondedForce->getNonbondedMethod() == NonbondedForce::NoCutoff)
                gbsa->setNonbondedMethod(GBSAOBCForce::NoCutoff);
            else if (nonbondedForce->getNonbondedMethod() == NonbondedForce::CutoffNonPeriodic)
                gbsa->setNonbondedMethod(GBSAOBCForce::CutoffNonPeriodic);
            else if (nonbondedForce->getNonbondedMethod() == NonbondedForce::CutoffPeriodic)
                gbsa->setNonbondedMethod(GBSAOBCForce::CutoffPeriodic);
            else
                gmx_fatal(FARGS,"OpenMM supports only Reaction-Field electrostatics with OBC/GBSA.\n");

            for (int i = 0; i < numAtoms; ++i)
                gbsa->addParticle(charges[i],
                                  top_global->atomtypes.gb_radius[atoms.atom[i].type],
                                  top_global->atomtypes.S_hct[atoms.atom[i].type]);
        }

        // Set constraints.

        const int* constraintAtoms = (int*) idef.il[F_CONSTR].iatoms;
        offset = 0;
        for (int i = 0; i < numConstraints; ++i)
        {
            int type = constraintAtoms[offset++];
            int atom1 = constraintAtoms[offset++];
            int atom2 = constraintAtoms[offset++];
            sys->addConstraint(atom1, atom2, idef.iparams[type].constr.dA);
        }
        const int* settleAtoms = (int*) idef.il[F_SETTLE].iatoms;
        offset = 0;
        for (int i = 0; i < numSettle; ++i)
        {
            int type = settleAtoms[offset++];
            int oxygen = settleAtoms[offset++];
            sys->addConstraint(oxygen, oxygen+1, idef.iparams[type].settle.doh);
            sys->addConstraint(oxygen, oxygen+2, idef.iparams[type].settle.doh);
            sys->addConstraint(oxygen+1, oxygen+2, idef.iparams[type].settle.dhh);
        }

        // Create an integrator for simulating the system.

        real friction = (ir->opts.tau_t[0] == 0.0 ? 0.0 : 1.0/ir->opts.tau_t[0]);
        Integrator* integ;
        if (ir->eI == eiMD || ir->eI == eiVV || ir->eI == eiVVAK)
        {
            integ = new VerletIntegrator(ir->delta_t);
            if ( ir->etc != etcNO)
            {
                real collisionFreq = ir->opts.tau_t[0] / 1000; /* tau_t (ps) / 1000 = collisionFreq (fs^-1) */
                AndersenThermostat* thermostat = new AndersenThermostat(ir->opts.ref_t[0], friction); 
                sys->addForce(thermostat);
            }
        }
        else if (ir->eI == eiBD)
        {
            integ = new BrownianIntegrator(ir->opts.ref_t[0], friction, ir->delta_t);
            static_cast<BrownianIntegrator*>(integ)->setRandomNumberSeed(ir->ld_seed); 
        }
        else if (EI_SD(ir->eI))
        {
            integ = new LangevinIntegrator(ir->opts.ref_t[0], friction, ir->delta_t);
            static_cast<LangevinIntegrator*>(integ)->setRandomNumberSeed(ir->ld_seed); 
        }
        integ->setConstraintTolerance(ir->shake_tol);

        // Create a context and initialize it.
        Context* context = NULL;
        if (platformOptStr == NULL || platformOptStr == "")
        {
            context = new Context(*sys, *integ);
        }
        else
        {
            // Find which platform is it.
            for (int i = 0; i < Platform::getNumPlatforms() && context == NULL; i++)
            {
                if (isStringEqNCase(opt->getOptionValue("platform"), Platform::getPlatform(i).getName()))
                {
                    Platform& platform = Platform::getPlatform(i);
                    // set standard properties
                    platform.setPropertyDefaultValue("CudaDevice", opt->getOptionValue("deviceid"));
                    // TODO add extra properties
                    context = new Context(*sys, *integ, platform);
                }
            }
            if (context == NULL)
            {
                gmx_fatal(FARGS, "The requested platform \"%s\" could not be found.\n", 
                        opt->getOptionValue("platform").c_str());
            }
        }

        Platform& platform = context->getPlatform();
        fprintf(fplog, "Gromacs will use the OpenMM platform: %s\n", platform.getName().c_str());

        const vector<string>& properties = platform.getPropertyNames();
        if (debug)
        {
            for (int i = 0; i < properties.size(); i++)
            {
                printf(">> %s: %s\n", properties[i].c_str(), 
                        platform.getPropertyValue(*context, properties[i]).c_str());
                fprintf(fplog, ">> %s: %s\n", properties[i].c_str(), 
                        platform.getPropertyValue(*context, properties[i]).c_str());
            }
        }

        int devId;
        if (!from_string<int>(devId, platform.getPropertyValue(*context, "CudaDevice"), std::dec))
        {
            gmx_fatal(FARGS, "Internal error: couldn't determine the device selected by OpenMM");
        }

        /* check GPU compatibility */
        char gpuname[STRLEN];
        if (!is_supported_cuda_gpu(-1, gpuname))
        {
            if (!gmx_strcasecmp(opt->getOptionValue("force-device").c_str(), "yes"))
            {
                sprintf(warn_buf, "Non-supported GPU selected (#%d, %s), forced continuing.\n"
                        "Note, that the simulation can be slow or it migth even crash.", 
                        devId, gpuname);
                fprintf(fplog, "%s", warn_buf);
                gmx_warning(warn_buf);
            }
            else
            {
                gmx_fatal(FARGS, "The selected GPU (#%d, %s) is not supported by Gromacs! "
                          "Most probably you have a low-end GPU which would not perform well, " 
                          "or new hardware that has not been tested yet with Gromacs-OpenMM. "
                          "If you still want to try using the device, use the force=on option.", 
                          devId, gpuname);
            }
        }
        else
        {
            fprintf(fplog, "Gromacs will run on the GPU #%d (%s).\n", devId, gpuname);
        }

        /* do the pre-simulation memtest */
        run_memtest(fplog, -1, "Pre", opt);

        vector<Vec3> pos(numAtoms);
        vector<Vec3> vel(numAtoms);
        for (int i = 0; i < numAtoms; ++i)
        {
            pos[i] = Vec3(state->x[i][0], state->x[i][1], state->x[i][2]);
            vel[i] = Vec3(state->v[i][0], state->v[i][1], state->v[i][2]);
        }
        context->setPositions(pos);
        context->setVelocities(vel);

        // Return a structure containing the system, integrator, and context.
        OpenMMData* data = new OpenMMData();
        data->system = sys;
        data->integrator = integ;
        data->context = context;
        data->removeCM = (ir->nstcomm > 0);
        data->platformOpt = opt;
        return data;

    }
    catch (std::exception& e)
    {
        gmx_fatal(FARGS, "OpenMM exception caught while initializating: %s\n", e.what());
    }
}

/**
 * Integrate one step.
 *
 * @param data the OpenMMData object created by openmm_init().
 */
void openmm_take_one_step(void* data)
{
    // static int step = 0; printf("----> taking step #%d\n", step++);
    try
    {
        static_cast<OpenMMData*>(data)->integrator->step(1);
    }
    catch (std::exception& e)
    {
        gmx_fatal(FARGS, "OpenMM exception caught while taking a step: %s\n", e.what());
    }
}

/**
 * Integrate n steps.
 *
 * @param data the OpenMMData object created by openmm_init().
 */
void openmm_take_steps(void* data, int nstep)
{
    try
    {
        static_cast<OpenMMData*>(data)->integrator->step(nstep);
    }
    catch (std::exception& e)
    {
        gmx_fatal(FARGS, "OpenMM exception caught while taking a step: %s\n", e.what());
    }
}

/**
 * Clean up all the data structures used by OpenMM.
 *
 * @param log file pointer
 * @param data the OpenMMData object created by openmm_init().
 */
void openmm_cleanup(FILE* fplog, void* data)
{
    OpenMMData* d = static_cast<OpenMMData*>(data);
    run_memtest(fplog, -1, "Post", d->platformOpt);
    delete d->system;
    delete d->integrator;
    delete d->context;
    delete d->platformOpt;
    delete d;
}

/**
 * Copy the current state information (positions, velocities, and forces) into the data structures used
 * by Gromacs.
 *
 * @param data the OpenMMData object created by openmm_init().
 */
void openmm_copy_state(void *data,
                       t_state *state, double *time,
                       rvec f[], gmx_enerdata_t *enerd,
                       bool includePos, bool includeVel, bool includeForce, bool includeEnergy)
{
    int types = 0;
    if (includePos)
        types += State::Positions;
    if (includeVel)
        types += State::Velocities;
    if (includeForce)
        types += State::Forces;
    if (includeEnergy)
        types += State::Energy;
    if (types == 0)
        return;
    try
    {
        State currentState = static_cast<OpenMMData*>(data)->context->getState(types);
        int numAtoms =  static_cast<OpenMMData*>(data)->system->getNumParticles();
        if (includePos)
        {
            for (int i = 0; i < numAtoms; i++)
            {
                Vec3 x = currentState.getPositions()[i];
                state->x[i][0] = x[0];
                state->x[i][1] = x[1];
                state->x[i][2] = x[2];
            }
        }
        if (includeVel)
        {
            for (int i = 0; i < numAtoms; i++)
            {
                Vec3 v = currentState.getVelocities()[i];
                state->v[i][0] = v[0];
                state->v[i][1] = v[1];
                state->v[i][2] = v[2];
            }
        }
        if (includeForce)
        {
            for (int i = 0; i < numAtoms; i++)
            {
                Vec3 force = currentState.getForces()[i];
                f[i][0] = force[0];
                f[i][1] = force[1];
                f[i][2] = force[2];
            }
        }
        if (includeEnergy)
        {
            int numConstraints = static_cast<OpenMMData*>(data)->system->getNumConstraints();
            int dof = 3*numAtoms-numConstraints;
            if (static_cast<OpenMMData*>(data)->removeCM)
                dof -= 3;
            enerd->term[F_EPOT] = currentState.getPotentialEnergy();
            enerd->term[F_EKIN] = currentState.getKineticEnergy();
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];
            enerd->term[F_TEMP] = 2.0*enerd->term[F_EKIN]/dof/BOLTZ;
        }
        *time = currentState.getTime();
    }
    catch (std::exception& e)
    {
        gmx_fatal(FARGS, "OpenMM exception caught while retrieving state information: %s\n", e.what());
    }
}
