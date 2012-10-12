#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "gmx_babelio.h"

// Include Open Babel classes for OBMol and OBConversion
#ifdef HAVE_LIBOPENBABEL2
#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;

int read_babel(int argc,char **argv)
{
  // Read from STDIN (cin) and Write to STDOUT (cout)
  OpenBabel::OBConversion conv(&cin,&cout);
  
  // Try to set input format to MDL SD file
  // and output to SMILES
  if (conv.SetInAndOutFormats("SDF","SMI"))
    {
      OpenBabel::OBMol mol;
      if (conv.Read(&mol))
        {
          //  ...manipulate molecule
          cerr << " Molecule has: " << mol.NumAtoms()
               << " atoms." << endl;
        }
      
      // Write SMILES to the standard output
      conv.Write(&mol);
    }
  return 0; // exit with success
}

#else

#include "gmx_fatal.h"

int read_babel(int argc,char **argv)
{
  gmx_fatal(FARGS,"No support for OpenBabel");
  return -1;
}
#endif
