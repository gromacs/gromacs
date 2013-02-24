#include <string>
#include <vector>
#include "molprop.hpp"
        
namespace alexandria {
  void MolecularComposition::ReplaceAtom(std::string oldatom,std::string newatom)
  {
    AtomNumIterator i;
    
    for(i=BeginAtomNum(); (i<=EndAtomNum()); i++) 
      {
        if (oldatom == i->GetAtom())
          {
            i->SetAtom(newatom);
            break;
          }
      }
  }
  
  int MolecularComposition::CountAtoms(std::string atom)
  {
    AtomNumIterator i;
    
    for(i=BeginAtomNum(); (i<=EndAtomNum()); i++) 
      {
        if (atom == i->GetAtom())
          return i->GetNumber();
      }
    return 0;
  }
  
  int MolecularComposition::CountAtoms(const char *atom)
  {
    std::string str(atom);
    
    return CountAtoms(str);
  }
  
  int MolecularComposition::CountAtoms()
  {
    int nat = 0;
    AtomNumIterator i;
    
    for(i=BeginAtomNum(); (i<=EndAtomNum()); i++) 
      {
        nat += i->GetNumber();
      }
    return nat;
  }

  void MolProp::CheckConsistency() 
  {
  }
  
  int MolProp::SearchCategory(std::string catname)
  {
    std::vector<std::string>::iterator i;
    
    for(i=BeginCategory(); (i < EndCategory()); i++)
      {
        if (*i == catname)
          return 1;
      }
    return 0;
  }
  
  void MolProp::DeleteComposition(std::string compname)
  {
    MolecularCompositionIterator i;
    
    for(i=BeginMolecularComposition(); (i<=EndMolecularComposition()); i++) 
      {
        if (i->CompName() == compname) 
          {
            break;
          }
      }
    if (i < EndMolecularComposition())
      _mol_comp.erase(i);
  }
}

