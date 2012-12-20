class MolecularComposition {
  const char *GetAtom(int *cnumber) 
  {
    std::string str = GetAtom(cnumber);
    
    return str.c_str();
  }

  const std::string GetAtom(int *cnumber) 
  {
    if (_catom_index < _catom.size())
      { 
        *cnumber = _cnumber[_catom_index]; 
        return _catom[_catom_index++]; 
      } 
    else
      {
        Reset();
        *cnumber = 0;
        return NULL;
      }
  }
  
  void AddAtom(std::string catom,const int cnumber)
  {
    int i;
    
    for(i=0; (i<_catom.size()); i++) 
      {
        if (_catom[i] == catom) 
          {
            _cnumber[i] += cnumber;
            break;
          }
      }
    if (i == _catom.size()) 
      {
        _catom.push_back(catom);
        _cnumber.push_back(cnumber);
      }
  
  }
  
  void AddAtom(const char *catom,const int cnumber)
  {
    std::string str(catom);
    
    AddAtom(str,cnumber);
  }
  
  void DeleteAtom(std::string catom)
  {
    int i;
    
    for(i=0; (i<_catom.size()); i++)
      {
        if (_catom[i].c_str() == catom)
          {
            _catom.erase(_catom.begin()+i-1);
            _cnumber.erase(_cnumber.begin()+i-1);
          }
      }
  }
  
  void DeleteAtom(const char *catom)
  {
    std::string str(catom);
    
    DeleteAtom(str);
  }
  
  void ReplaceAtom(std::string oldatom,std::string newatom)
  {
    int i;
    
    for(i=0; (i<_catom.size); i++) 
      {
        if (oldatom == catom[i])
          {
            catom[i].assign(newatom);
            break;
          }
      }
  }
  
  void ReplaceAtom(const char *oldatom,const char *newatom)
  {
    std::string oa(oldatom);
    std::string na(newatom);
    
    ReplaceAtom(oa,na);
  }
  
  int CountAtoms(std::string atom)
  {
    int i;
    
    for(i=0; (i<_catom.size()); i++)
      if (_catom[i] == atom)
        return _cnumber[i];
    return 0;
  }
  
  int CountAtoms(const char *atom)
  {
    std::string str(atom);
    
    return CountAtoms(str);
  }
}

class MolProp {
  const char *GetCategory() 
  { 
    if (_category_index < _category.size()) 
      { 
        return _category[_category_index++].c_str(); 
      }
    else {
      ResetCategory();
      return NULL; 
    }
  }
  std::string GetCategory() 
  { 
    if (_category_index < _category.size()) 
      { 
        return _category[_category_index++]; 
      }
    else {
      ResetCategory();
      return NULL; 
    }
  }

  int SearchCategory(const char *catname)
  {
    std::string str(catname);
    
    return SearchCategory(str);
  }
  
  int SearchCategory(std::string catname)
  {
    int i;
    
    for(i=0; (i < _category.size()); i++)
      if (_category[i] == catname)
        return 1;
    
    return 0;
  }
  
}
