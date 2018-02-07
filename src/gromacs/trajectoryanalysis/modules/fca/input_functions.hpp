#ifndef INPUT_FUNCTIONS_HPP
#define INPUT_FUNCTIONS_HPP

#include <vector>
#include <cstdio>
#include <cassert>

#include "utils.hpp"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/cstringutil.h"

namespace FCA{

namespace input_functions{

/////////////////////////////////////////////////////////////////////////////////////
/// File loading functions
/////////////////////////////////////////////////////////////////////////////////////

inline int read_conffile(const char* confin, char* title, rvec* x[]) {
    /* read coordinates out of STX file  */
    int natoms;
    t_atoms confat;
    matrix box;

    printf("read coordnumber from file %s\n", confin);
    get_stx_coordnum(confin, &natoms);
    printf("number of coordinates in file %d\n", natoms);
    /*  if (natoms != ncoords)
     fatal_error(0,"number of coordinates in coordinate file (%s, %d)\n"
     "             does not match topology (= %d)",
     confin,natoms,ncoords);
     else {*/
    /* make space for coordinates and velocities */
    init_t_atoms(&confat, natoms, FALSE);
    printf("init_t\n");
    snew(*x, natoms);
    read_stx_conf(confin, title, &confat, *x, nullptr, nullptr /*epbcXYZ*/, box);
    return natoms;
}

inline void filter2x(rvec* x, int nindex, int index[], int ngro, int igro[],
              rvec* xori, const char* structure) {
    /* filter2edx copies coordinates from x to edx which are given in index
     */
    int pos, i;
    for(i = 0; i < nindex; i++) {
        for(pos = 0; pos < ngro - 1 && igro[pos] != index[i]; ++pos) {
        } /*search element in igro*/
        if(igro[pos] != index[i]) {
            fprintf(stderr, "Couldn't find atom with index %d in structure %s",
                    index[i], structure);
            exit(1);
        }
        copy_rvec(xori[pos], x[i]);
    }
}

inline bool get_structure(t_atoms* atoms, const char* IndexFile,
                   const char* StructureFile, rvec* x, int natoms, int index[]) {
    int* igro; /*index corresponding to target or origin structure*/
    int ngro;
    rvec* xtar;
    char title[STRLEN];
    char* grpname;

    const int ntar = read_conffile(StructureFile, title, &xtar);
    get_index(atoms, IndexFile, 1, &ngro, &igro, &grpname);
    filter2x(x, natoms, index, ngro, igro, xtar, StructureFile);
    if(ngro != ntar){
       return false;
    }
    return true;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Reading frames
/////////////////////////////////////////////////////////////////////////////////////

inline int count_number_columns(FILE* fp) {
    constexpr size_t BUFFSIZE = 20000;
    char buffer[BUFFSIZE];
    const char* ptrret = fgets(buffer, BUFFSIZE, fp);
    assert(ptrret != nullptr);

    // skip first spaces
    int idxRead = 0;
    while(isspace(buffer[idxRead]) && buffer[idxRead]) {
        ++idxRead;
    } //load_frame_step spaces

    int cptColumns   = 0;
    while(buffer[idxRead]) {
        // if it is a char count it as a column
        if(!isspace(buffer[idxRead])){
            ++cptColumns;
        }
        // skip following char
        while(!isspace(buffer[idxRead]) && buffer[idxRead]) {
            ++idxRead;
        } //load_frame_step rest of number
        // skip following space
        while(isspace(buffer[idxRead]) && buffer[idxRead]) {
            ++idxRead;
        } //load_frame_step spaces
    }
    fprintf(stderr, "analyze line %s\n", buffer);
    fprintf(stderr, "found %d cols\n", cptColumns);
    return cptColumns;
}

inline std::vector< std::unique_ptr< real[] > > read_fca_proj(const char pos_file[], int* dim) {
    FILE* fp   = gmx_ffopen(pos_file, "r");
    (*dim) = count_number_columns(fp);
    rewind(fp);
    fprintf(stderr, "read vectors of dimensionality %d\n", *dim);

    std::vector< std::unique_ptr< real[] > > projx;
    while(true){
        std::unique_ptr< real[] > nextFrame(new real[*dim]);

        int cptRead = 0;
        while(cptRead < (*dim) && fscanf(fp, "%f", &(nextFrame[cptRead])) > 0){
            cptRead += 1;
        }

        if(cptRead != (*dim)){
            break;
        }

        projx.emplace_back(std::move(nextFrame));
    }
    gmx_ffclose(fp);
    return projx;
}


}
}

#endif
