#include <stdio.h>
#include "typedefs.h"
#include "types/simple.h"
#include "vec.h"
#include "pbc.h"
#include "localpressure.h"

#define SPACINGFACTOR 4.0


/*
 Basic idea:
 Take grid, line between two point sources, virial tensor
 Distance between point sources = k
 Divide into n point sources, where n/k < grid spacing
 For now, choose n/k << grid spacing, simply sum fractional contribution onto grid
 
 A more advanced approach would be to choose n/k ~ grid spacing and have each point
 source contribute a Gaussian blur onto grid, perhaps by FFT, Gaussian transform, FFT-1
 
 */

/*private functions: */
void 
grid_lookup(gmx_localp_grid_t *grid,rvec pt, int *i, int *j, int *k)
{
    int tx,ty,tz,nx,ny,nz;
    real rxx,ryx,ryy,rzx,rzy,rzz;
    
	/*Looks up indices for grid:
     fractional indices are invbox * coordinate; 
     grid indices are then nx*f_ind[XX], etc.
     */
    rxx = grid->invbox[XX][XX];
    ryx = grid->invbox[YY][XX];
    ryy = grid->invbox[YY][YY];
    rzx = grid->invbox[ZZ][XX];
    rzy = grid->invbox[ZZ][YY];
    rzz = grid->invbox[ZZ][ZZ];
    
    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;
    
    tx = nx * ( 2.0 + pt[XX] * rxx + pt[YY] * ryx + pt[ZZ] * rzx );
    ty = ny * ( 2.0 +                pt[YY] * ryy + pt[ZZ] * rzy );
    tz = nz * ( 2.0 +                               pt[ZZ] * rzz );
    
    *i = tx%nx;
    *j = ty%ny;
    *k = tz%nz;
}

/*Looks up indices for grid*/


void
gmx_spread_local_virial_on_grid(gmx_localp_grid_t *    grid,
                                real                   ix, /*i is shifted*/
                                real                   iy,
                                real                   iz,
                                real                   jx,
                                real                   jy,
                                real                   jz,
                                real                   fx,
                                real                   fy,
                                real                   fz)
{
	/*takes reals for points i and j, delta-F between points.
     Calculates local virial and spreads onto grid
     */
    
	rvec linespread,curpt,coord1,coord2;
	int n,ctr,i,j,k;
    real k_line;
    /*	Use following two variables for 
     real grid_spacing;
     rvec  gridsize;
     */
	matrix gridres;
	matrix *pgrid;
    
	/*
     If need to calculate grid spacing on the fly, use the following:
     Think unit cell in 3D is [1/nx 1/ny 1/nz] * box
     Then take (min(coord[XX],coord[YY]),coord[ZZ])
     curpt[XX]=1.0/grid->nz;
     curpt[YY]=1.0/grid->ny;
     curpt[ZZ]=1.0/grid->nz;
     mvmul(grid->box,curpt,gridsize);
     grid_spacing=min(min(gridsize[XX],gridsize[YY]),gridsize[ZZ]);
     */
	
	/*virial tensor is Rij x delta-F (delta-F is passed to this function)*/
	matrix virial;

    if(grid==NULL)
        return;

    pgrid = (grid->bLR==TRUE) ? grid->longrange_grid : grid->current_grid;
        
	/*virial tensor is Rij x delta-F (delta-F is passed to this function)*/
	coord1[XX]=ix-jx;
	coord1[YY]=iy-jy;
	coord1[ZZ]=iz-jz;
    
	virial[XX][XX]=-0.5*coord1[XX]*fx;
	virial[XX][YY]=-0.5*coord1[XX]*fy;
	virial[XX][ZZ]=-0.5*coord1[XX]*fz;
	
	virial[YY][XX]=-0.5*coord1[YY]*fx;
	virial[YY][YY]=-0.5*coord1[YY]*fy;
	virial[YY][ZZ]=-0.5*coord1[YY]*fz;
	
	virial[ZZ][XX]=-0.5*coord1[ZZ]*fx;
	virial[ZZ][YY]=-0.5*coord1[ZZ]*fy;
	virial[ZZ][ZZ]=-0.5*coord1[ZZ]*fz;
	
	/*map input reals into coord1,coord2*/
	coord1[XX]=ix;
	coord1[YY]=iy;
	coord1[ZZ]=iz;
	coord2[XX]=jx;
	coord2[YY]=jy;
	coord2[ZZ]=jz;
	
	rvec_sub(coord2,coord1,linespread); /*linespread = coord2 - coord1 in PBC reflection s.t. distance is minimized*/
    
	k_line=norm(linespread); /*k=length of linespread)*/
	/* n/k = grid->spacing/SPACINGFACTOR */
	n= (SPACINGFACTOR * k_line / grid->spacing);
    
	msmul(virial,1.0/(1.0+n),gridres);

	for (ctr=0;ctr<=n;++ctr) {
		/*curpt = coord1 + linespread/n*/
		curpt[XX]=coord1[XX]+linespread[XX]/n*ctr;
		curpt[YY]=coord1[YY]+linespread[YY]/n*ctr;
		curpt[ZZ]=coord1[ZZ]+linespread[ZZ]/n*ctr;
		//put_atom_in_box(grid->box,curpt);
		/*Get indices for grid*/
		grid_lookup(grid,curpt,&i,&j,&k);
		/*distribute pressure--divide by n+1 because doing range 0..n*/
		/*clean way would be to have tmp variable gridres, do following*/
		m_sub(pgrid[i*grid->nz*grid->ny+j*grid->nz+k],gridres,pgrid[i*grid->nz*grid->ny+j*grid->nz+k]);
	}    
}

void
gmx_spread_local_virial_on_grid_mat(gmx_localp_grid_t *    grid,
                                    real                   ix, /*i is shifted*/
                                    real                   iy,
                                    real                   iz,
                                    matrix				   virial)
{
	/*This version is for single-point assignment given a pressure matrix
     */
    
	rvec curpt;
	int i,j,k;
    /*	Use following two variables for 
     real grid_spacing;
     rvec  gridsize;
     */
	matrix gridres;
	matrix *pgrid;

    if(grid==NULL)
        return;

    pgrid = (grid->bLR==TRUE) ? grid->longrange_grid : grid->current_grid;

	/*
     If need to calculate grid spacing on the fly, use the following:
     Think unit cell in 3D is [1/nx 1/ny 1/nz] * box
     Then take (min(coord[XX],coord[YY]),coord[ZZ])
     curpt[XX]=1.0/grid->nz;
     curpt[YY]=1.0/grid->ny;
     curpt[ZZ]=1.0/grid->nz;
     mvmul(grid->box,curpt,gridsize);
     grid_spacing=min(min(gridsize[XX],gridsize[YY]),gridsize[ZZ]);
     */
	/*virial tensor is Rij x delta-F (delta-F is passed to this function)*/
        
	/*map input reals into curpt*/
	curpt[XX]=ix;
	curpt[YY]=iy;
	curpt[ZZ]=iz;
	
	
	//put_atom_in_box(grid->box,curpt);
	/*Get indices for grid*/
	grid_lookup(grid,curpt,&i,&j,&k);
	m_sub(pgrid[i*grid->nz*grid->ny+j*grid->nz+k],virial,pgrid[i*grid->nz*grid->ny+j*grid->nz+k]);
    
}

