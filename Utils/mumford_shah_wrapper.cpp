/**
 * Emanuele Rodolà
 * TU Munich
 * Mar 2015
 *
 * mex mumford_shah_wrapper.cpp COMPFLAGS="/openmp $COMPFLAGS"
 */
#include "mex.h"
#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include "mesh.h"
#include <climits>

#include "class_handle.hpp"

#define ALMOST_ZERO 1e-10


void initialize(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 2 || nlhs > 1)
      mexErrMsgTxt("Usage: [mesh] = init_mesh(vertices, triangles).");

	const double* const pts = mxGetPr(prhs[0]);
	const int np = int( mxGetN(prhs[0]) );
	const double* const tri = mxGetPr(prhs[1]);
	const int nt = int( mxGetN(prhs[1]) );
	
	if (np == 3)
		mexErrMsgTxt("It seems like you only have 3 vertices. Please try to transpose the input matrix.");
		
	if (nt == 3)
		mexErrMsgTxt("It seems like you only have 3 triangles. Please try to transpose the input matrix.");
	
	// Load the mesh   
	mesh_t* mesh = new mesh_t();
	
	std::vector< vec3d<double> > vertices(np);
	
	for (int i=0; i<np; ++i)
	{
		vec3d<double>& pt = vertices[i];
		pt.x = pts[i*3];
		pt.y = pts[i*3+1];
		pt.z = pts[i*3+2];
		//std::cout << pt.x << " " << pt.y << " " << pt.z << std::endl;
	}
	
	mesh->put_vertices(vertices);
	
	for (int i=0; i<nt; ++i)
	{
		int a, b, c;
		a = (int)tri[i*3] - 1; // 1-based to 0-based
		b = (int)tri[i*3+1] - 1;
		c = (int)tri[i*3+2] - 1;
		mesh->add_triangle(a,b,c);
		//std::cout << a << " " << b << " " << c << std::endl;
	}
	
    int dims[2]; dims[0] = 1; dims[1] = 1;
    plhs[0] = convertPtr2Mat<mesh_t>(mesh);
	
    for (int j=0; j<nt; ++j)
	{
		mesh_t::triangle_data& tri = mesh->triangles[j];
		
		const vec3d<double>& xj1 = mesh->vertices[tri.p0].p;
		const vec3d<double>& xj2 = mesh->vertices[tri.p1].p;
		const vec3d<double>& xj3 = mesh->vertices[tri.p2].p;
		
		tri.E = dot_product(xj2-xj1, xj2-xj1);
		tri.F = dot_product(xj2-xj1, xj3-xj1);
		tri.G = dot_product(xj3-xj1, xj3-xj1);
    }
//     mexLock();
//     mexPrintf("Initialized mesh %d -\n",mesh);	
}

void compute_cost(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 6 || nlhs > 1)
      mexErrMsgTxt("Usage: [cost_per_tri] = calc_cost_mumford_shah(mesh_ptr, scalar_fun, sigma, area_weighted, mean, perturb).");

	const mesh_t* mesh = convertMat2Ptr<mesh_t>(prhs[0]);
    
    /*if(mesh->mesh_class_id != MESH_CLASS_ID){
        mexPrintf("Mesh id: %d\n",mesh->mesh_class_id);
        mexErrMsgTxt("COST: Wrong mesh pointer\n");
    }*/

	//std::cout << np << " vertices, " << nt << " triangles" << std::endl;
	
	const double* const u = mxGetPr(prhs[1]);
	const int nvert = int( mxGetM(prhs[1]) );
    
	if (int( mxGetN(prhs[1]) ) != 1 && int( mxGetM(prhs[1]) ) != 1)
		mexErrMsgTxt("The scalar function must be either a row or a column vector.");
	
	const double sigma = *mxGetPr(prhs[2]);
	const bool area_weighted = *mxGetPr(prhs[3]);
	const double mean = *mxGetPr(prhs[4]);
	
	const double* const perturb = mxGetPr(prhs[5]);
    
    int nt = mesh->triangles.size();
    int np = mesh->vertices.size();
	
    if (nvert!=np)
         mexErrMsgTxt("COST: wrong number of vertices");

	plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
	double* cost_ms = mxGetPr(plhs[0]);
	
	const double var = 2*sigma*sigma;
	
//     mexPrintf("Mesh id: %d\n",mesh->mesh_class_id);
     
	 double* h = new double[np];
	 
	#pragma omp parallel for
	for (int i=0; i<np; ++i)
		h[i] = std::exp( - (u[i]-mean)*(u[i]-mean) / var );
	 
    #pragma omp parallel for
	for (int j=0; j<nt; ++j)
	{
		const mesh_t::triangle_data& tri = mesh->triangles[j];
						
		const double E = tri.E;
		const double F = tri.F;
		const double G = tri.G;
		
		const double u1 = u[tri.p0];
		const double u2 = u[tri.p1];
		const double u3 = u[tri.p2];
		
		const double ua = u2 - u1;
		const double ub = u3 - u1;
		
		double norm_grad_u = std::sqrt(ua*ua*G + ub*ub*E - 2*ua*ub*F);
		
		if (area_weighted)
		{
			const double det = std::sqrt( fabs(E*G-F*F) );
			if (det <= ALMOST_ZERO)
				continue;
			norm_grad_u /= det;
		}
		
		const double h1 = perturb[tri.p0] * h[tri.p0];
		const double h2 = perturb[tri.p1] * h[tri.p1];
		const double h3 = perturb[tri.p2] * h[tri.p2];
		
		const double hsum = h1+h2+h3;
		
		if (hsum > ALMOST_ZERO)
			cost_ms[j] = hsum * norm_grad_u / 6;
	}
	
	delete[] h;
}

void compute_gradient(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 6 || nlhs > 1)
      mexErrMsgTxt("Usage: [grad_ms] = calc_grad_mumford_shah(mesh_ptr, scalar_fun, sigma, area_weighted, mean, perturb).");

	const mesh_t* mesh = convertMat2Ptr<mesh_t>(prhs[0]);
	
    //if(mesh->mesh_class_id != MESH_CLASS_ID)
    //     mexErrMsgTxt("GRADIENT: Wrong mesh pointer");

	//std::cout << np << " vertices, " << nt << " triangles" << std::endl;
	
	const double* const u = mxGetPr(prhs[1]);
	
	if (int( mxGetN(prhs[1]) ) != 1 && int( mxGetM(prhs[1]) ) != 1)
		mexErrMsgTxt("The scalar function must be either a row or a column vector.");
		
	const double sigma = *mxGetPr(prhs[2]);
	const bool area_weighted = *mxGetPr(prhs[3]);
	const double mean = *mxGetPr(prhs[4]);
	
	const double* const perturb = mxGetPr(prhs[5]);
    
    int nt = mesh->triangles.size();
    int np = mesh->vertices.size();
    
	plhs[0] = mxCreateDoubleMatrix(np, 1, mxREAL);
	double* grad_ms = mxGetPr(plhs[0]);
	
	const double var = sigma*sigma;
	
	double* h = new double[np];
	 
	#pragma omp parallel for
	for (int i=0; i<np; ++i)
		h[i] = std::exp( - (u[i]-mean)*(u[i]-mean) / (2*var) );
	
    #pragma omp parallel for
	for (int k=0; k<np; ++k)
	{
		const mesh_t::vertex_data& v = mesh->vertices[k];
		grad_ms[k] = 0.;
		
		for (std::list<int>::const_iterator it = v.tri_indices.begin(); it != v.tri_indices.end(); ++it)
		{
			const int j = *it;
			const mesh_t::triangle_data& tri = mesh->triangles[j];
			
			int p1, p2, p3;
			if (k==tri.p0)
			{
				p1 = tri.p0;
				p2 = tri.p1;
				p3 = tri.p2;
			}
			else if (k==tri.p1)
			{
				p1 = tri.p1;
				p2 = tri.p2;
				p3 = tri.p0;
			}
			else if (k==tri.p2)
			{
				p1 = tri.p2;
				p2 = tri.p0;
				p3 = tri.p1;
			}
			else
				mexErrMsgTxt("Unexpected condition: Cannot find vertex among its triangles.");
			
			
			const vec3d<double>& xj1 = mesh->vertices[p1].p;
			const vec3d<double>& xj2 = mesh->vertices[p2].p;
			const vec3d<double>& xj3 = mesh->vertices[p3].p;
			
			const double E = dot_product(xj2-xj1, xj2-xj1);
			const double F = dot_product(xj2-xj1, xj3-xj1);
			const double G = dot_product(xj3-xj1, xj3-xj1);
			
			const double uk = u[p1];
			const double u2 = u[p2];
			const double u3 = u[p3];
			
			const double ua = u2 - uk;
			const double ub = u3 - uk;
			
			const double Dj = std::sqrt( ua*ua*G + ub*ub*E - 2*ua*ub*F );
			
			if (Dj <= ALMOST_ZERO)
				continue;
			
			const double hk = perturb[p1] * h[p1];
			const double dh = -(uk-mean)*hk/var;
			
			const double h2 = perturb[p2] * h[p2];
			const double h3 = perturb[p3] * h[p3];
			
			const double weight = (area_weighted ? std::sqrt(fabs(E*G-F*F)) : 1.);
			
			grad_ms[k] += ( (hk+h2+h3)*( (-ua)*(G-F) + (-ub)*(E-F) )/Dj + Dj*dh ) / weight;
		}
		grad_ms[k] /= 6;
	}
	
	delete[] h;
}

void destroy(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 1 || nlhs > 0)
      mexErrMsgTxt("Usage: destroy_mesh(mesh_ptr).");
    
    destroyObject<mesh_t>(prhs[0]);
//     mexUnlock();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs < 1 || nlhs > 1)
      mexErrMsgTxt("Usage: init_mesh(action, ...) where actions={0:initialize, -1:destroy, 1:compute_cost, 2:compute_gradeint} .");
    
    if(mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
         mexErrMsgTxt("action should be a double value (double)");
    int action = *mxGetPr(prhs[0]);
    
    switch(action)
    {
        case 0: initialize( nlhs, plhs, nrhs-1, prhs+1); break;
        case -1: destroy( nlhs, plhs, nrhs-1, prhs+1); break;
        case 1: compute_cost( nlhs, plhs, nrhs-1, prhs+1); break;
        case 2: compute_gradient( nlhs, plhs, nrhs-1, prhs+1); break;
        default:
             mexErrMsgTxt("Usage: init_mesh(action, ...) where actions={0:initialize, -1:destroy, 1:compute_cost, 2:compute_gradeint} .");
    }
}
