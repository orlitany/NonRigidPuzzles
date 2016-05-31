/*
 * Emanuele Rodola <rodola@in.tum.de>
 * Mar 2014
 *
 * to compile:
 * mex -v -O conn_comp.cpp ../cvlab/cvlab/math3d.cpp ../cvlab/cvlab/mesh.cpp ../cvlab/cvlab/cloud3d.cpp ../cvlab/cvlab/spin_image.cpp ../cvlab/cvlab/string_utilities.cpp -I../cvlab/ -I../boost_1_54_0
 */

#include "mex.h"
#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include <cvlab/mesh.h>
#include <boost/cstdint.hpp>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 3 || nlhs != 3)
      mexErrMsgTxt("Usage: [components, nc, outliers] = conn_comp(vertices, triangles, mask).");

	const double* const pts = mxGetPr(prhs[0]);
	const int np = int( mxGetN(prhs[0]) );
	const double* const tri = mxGetPr(prhs[1]);
	const int nt = int( mxGetN(prhs[1]) );
	const double* const mask = mxGetPr(prhs[2]);
	
	if (np == 3)
		mexErrMsgTxt("It seems like you only have 3 vertices. Please try to transpose the input matrix.");
		
	if (nt == 3)
		mexErrMsgTxt("It seems like you only have 3 triangles. Please try to transpose the input matrix.");
	
	//std::cout << np << " vertices, " << nt << " triangles" << std::endl;
	
	const int mm = int( mxGetM(prhs[2]) );
	const int mn = int( mxGetN(prhs[2]) );
	
	if (std::min(mm,mn) != 1)
		mexErrMsgTxt("The segment mask must be either a row or a column vector.");
	
	if (nt != std::max(mm,mn))
		mexErrMsgTxt("The mask should have the same number of elements as there are triangles.");

	// Load the mesh

	cvlab::mesh mesh;
	
	std::vector<cvlab::point3d> vertices(np);
	
	for (int i=0; i<np; ++i)
	{
		cvlab::point3d& pt = vertices[i];
		pt.x = pts[i*3];
		pt.y = pts[i*3+1];
		pt.z = pts[i*3+2];
		//std::cout << pt << std::endl;
	}
	
	mesh.put_vertices(vertices);
	
	for (int i=0; i<nt; ++i)
	{
		int a, b, c;
		a = tri[i*3] - 1; // 1-based to 0-based
		b = tri[i*3+1] - 1;
		c = tri[i*3+2] - 1;
		mesh.add_triangle(a,b,c);
		//std::cout << a << " " << b << " " << c << std::endl;
	}
	
	// Detect connected components for the mask given in input
	
	plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
	double* label = mxGetPr(plhs[0]);
	
	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	double* nc = mxGetPr(plhs[1]);
	
	plhs[2] = mxCreateDoubleMatrix(nt, 1, mxREAL);
	double* outliers = mxGetPr(plhs[2]);
	
	//std::cout << "Detecting connected components... " << std::flush;
	
	for (int i=0; i<nt; ++i)
	{
		label[i] = 0;
		outliers[i] = 0;
	}
	
	int cur_label = 0;
	
	// NOTE:
	// mask will be tested with 0.9 instead of true/false because it is of class double coming from matlab
	
	while(true)
	{
		++cur_label;
		
		// Find next connected component
		
		int seed = -1;
		for (int i=0; i<nt; ++i)
		{
			if (mask[i]>0.9 && label[i]==0 && outliers[i]==0)
			{
				seed=i;
				break;
			}
		}
		
		if (seed == -1)
			break;
		
		//std::cout << "starting new component from triangle: " << seed << std::endl;
		
		//
		
		std::queue<int> seeds;
		seeds.push(seed);
		
		while (!seeds.empty())
		{
			//std::cout << "queue size: " << seeds.size() << std::endl;
			
			const int s = seeds.front();
			seeds.pop();
			
			std::vector<int> neighs;
			mesh.get_tri_imm_neighbors(s, neighs);
		
			const int n_neighs = neighs.size();
		
			// outliers do not increase the final label count
			if (n_neighs==0)
			{
				//std::cout << "[WARNING] The current seed has no neighboring triangles in the mesh!" << std::endl;
				outliers[s] = 1;
				--cur_label;
				continue;
			}
			
			//std::cout << "new seed triangle: " << s << std::endl;
			label[s] = cur_label;
		
			for (int k=0; k<n_neighs; ++k)
			{
				const int cur_neigh = neighs[k];
				
				if (s==cur_neigh)
				{
					std::cout << "[WARNING] It seems like triangle seeds are included among their neighbors." << std::endl;
					continue;
				}
				
				if (mask[cur_neigh]>0.9 && label[cur_neigh] == 0)
				{
					//std::cout << "seed " << s << ", pushing " << cur_neigh << std::endl;
					seeds.push(cur_neigh);
					label[cur_neigh] = cur_label; // NOTE this is a workaround due to get_tri_imm_neighbors() containing repeated indices
				}
			} // next neighboring triangle
		} // next seed
	} // next connected component
	
	*nc = cur_label-1;
	//std::cout << "done." << std::endl;
}
