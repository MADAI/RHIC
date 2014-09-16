/*  hydroDef.h
 *  Created by Joshua Vredevoogd on 2/16/09.
 */

// *********************** PLEASE NOTE ************************* //
// If anything changes in this file, the code must be recompiled //
// *********************** PLEASE NOTE ************************* //

// constants
#define ROOT2 1.41421356237309504880168872420969807856967187537694807317667973799
#define ROOT3 1.73205080756887729352744634150587236694280525381038062805580
#define JPI 3.14159265358979323846264338327950288419716939937510
#define JFOURPI (4.*JPI)
#define JIFPI (1./JFOURPI)
#define JHBARC (.197326968) // in GeV*fm

// scaling factor for the israel-stewart terms must be:
// 'c' for constant
// 'e' for energy density
// 's' for entropy density
// 'j' for (correct) standard (e+p)
#define ISSCALE 'c' 
		
	//#define FLAT_EOS

// instead of the rescale it adds the necessary divergence of u
// do not use unless above = c
#define ISRESCALE true

// maximum grid sizes in each direction
#define XSIZE 200		// number of mesh points in radial directions
#define YSIZE 200		// number of mesh points in radial directions
#define NSIZE 100		// number of mesh points in eta direction

	//#define HYDRO_XDMF_OUTPUT
	
	//#define HYDRO_BOOST_THREADS
#define NTHREADS 4
