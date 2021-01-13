#ifndef MULTEM_ATOM_DATA_API
#define MULTEM_ATOM_DATA_API

#include "math.cuh"
#include "safe_types.cuh"
#include "lin_alg_def.cuh"

namespace mt {

  #ifdef __CUDACC__
	template <class T, eDevice dev>
	struct Atom_Type;
  #endif

	template <class T>
	class Atom_Data
	{
		public:
			using value_type = T;
			using size_type = std::size_t;
			using TVector_r = std::vector<T>;

			Atom_Data(): l_x(0), l_y(0), l_z(0), dz(0.25), 
				ct_na(1), ct_nb(1), ct_nc(1), ct_a(0), 
				ct_b(0), ct_c(0), ct_x0(0), ct_y0(0), 
				Z_min(0), Z_max(0), x_min(0), x_max(0), 
				y_min(0), y_max(0), z_min(0), z_max(0), 
				sigma_min(0), sigma_max(0), occ_min(0), 
				occ_max(0), region_min(0), region_max(0), 
				R_int_min(0), R_int_max(0), 
				x_mean(0), y_mean(0), z_mean(0), x_std(0), 
				y_std(0), z_std(0), s_x(0), s_y(0), s_z(0), 
				x_int_min(0), x_int_max(0), y_int_min(0), 
				y_int_max(0), z_int_min(0), z_int_max(0), 
				s_x_int(0), s_y_int(0), s_z_int(0), 
				l_x_int(0), l_y_int(0), l_z_int(0){}

			size_type size() const
			{
				return Z.size();
			}

			bool empty() const
			{
				return size() == 0;
			}

			void clear()
			{
				l_x = l_y = l_z = 0; 
				dz = 0.25; 
				ct_na = ct_nb = ct_nc = 1;
				ct_a = ct_b = ct_c = 0;
				ct_x0 = ct_y0 = 0;
				Z_min = Z_max = 0;
				x_min = x_max = 0; 
				y_min = y_max = 0;
				z_min = z_max = 0;
				sigma_min = sigma_max = 0;
				occ_min = occ_max = 0; 
				region_min = region_max = 0;
				R_int_min = R_int_max = 0;
				x_mean = y_mean = 0;
				z_mean = x_std = 0;
				y_std = z_std = 0; 
				s_x = s_y = s_z = 0;
				x_int_min = x_int_max = 0;
				y_int_min = y_int_max = 0;
				z_int_min = z_int_max = 0;
				s_x_int = s_y_int = s_z_int = 0;
				l_x_int =l_y_int = l_z_int = 0;

				Z.clear();
				x.clear();
				y.clear();
				z.clear();
				sigma.clear();
				occ.clear();
				region.clear();
				charge.clear();

				amorp_lay_info.clear();
				Z_unique.clear();

			}

			// resize number of atoms
			void resize(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				Z.resize(new_size);
				x.resize(new_size);
				y.resize(new_size);
				z.resize(new_size);
				sigma.resize(new_size);
				occ.resize(new_size);
				region.resize(new_size);
				charge.resize(new_size);
			}

			template <class TAtom_Data>
			void assign(TAtom_Data &atoms)
			{
				this->l_x = atoms.l_x;
				this->l_y = atoms.l_y;
				this->l_z = atoms.l_z;
				this->dz = atoms.dz;

				this->ct_na = atoms.ct_na;
				this->ct_nb = atoms.ct_nb;
				this->ct_nc = atoms.ct_nc; 

				this->ct_a = atoms.ct_a;
				this->ct_b = atoms.ct_b;
				this->ct_c = atoms.ct_c;

				this->ct_x0 = atoms.ct_x0;
				this->ct_y0 = atoms.ct_y0;		

				this->ct_x0 = atoms.ct_x0;
				this->ct_y0 = atoms.ct_y0;	

				this->amorp_lay_info.resize(atoms.amorp_lay_info.size());
				for(auto ik=0; ik<atoms.amorp_lay_info.size(); ik++)
					this->amorp_lay_info[ik] = atoms.amorp_lay_info[ik];

				this->Z.assign(atoms.Z.begin(), atoms.Z.end());
				this->x.assign(atoms.x.begin(), atoms.x.end());
				this->y.assign(atoms.y.begin(), atoms.y.end());
				this->z.assign(atoms.z.begin(), atoms.z.end());
				this->sigma.assign(atoms.sigma.begin(), atoms.sigma.end());
				this->occ.assign(atoms.occ.begin(), atoms.occ.end());
				this->region.assign(atoms.region.begin(), atoms.region.end());
				this->charge.assign(atoms.charge.begin(), atoms.charge.end());

				this->Z_unique.assign(atoms.Z_unique.begin(), atoms.Z_unique.end());

				this->Z_min = atoms.Z_min;
				this->Z_max = atoms.Z_max;

				this->x_min = atoms.x_min;
				this->x_max = atoms.x_max;

				this->y_min = atoms.y_min;
				this->y_max = atoms.y_max;

				this->z_min = atoms.z_min;
				this->z_max = atoms.z_max;

				this->sigma_min = atoms.sigma_min;
				this->sigma_max = atoms.sigma_max;

				this->occ_min = atoms.occ_min;
				this->occ_max = atoms.occ_max;

				this->region_min = atoms.region_min;
				this->region_max = atoms.region_max;

				this->R_int_min = atoms.R_int_min;
				this->R_int_max = atoms.R_int_max;

				this->x_mean = atoms.x_mean;
				this->y_mean = atoms.y_mean;
				this->z_mean = atoms.z_mean;

				this->x_std = atoms.x_std;
				this->y_std = atoms.y_std;
				this->z_std = atoms.z_std;

				this->s_x = atoms.s_x;
				this->s_y = atoms.s_y;
				this->s_z = atoms.s_z;

				this->x_int_min = atoms.x_int_min;
				this->x_int_max = atoms.x_int_max;

				this->y_int_min = atoms.y_int_min;
				this->y_int_max = atoms.y_int_max;

				this->z_int_min = atoms.z_int_min;
				this->z_int_max = atoms.z_int_max;

				this->s_x_int = atoms.s_x_int;
				this->s_y_int = atoms.s_y_int;
				this->s_z_int = atoms.s_z_int;

				this->l_x_int = atoms.l_x_int;
				this->l_y_int = atoms.l_y_int;
				this->l_z_int = atoms.l_z_int;
			}

			template <class TAtom_Data> 
			Atom_Data<T>& operator=(TAtom_Data &atoms)
			{
				assign(atoms);
				return *this; 
			}

			void shrink_to_fit()
			{
				Z.shrink_to_fit();
				x.shrink_to_fit();
				y.shrink_to_fit();
				z.shrink_to_fit();
				sigma.shrink_to_fit();
				occ.shrink_to_fit();
				region.shrink_to_fit();
				charge.shrink_to_fit();			
			}

			// set xtl_build parameters
			void set_crystal_parameters(int ct_na_i = 1, int ct_nb_i = 1, int ct_nc_i = 1, 
			T ct_a_i = 0, T ct_b_i = 0, T ct_c_i = 0, T ct_x0_i = 0, T ct_y0_i = 0)
			{
				ct_na = max(1, ct_na_i);
				ct_nb = max(1, ct_nb_i);
				ct_nc = max(1, ct_nc_i); 
				ct_a = max(T(0), ct_a_i); 
				ct_b = max(T(0), ct_b_i); 
				ct_c = max(T(0), ct_c_i); 
				ct_x0 = max(T(0), ct_x0_i);
				ct_y0 = max(T(0), ct_y0_i);
			}
			
			// set amorphous parameters
			void set_amorphous_parameters(std::vector<Amorp_Lay_Info<T>> &amorp_lay_info_i)
			{
				amorp_lay_info = amorp_lay_info_i;
			}

			// set atoms
			void set_atoms(size_type nr_atoms_i, size_type nc_atoms_i, 
			double *atoms_i, T l_x_i = 0, T l_y_i = 0, T l_z_i = 0, T dz_i= 0.25)
			{
				resize(nr_atoms_i);

				l_x = l_x_i;
				l_y = l_y_i;
				l_z = l_z_i;
				dz = dz_i;

				size_type iatoms_c = 0;
				for(auto iatoms = 0; iatoms < nr_atoms_i; iatoms++)
				{
					auto atom = read_atom(nr_atoms_i, nc_atoms_i, atoms_i, iatoms);
					//if(atom.is_xy_positive())
					Z[iatoms_c] = atom.Z; 					// Atomic number
					x[iatoms_c] = atom.x; 					// x-position
					y[iatoms_c] = atom.y; 					// y-position
					z[iatoms_c] = atom.z; 					// z-position
					sigma[iatoms_c] = atom.sigma;			// standard deviation
					occ[iatoms_c] = atom.occ; 				// Occupancy
					region[iatoms_c] = abs(atom.region); 	// Region
					charge[iatoms_c] = atom.charge; 		// charge

					iatoms_c++;
				}

				resize(iatoms_c);

				get_statistic();
			}

			// set atoms
      #ifdef __CUDACC__
			template<class X>
			void set_atoms(const Atom_Data<X> &atoms, bool pbc_xy_i = false, 
			std::vector<Atom_Type<T, e_host>> *atom_type = 0, bool b_statistic=true)
			{
				resize(atoms.size());

				l_x = atoms.l_x;
				l_y = atoms.l_y;
				l_z = atoms.l_z;
				dz = atoms.dz;

				ct_na = atoms.ct_na;
				ct_nb = atoms.ct_nb;
				ct_nc = atoms.ct_nc; 

				ct_a = atoms.ct_a;
				ct_b = atoms.ct_b;
				ct_c = atoms.ct_c;

				ct_x0 = atoms.ct_x0;
				ct_y0 = atoms.ct_y0;

				amorp_lay_info.resize(atoms.amorp_lay_info.size());
				for(auto ik=0; ik<atoms.amorp_lay_info.size(); ik++)
					amorp_lay_info[ik] = atoms.amorp_lay_info[ik];

				T dl = 1e-04;
				T lx_b = l_x - dl;
				T ly_b = l_y - dl;

				size_type iatoms_c = 0;
				for(auto iatoms = 0; iatoms < atoms.size(); iatoms++)
				{
					if((!pbc_xy_i)||((atoms.x[iatoms]<lx_b) && (atoms.y[iatoms]<ly_b)))
					{
						Z[iatoms_c] = atoms.Z[iatoms];
						x[iatoms_c] = atoms.x[iatoms];
						y[iatoms_c] = atoms.y[iatoms];
						z[iatoms_c] = atoms.z[iatoms];
						sigma[iatoms_c] = atoms.sigma[iatoms];
						occ[iatoms_c] = atoms.occ[iatoms];
						region[iatoms_c] = abs(atoms.region[iatoms]);
						charge[iatoms_c] = atoms.charge[iatoms];

						iatoms_c++;
					}
				}

				resize(iatoms_c);

				if(b_statistic)
				{
					get_statistic(atom_type);
				}
			}
      #endif
		
			int  get_Z(int Z)
			{
				return (Z % 1000);
			}

			void get_Z_unique()
			{
				Z_unique.resize(c_nAtomsTypes);
				std::fill(Z_unique.begin(), Z_unique.end(), 0);

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					auto iZ = get_Z(Z[iatoms])-1;
					Z_unique[iZ] = iZ + 1;
				}

				Z_unique.erase(std::remove_if(Z_unique.begin(), Z_unique.end(), [](int v){ return v<=0; }), Z_unique.end());
			}

			// get statistic
      #ifdef __CUDACC__
			void get_statistic(std::vector<Atom_Type<T, e_host>> *atom_type_ptr);
      #endif
			DLL_PUBLIC void get_statistic();

			// Sort atoms along z-axis.
			DLL_PUBLIC void sort_by_z();

			// max z value within a region
			DLL_PUBLIC void minmax_z_by_region(T region_v, T &zr_min, T &zr_max);

			DLL_PUBLIC void validate_amorphous_parameters();

			DLL_PUBLIC void xy_recenter(T lx=0, T ly=0);

			inline
			r3d<T> to_r3d(const int &iatoms) const
			{
				return r3d<T>(x[iatoms], y[iatoms], z[iatoms]);
			}

			inline
			T norm(const int &iatoms, r3d<T> r) const
			{
				auto r_i = to_r3d(iatoms)-r;
				return r_i.norm();
			}

			inline
			T norm_pbc_xy(const int &iatoms, r3d<T> r) const
			{
				auto x_d = fabs(x[iatoms]-r.x);
				auto y_d = fabs(y[iatoms]-r.y); 
				auto z_d = z[iatoms]-r.z; 

				x_d = ::fmin(x_d, fabs(x_d-l_x));
				y_d = ::fmin(y_d, fabs(y_d-l_y));

				return (x_d*x_d + y_d*y_d + z_d*z_d);
			}

			inline
			T norm(const int &iatoms, T x_i, T y_i, T z_i) const
			{
				auto x_d = x[iatoms]-x_i; 
				auto y_d = y[iatoms]-y_i; 
				auto z_d = z[iatoms]-z_i;
				return (x_d*x_d + y_d*y_d + z_d*z_d);
			}

			inline
			T distance(const int &iatoms, T x_i, T y_i, T z_i)
			{
				return sqrt(norm(iatoms, x_i, y_i, z_i));
			}

			T l_x; 												// box length along x direction (Å)
			T l_y; 												// box length along y direction (Å)
			T l_z; 												// box length along z direction (Å)
			T dz;												// slice thickness (Å)

			int ct_na;											// number of unit cell along a
			int ct_nb;											// number of unit cell along b
			int ct_nc;											// number of unit cell along c

			T ct_a;												// length along a (Å)
			T ct_b;												// length along b (Å)
			T ct_c;												// length along c (Å)

			T ct_x0;											// reference position along x direction (Å)
			T ct_y0;											// reference position along y direction (Å)

      std::vector<Amorp_Lay_Info<T>> amorp_lay_info;	// amorphous layer information

      std::vector<int> Z;
			TVector_r x;
			TVector_r y;
			TVector_r z;
      std::vector<float> sigma;
      std::vector<float> occ;
      std::vector<int> region;
      std::vector<int> charge;

      std::vector<int> Z_unique;

			int Z_min;
			int Z_max;

			T x_min;
			T x_max;

			T y_min;
			T y_max;

			T z_min;
			T z_max;

			float sigma_min;
			float sigma_max;

			float occ_min;
			float occ_max;

			int region_min;
			int region_max;

			T R_int_min;
			T R_int_max;

			T x_mean;
			T y_mean;
			T z_mean;

			T x_std;
			T y_std;
			T z_std;

			T s_x; 			// size-x
			T s_y; 			// size-y
			T s_z; 			// size-z

			T x_int_min;
			T x_int_max;

			T y_int_min;
			T y_int_max;

			T z_int_min;
			T z_int_max;

			T s_x_int;
			T s_y_int;
			T s_z_int;

			T l_x_int;
			T l_y_int;
			T l_z_int;

		private:

			struct sort_atoms_by_z;

			struct Atom
			{
				int Z;
				T x;
				T y;
				T z;
				T sigma;
				T occ;
				int region;
				int charge;

				Atom():Z(0), x(0), y(0), z(0), sigma(0), occ(0), charge(0){};

				// check if atomic position are in the first quadrant
				bool is_xy_positive() const
				{
					const T ee = 1e-4;
					return (-ee<x) && (-ee<y);
				}
			};

			template <class X>
			Atom read_atom(const int &nr, const int &nc, X *atoms, const int &iatoms)
			{
				Atom atom;
				atom.Z = static_cast<int>(atoms[0*nr + iatoms]); 						// Atomic number
				atom.x = atoms[1*nr + iatoms]; 											// x-position
				atom.y = atoms[2*nr + iatoms]; 											// y-position
				atom.z = atoms[3*nr + iatoms]; 											// z-position
				atom.sigma = static_cast<float>((nc>4)?(atoms[4*nr + iatoms]):0.085);	// Standard deviation
				atom.occ = static_cast<float>((nc>5)?(atoms[5*nr + iatoms]):1.0); 		// Occupancy
				atom.region = static_cast<int>((nc>6)?(atoms[6*nr + iatoms]):0); 		// Region
				atom.charge = static_cast<int>((nc>7)?(atoms[7*nr + iatoms]):0);		// charge

				return atom;
			}
	};

}

#endif
