#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iomanip>
using namespace std;
#include "random.h"
#include "util.h"
#include "ellipsoid.h"
#include "omp.h"

class Ellipsoid : public Ellipsoid_basic {
	vector f, fu;
public:
	Ellipsoid (Species* k, vector box) : Ellipsoid_basic(k, box) {}
	void set_force(vector ff = vector()) { f = ff; }
	void set_forceu(vector ffu = vector()) { fu = ffu; }
	vector& get_f() { return f; }
	vector& get_fu() { return fu; }
};

class Packing {
	bool Lend, Lreslt_x;
	ifstream pDatin;
	ofstream pOutput;
	ofstream pResltx;

	int No_parts, No_species;
	double Epsilon, Epsilon_scl, Eps_rot;
	Ellipsoid **parts;
	Species **species;

	double Force_step;
	double Din, Dout0, Dout, Dout2;
	double Pactual, Pnom0, Pnomin;

	int No_cells_x, No_cells_y, No_cells_z;
	int No_cells, Ncell_min;
	int Max_steps, Nstep, Ntau;
	double Power, Sphere_vol, Diam_dens;
	double Rc_max, Rcut, Rcut2, Relax;
	double Diam_incr;
	vector Box, Half;

	int *Link_head, *Link_list;
	int *Ncell_bound_x, *Ncell_bound_y, *Ncell_bound_z;
	double *Pbc_x, *Pbc_y, *Pbc_z;

	bool Lkill, Lprec_chan, Leq_vol;
	int Nfig, Nsf;
	int Nlines, No_end_page, No_page, Npage_len;
	int Nprint_step;
	int Nrslt_step;
	int Nrot_step;
	

	void data();
	void auxiliary_val();
	void stepon();
	void forces();
	void motion();
	void force_part(int ipart_p);
	void force_all(int ipart_p);
	void output(int kind_p);
	void save_x();
	double zeroin(Ellipsoid_2& ell, double ax, double bx);

	double calc_forces(Ellipsoid_2& e, Ellipsoid& ei, Ellipsoid& ej) {
//		double lambda = e.newton();
		double lambda = zeroin(e, 0., 1.);
		e.get_lambda() = lambda;
//		cout << "lambda = " << lambda << " " << e.f_AB(lambda) << " " << e.f_AB_d(lambda) << endl;
		if (fabs(lambda) < 1e-10) {
			cout << "e1" << endl << ei << endl;
			cout << "e2" << endl << ej << endl;
			exit(1);
		}
		double f_AB_scl = 4 * e.f_AB(lambda);
		double f_AB = f_AB_scl / Dout2;
		vector n = e.count_n();
		if (f_AB >= 1) return f_AB;
		if (f_AB <= 0) {
			cout << "fab " << f_AB << endl;
			cout << ei << endl;
			cout << ej << endl;
			cout << e << endl;
			cout << "-----------------------------------------------" << endl;
		}
//	shift
		if (f_AB_scl < Din) Din = f_AB_scl;
		vector f = (1 - f_AB) * norm(n);
		ei.get_f() -= f;
		ej.get_f() += f;
/*
		if (&ei == parts[0] || &ej == parts[0]) {
//			cout << "ei " << ei;
//			cout << "ej " << ej;
//			cout << "fAB: " << f_AB << endl;
//			cout << "r: " << e.get_r() << endl;
			cout << "f: " << f << endl;
		}
		double psph = 1 - e.get_r()*e.get_r() / Dout2;
		vector fsph = psph * e.get_r().norm();
		vector d = f - fsph;
		if (d*d > 1e-4) {
			cout << "ei: " << ei << endl;
			cout << "ej: " << ej << endl;
			cout << "f: " << f << endl;
			cout << "lambda: " << lambda << endl;
			cout << "fAB: " << f_AB << endl;
			cout << "psph: " << psph << endl;
			cout << "fsph: " << fsph << endl;
		}
*/
//	rotation
		if (Eps_rot < 1e-10 || Nstep % Nrot_step) return f_AB;
//		cout << "RAC " << e.count_rac() << endl;
//		cout << "RBC " << e.count_rbc() << endl;
//		cout << "n " << n << endl;
//		cout << "e " << e << endl;
		if (!ei.get_k()->gets())
			ei.get_fu() -= (1 - f_AB) * norm(prod(e.count_rac(), n));
		if (!ej.get_k()->gets())
			ej.get_fu() += (1 - f_AB) * norm(prod(e.count_rbc(), n));
		return f_AB;
	}
public:
	Packing (int, char **);
	~Packing();
	void control();
};

Packing::Packing (int argc_p, char *argv_p[]) {
	if (argc_p < 3) {
		cout << "Usage: " << argv_p[0] << " input output [x-output]";
		cout << endl;
		exit(1);
	}
	pDatin.open(*++argv_p);
	if (!pDatin) open_failure(cerr, *argv_p);
	data();
	pOutput.open (*++argv_p);
	if (!pOutput) open_failure(cerr, *argv_p);
	pOutput.setf (ios::fixed, ios::floatfield);
	cout.setf (ios::fixed, ios::floatfield);
	if (Lreslt_x) {
		if (argc_p < 4) {
			cout << "Usage: " << argv_p[0] << " input output x-output";
			cout << endl;
			exit(1);
		}
		pResltx.open(*++argv_p);
		if (!pResltx) open_failure(cerr, *argv_p);
	}
	Din = 0;
	Lend = 0;
	Lkill = 0;
	Lprec_chan = 0;
	Nfig = 8;
	Nlines = 0;
	No_end_page = 0;
	No_page = 0;
	Nsf = 1;
	Power = 1./3.;
	Sphere_vol = M_PI / 6.;
	Random::init();
}

Packing::~Packing() {
	delete[] Link_head;
	delete[] Link_list;
	delete[] Ncell_bound_x;
	delete[] Ncell_bound_y;
	delete[] Ncell_bound_z;
	delete[] Pbc_x;
	delete[] Pbc_y;
	delete[] Pbc_z;
	for (int i = 0; i < No_parts; i++) delete parts[i];
	for (int i = 0; i < No_species; i++) delete species[i];
	delete[] parts;
	delete[] species;
	if (Lreslt_x) pResltx.close();
	pOutput.close ();
}

void Packing::control() {
	cout << "control begin" << endl;
	auxiliary_val();
	cout << "aux end" << endl;
	forces();
	cout << "forces end" << endl;
	double start = omp_get_wtime( );
	Pactual = Din*Din*Din / Diam_dens;
	output(1);
	for (Nstep=0; Nstep < Max_steps && !Lend; Nstep++) {
		stepon();
//		if (Lreslt_x && Nstep % Nrslt_step == 0) output(4);
		if (Nstep % Nprint_step == 0) output(2);
	}
	double end = omp_get_wtime( );
	cout<<"start "<<start << endl<< "koniec " <<end<<endl<<"czas pomiaru "<<end - start<<endl; 
	output(3);
	if (Lreslt_x) output(4);
}

void Packing::data() {
	omit_line(pDatin);
	Lreslt_x = get_bool(pDatin);
	omit_line(pDatin);
	No_parts = get_int(pDatin);
	No_species = get_int(pDatin);
	Epsilon = get_double(pDatin);
	Eps_rot = get_double(pDatin);
	Diam_incr = get_double(pDatin);
//	double dummy = get_double(pDatin);
	No_cells_x = get_int(pDatin);
	No_cells_y = get_int(pDatin);
	No_cells_z = get_int(pDatin);
	Ntau = get_int(pDatin);
	Pnom0 = get_double(pDatin);
	Max_steps = get_int(pDatin);
	Leq_vol = get_bool(pDatin);
	if (No_parts <= 0 || No_cells_x < 0 ||
		No_cells_y < 0 || No_cells_z < 0 ||
		Max_steps < 0) error(cerr, "Wrong input data -1-");
	species = new Species*[No_species];
	int np = 0;
	for (int i = 0; i < No_species; i++) {
		int n = get_int(pDatin);
		double r0 = get_double(pDatin);
		double r1 = get_double(pDatin);
		double r2 = get_double(pDatin);
		species[i] = new Species (n, vector(r0, r1, r2));
		np += n;
	}
	if (np != No_parts) error (cerr, "Wrong input data -2-");
	omit_line(pDatin);
	Npage_len = get_int(pDatin);
	Nprint_step = get_int(pDatin);
	Nrslt_step = get_int(pDatin);
	Nrot_step = get_int(pDatin);
	if (Npage_len <= 0 || Nprint_step <= 0 || 
		Nrslt_step <= 0) error(cerr, "Wrong input data -3-");
	pDatin.close();
}

void Packing::auxiliary_val() {
	int imult;
	int ndim_x, ndim_y, ndim_z;
	No_cells = No_cells_x * No_cells_y * No_cells_z;
	Ncell_min = (No_cells_x < No_cells_y) ?  No_cells_x : No_cells_y;
	Ncell_min = (No_cells_z < Ncell_min) ? No_cells_z : Ncell_min;
	Box = vector(No_cells_x, No_cells_y, No_cells_z);
	Half = 0.5 * Box;
	Diam_dens = No_cells / Sphere_vol;
	double corr = 0;
	vector r;
	Rc_max = 0;
	for (int is = 0; is < No_species; is++) {
		r = species[is]->getr();
		corr += species[is]->getn() * r[0]*r[1]*r[2];
		if (r[0] > Rc_max) Rc_max = r[0];
	}
	Diam_dens /= corr;
	Pnomin = Pnom0;
	Dout0 = pow (Diam_dens * Pnom0 , Power);
	Dout = Dout0;
	Relax = 0.5 * Dout / Ntau;
	Dout2 = Dout * Dout;

	parts = new Ellipsoid*[No_parts];
	int ipart = 0;
	for (int is = 0; is < No_species; is++)
		for (int ip = 0; ip < species[is]->getn(); ip++)
			parts[ipart++] = new Ellipsoid (species[is], Box);

	Link_head = new int [No_cells];
	Link_list = new int [No_parts];

	ndim_x = 3 * No_cells_x - 1;
	ndim_y = 3 * No_cells_y - 1;
	ndim_z = 3 * No_cells_z - 1;

	Ncell_bound_x = new int [ndim_x];
	Ncell_bound_y = new int [ndim_y];
	Ncell_bound_z = new int [ndim_z];

	Pbc_x = new double [ndim_x];
	Pbc_y = new double [ndim_y];
	Pbc_z = new double [ndim_z];

	imult = 1;
	corr = -No_cells_x;
	for (int i = 1; i < ndim_x; i++) {
		Ncell_bound_x[i] = imult * No_cells_y * No_cells_z;
		Pbc_x[i] = corr;
		if (++imult >= No_cells_x) {
			imult = 0;
			corr += No_cells_x;
		}
	}
	imult = 1;
	corr = -No_cells_y;
	for (int i = 1; i < ndim_y; i++) {
		Ncell_bound_y[i] = imult * No_cells_z;
		Pbc_y[i] = corr;
		if (++imult >= No_cells_y) {
			imult = 0;
			corr += No_cells_y;
		}
	}
	imult = 1;
	corr = -No_cells_z;
	for (int i = 1; i < ndim_z; i++) {
		Ncell_bound_z[i] = imult;
		Pbc_z[i] = corr;
		if (++imult >= No_cells_z) {
			imult = 0;
			corr += No_cells_z;
		}
	}

	if (Lreslt_x) {
		pResltx.write((char*)&No_parts, sizeof(int));
		pResltx.write((char*)&Box, sizeof(vector));
		pResltx.write((char*)&Pactual, sizeof(double));
		pResltx.write((char*)&Din, sizeof(double));
		pResltx.write((char*)&No_species, sizeof(int));
		pResltx.write((char*)species, sizeof(Species*));
		for (int ip = 0; ip < No_species; ip++)
			pResltx.write((char *)species[ip], sizeof(Species));
		flush (pResltx);
	}
}

void Packing::stepon() {
	int iprec;
	Dout -= Relax;
	Dout2 = Dout * Dout;
	motion();
	forces();
	Pnomin = Dout * Dout * Dout / Diam_dens;
	Pactual = Din * Din * Din / Diam_dens;
	Lend = (Dout <= Din);
	if (Lend) {
		double dsave = Dout;
		Dout = 1.1 * Din;
		Dout2 = Dout*Dout;
		output(2);
		forces();
		Pactual = Din * Din * Din / Diam_dens;
		Dout = dsave;
		Dout2 = Dout*Dout;
		output(2);
		return;
	}
	if (Lkill) return;
	iprec = (int) (-log10 (Pnomin - Pactual));
	if (iprec >= Nsf) {
		Relax *= 0.5;
		Lprec_chan = true;
		blines (pOutput);
		Nlines++;
		if (++Nsf >= Nfig) Lkill = true;
	}
}

void Packing::forces() {
//	cout << "forces" << endl;
	Din = Dout * Dout; 
	int i;
	// 3 pragma gave nothing
	#pragma omp parallel  for private(i) 
	for (i = 0; i < No_cells; i++) 
		Link_head[i] = -1;
	// 2 pragma gave nothing	
	#pragma omp parallel  for private(i) 
	for (i = 0; i < No_parts; i++) {
		parts[i]->set_force();
		parts[i]->set_forceu();
	}
	// 1 pragma gave 3 min 38 sek
	int ipart;
	
	#pragma omp parallel  for private(ipart) schedule(static) 
	for (ipart = 0; ipart < No_parts; ipart++) {
		Species* k = parts[ipart]->get_k();
		Rcut = 0.55 * Dout * ((k->getr())[0] + Rc_max);
		Rcut2 = Rcut * Rcut;
//		cout << "i = " << ipart << endl;
		if ((int) (2.0 * Rcut) < Ncell_min - 1) force_part(ipart);
		else force_all(ipart);
	}
	Din = sqrt(Din);
}

void Packing::motion() {
//	cout << "motion" << endl;
	vector buff;
	Force_step = 0;
	Epsilon_scl = Epsilon * Dout0;
	int i;
	//4 pragma gave 3 min 24 sec
	#pragma omp parallel for private(i,buff)   
	for (i = 0; i < No_parts; i++) {
		Ellipsoid *p = parts[i];
		Force_step += sqrt(p->get_f()*p->get_f());
//	shift
		buff = p->get_pos() + p->get_f() * Epsilon_scl;
		buff.pbc(Box);
		p->get_pos() = buff;
//		if (i == 0) cout << p->get_f() << endl << p->get_pos() << endl;

//	rotation
		if (Nstep % Nrot_step) continue;
		if (fabs(Eps_rot) < 1e-5) continue;
		if (p->get_k()->gets()) continue;	//	sphere
//		cout << i << " " << p->get_fu() << endl;
		double len = Eps_rot * sqrt(p->get_fu()*p->get_fu());
		double cp = 1 / sqrt(1 + len*len);
		double sp = len * cp;
		double cp2 = sqrt(0.5 * (1 + cp));
		double sp2 = 0.5 * sp / cp2;
		Quaternion q(cp2, sp2*p->get_fu());
		p->rotate(q);
	}
//	cout << *parts[0] << endl;
	Force_step *= Epsilon_scl / (Dout0 * No_parts);
}

void Packing::force_part(int ipart_p) {
	int	jpart, ifirst;
	int	icell_x, icell_xy, icell_y, icell_z, icell;
	int leap_x, leap_y, leap_z;
	Ellipsoid *pi = parts[ipart_p], *pj;
	vector rij;
	vector pos = pi->get_pos();

	int low_cell_x = (int) (pos[0] + No_cells_x - Rcut);
	int low_cell_y = (int) (pos[1] + No_cells_y - Rcut);
	int low_cell_z = (int) (pos[2] + No_cells_z - Rcut);
	int lim_cell_x = (int) (pos[0] + No_cells_x + Rcut);
	int lim_cell_y = (int) (pos[1] + No_cells_y + Rcut);
	int lim_cell_z = (int) (pos[2] + No_cells_z + Rcut);
	int idif_x = (Ncell_bound_x[lim_cell_x] > Ncell_bound_x[low_cell_x]) ? 0 : 4;
	int idif_y = (Ncell_bound_y[lim_cell_y] > Ncell_bound_y[low_cell_y]) ? 0 : 2;
	int idif_z = (Ncell_bound_z[lim_cell_z] > Ncell_bound_z[low_cell_z]) ? 0 : 1;
	int idif = idif_x + idif_y + idif_z;
	switch (idif) {
	case 0:
	#pragma omp  parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 1:
	#pragma omp parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_z(Pbc_z[leap_z]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 2:
	#pragma omp parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_y(Pbc_y[leap_y]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 3:
	#pragma omp parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_y(Pbc_y[leap_y]);
						rij.pbc_z(Pbc_z[leap_z]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 4:
	#pragma omp parallel
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_x(Pbc_x[leap_x]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 5:
	#pragma omp parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_x(Pbc_x[leap_x]);
						rij.pbc_z(Pbc_z[leap_z]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 6: 
	#pragma omp parallel
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_x(Pbc_x[leap_x]);
						rij.pbc_y(Pbc_y[leap_y]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	case 7:
	#pragma omp parallel 
		for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++) {
			icell_x = Ncell_bound_x[leap_x];
			for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++) {
				icell_xy = Ncell_bound_y[leap_y] + icell_x;
				#pragma omp for private(leap_z)
				for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++) {
					icell = Ncell_bound_z[leap_z] + icell_xy;
					jpart = Link_head[icell];
					while (jpart != -1) {
						pj = parts[jpart];
						rij = pj->get_pos() - pi->get_pos();
						rij.pbc_x(Pbc_x[leap_x]);
						rij.pbc_y(Pbc_y[leap_y]);
						rij.pbc_z(Pbc_z[leap_z]);
						Ellipsoid_2 eij(*pi, *pj, rij);
						calc_forces(eij, *pi, *pj);
						jpart = Link_list[jpart];
					}
				}
			}
		}
		break;
	}
	icell_x = (int) pos[0];
	icell_y = (int) pos[1];
	icell_z = (int) pos[2];
	icell = No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
	ifirst = Link_head[icell];
	Link_list[ipart_p] = ifirst;
	Link_head[icell] = ipart_p;
}

void Packing::force_all(int ipart_p) {
	Ellipsoid *pi = parts[ipart_p], *pj;
	vector pos = pi->get_pos();
	vector rij;
	int jpart;
	// pragma 5 gave 3 min 13 sek

	#pragma omp parallel for private(jpart) 
	for (jpart = 0; jpart < ipart_p; jpart++) {
		pj = parts[jpart];
		rij = pj->get_pos() - pi->get_pos();
		rij.pbc_diff(Half, Box);
		Ellipsoid_2 eij(*pi, *pj, rij);
		calc_forces(eij, *pi, *pj);
	}
	int icell_x = (int) pos[0];
	int icell_y = (int) pos[1];
	int icell_z = (int) pos[2];
	int icell = No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
	int ifirst = Link_head[icell];
	Link_list[ipart_p] = ifirst;
	Link_head[icell] = ipart_p;
}

void Packing::output(int kind_p) {
	if (kind_p == 1) {
		page(pOutput);
		pOutput << "--------------------- SYSTEM SPECIFICATI";
		pOutput << "ON -------------------------------------";
		blines (pOutput, 2);
		ivar (pOutput, "No_parts", No_parts, "NUMBER OF SPHERES");
		ivar (pOutput, "No_species", No_species, "NUMBER OF SPECIES");
		rvar (pOutput, "Epsilon", Epsilon, "COEFFICIENT FOR POTENTIAL CALCULATION");
		rvar (pOutput, "Eps_rot", Eps_rot, "COEFFICIENT FOR ROTATION");
		rvar (pOutput, "Din", Din, "INNER DIAMETER");
		rvar (pOutput, "Dout0", Dout0, "INITIAL VALUE OF OUTER DIAMETER");
		rvar (pOutput, "Dout", Dout, "OUTER DIAMETER");
		rvar (pOutput, "Pactual", Pactual, "ACTUAL PACKING DENSITY");
		rvar (pOutput, "Pnom0", Pnom0, "INITIAL VALUE OF NOMINAL PACKING DENSITY");
		rvar (pOutput, "Pnomin", Pnomin, "NOMINAL PACKING DENSITY");
		ivar (pOutput, "No_cells_x", No_cells_x, "NUMBER OF CELLS IN X DIRECTION");
		ivar (pOutput, "No_cells_y", No_cells_y, "NUMBER OF CELLS IN Y DIRECTION");
		ivar (pOutput, "No_cells_z", No_cells_z, "NUMBER OF CELLS IN Z DIRECTION");
		lvar (pOutput, "Max_steps", Max_steps, "MAXIMAL NUMBER OF STEPS");
		lvar (pOutput, "Ntau", Ntau, "OUTER DIAMETER CONTRACTION RATE");
		for (int is = 0; is < No_species; is++) {
			pOutput << "SPECIES " << is << endl;
			ivar (pOutput, "Number", species[is]->getn(), "NUMBER OF PARTICLES");
			rvar (pOutput, "Diameter", (species[is]->getr())[0], "DIAMETER OF PARTICLES");
			rvar (pOutput, "Diameter", (species[is]->getr())[1], "DIAMETER OF PARTICLES");
			rvar (pOutput, "Diameter", (species[is]->getr())[2], "DIAMETER OF PARTICLES");
		}
		blines (pOutput);
		pOutput << "--------------------- I/O SPECIFICATION ";
		pOutput << "----------------------------------------";
		blines (pOutput, 2);
		ivar (pOutput, "Npage_len", Npage_len, "NUMBER OF LINES PER PAGE");
		ivar (pOutput, "Nprint_step", Nprint_step, "NUMBER OF STEPS BETWEEN DATA PRINTOUT");
		ivar (pOutput, "Nrot_step", Nrot_step, "NUMBER OF STEPS BETWEEN ROTATIONS");
		lvar (pOutput, "Nrslt_step", Nrslt_step, "NUMBER OF STEPS BETWEEN DATA STORAGE");
		blines (pOutput);
		date_time(pOutput);
		out_char (pOutput, '-', 80);
		blines (pOutput);

		date_time(cout);
		cout << endl;
		cout << "     STEP";
		cout << "        ACTUAL";
		cout << "       NOMINAL";
		cout << "         INNER";
		cout << "         OUTER";
		cout << "  FORCE MODULE";
		cout << endl;
		cout << "         ";
		cout << "       DENSITY";
		cout << "       DENSITY";
		cout << "      DIAMETER";
		cout << "      DIAMETER";
		cout << "  PER PARTICLE";
		cout << endl;
		cout << endl;
		return;
	}
	if (kind_p == 2) {
		if (Lprec_chan) {
			Lprec_chan = false;
			Nlines += 3;
			blines (pOutput);
			out_char (pOutput, '-',79);
			blines (pOutput, 2);
		}
		Nlines++;
		if (Nlines > No_end_page) {
			No_end_page += Npage_len;
			No_page++;
			Nlines += 7;
			cout.precision(10);
			cout << setw(9) << Nstep;
			cout << setw(14) << Pactual;
			cout << setw(14) << Pnomin;
			cout << setw(14) << Din;
			cout << setw(14) << Dout;
			cout << setw(14) << Force_step;
			cout << endl;
			page (pOutput);
			date_time(pOutput);
			pOutput << "V SEP/1986	MB	 ";
			out_char (pOutput, '-', 60);
			blines (pOutput);
			pOutput << "     STEP";
			pOutput << "        ACTUAL";
			pOutput << "       NOMINAL";
			pOutput << "         INNER";
			pOutput << "         OUTER";
			pOutput << "  FORCE MODULE";
			blines (pOutput);
			pOutput << "         ";
			pOutput << "       DENSITY";
			pOutput << "       DENSITY";
			pOutput << "      DIAMETER";
			pOutput << "      DIAMETER";
			pOutput << "  PER PARTICLE";
			blines (pOutput, 2);
		}
		if (Nsf <= 1) pOutput.precision(4);
		else if (Nsf <= 6) pOutput.precision(Nsf+3);
		else pOutput.precision(10);
		pOutput << setw(9) << Nstep;
		pOutput << setw(14) << Pactual;
		pOutput << setw(14) << Pnomin;
		pOutput << setw(14) << Din;
		pOutput << setw(14) << Dout;
		pOutput << setw(14) << Force_step;
		blines (pOutput);
		return;
	}
	if (kind_p == 3) {
		blines (pOutput);
		Nlines += 7;
		out_char (pOutput, '-', 80);
		blines (pOutput);
		pOutput.precision(13);
		pOutput << setw(35) << Pactual;
		pOutput << setw(19) << Din;
		pOutput << setw(26) << Nstep;
		blines (pOutput);
		pOutput << "    SIMULATION  ";
		out_char (pOutput, ' ',  5);
		out_char (pOutput, ':',  1);
		out_char (pOutput, ' ', 18);
		out_char (pOutput, ':',  1);
		out_char (pOutput, ' ', 38);
		out_char (pOutput, ':',  1);
		blines (pOutput);
		pOutput << "    STATUS    >>>>   : STATIC DENSITY";
		pOutput << "   : SPHERE DIAMETER           ";
		pOutput << "ITERATIONS :";
		blines (pOutput);
		out_char (pOutput, '-', 80);
		blines (pOutput, 3);
		out_char (pOutput, '-', 29);
		if (Lend) pOutput << " END  OF  F B C P 3 D ";
		else pOutput << " FBCP3D -- TERMINATED ";
		out_char (pOutput, '-', 29);
		blines (pOutput, 2);
		date_time (pOutput);
		blines	(cout, 2);
		out_char (cout, '-', 29);
		if (Lend) cout << " END  OF  F B C P 3 D ";
		else cout << " FBCP3D -- TERMINATED ";
		out_char (cout, '-', 29);
		blines	(cout, 2);
		date_time(cout);
		blines	(cout, 2);

//	prepare pov file 

		ofstream povf ("ellipsoid.pov");
		povf.setf (ios::fixed, ios::floatfield);
		povf.precision (6);
		for (int i = 0; i < No_parts; i++) {
			Ellipsoid *pi = parts[i];
			vector pos = pi->get_pos();
			Species* k = pi->get_k();
			vector rad = 0.5 * Din * k->getr();
			bool sph = k->gets();
			povf << "sphere {<0, 0, 0>, 1" << endl;
				povf << "scale <" <<
				setw(12) << rad[0] << ", " <<
				setw(12) << rad[1] << ", " << 
				setw(12) << rad[2] << "> " << endl;
			if (!sph) {	//	sphere
				matrix Q = pi->get_q().countQ();
				vector euler(atan2(Q(1,2),Q(2,2)), -asin(Q(0,2)), atan2(Q(0,1),Q(0,0)));
				euler *= 180 / M_PI;
				povf << "rotate " << setw(12) << euler[0] << "*x" << endl;
				povf << "rotate " << setw(12) << euler[1] << "*y" << endl;
				povf << "rotate " << setw(12) << euler[2] << "*z" << endl;
			}
			povf << "translate <" << 
				setw(12) << pos[0] << ", " <<
				setw(12) << pos[1] << ", " << 
				setw(12) << pos[2] << "> " << endl;
			povf << "}" << endl;
		}
		povf.close();
		return;
	}
	if (kind_p == 4) {
		save_x ();
		return;
	}
}

void copy(Ellipsoid_basic &eb, const Ellipsoid& e) {
	eb.set_k(e.get_k());
	eb.get_pos() = e.get_pos();
	eb.get_q() = e.get_q();
}

void Packing::save_x() {
	Ellipsoid_basic *p = new Ellipsoid_basic[No_parts];
	for (int i = 0; i < No_parts; i++) copy(p[i], *parts[i]);
//	if (Lend) {
	pResltx.seekp(sizeof(int) + 3*sizeof(double));
	pResltx.write ((char *) &Pactual, sizeof(double));
	pResltx.write ((char *) &Din, sizeof(double));
//	}
	pResltx.seekp(0l, ios::end);
	pResltx.write ((char *) &Nstep, sizeof(int));
	pResltx.write ((char *) p, sizeof(Ellipsoid_basic)*No_parts);
	flush(pResltx);
	delete[] p;
}

/*
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

double Packing::zeroin(Ellipsoid_2& ell, double ax, double bx) {
	double a, b, c;
	double fa, fb, fc;
	double EPSILON = 1e-15;
	double tol = 0.0;

	a = ax;
	b = bx;
	c = a;
	fa = ell.f_AB_d(a);
	fb = ell.f_AB_d(b);
	fc = fa;

	for(;;) {
		double prev_step = b-a;
		double tol_act;
		double p;
		double q;
		double new_step;
		if( fabs(fc) < fabs(fb) ) {
			a = b;
			b = c;
			c = a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol_act = 2*EPSILON*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if( fabs(new_step) <= tol_act || fb == (double)0 )
			return b;

		if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
			register double t1, cb, t2;
			cb = c-b;
			if (a==c) {
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			} else {
				q = fa/fc;
				t1 = fb/fc;
				t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if (p > 0.) q = -q;
			else p = -p;

			if (p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2))
				new_step = p/q;
		}

		if (fabs(new_step) < tol_act) {
			if (new_step > 0.) new_step = tol_act;
			else new_step = -tol_act;
		}
		a = b;
		fa = fb;
		b += new_step;
		fb = ell.f_AB_d(b);
		if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
			c = a;
			fc = fa;
		}
	}
}


int main(int argc, char *argv[]) {
	cout<<"asd3";
	Packing pack (argc,argv);
	pack.control();
	return 0;
}

