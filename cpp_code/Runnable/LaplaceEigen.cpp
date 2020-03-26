#include "LaplaceEigen.h"

LaplaceEigen::LaplaceEigen() : LaplaceEigen(32, 16, 32)
{
	// call with default values
}

LaplaceEigen::LaplaceEigen(int grid_resolution, int N, int density_grid_resolution)
{
	this->mx = this->my = grid_resolution;
	this->dmx = this->dmy = density_grid_resolution;
	this->N = N;

	// vfield is a 2D array with velocity matrices
	// TODO: why 2D?
	// Answer: One is for velocity at x direction and the other on y direction?
	// TODO: follow-up: if so, reformat so that velocity at each position is in one vector
	//this->vfield[0].resize(this->mx + 1, this->my + 1);
	//this->vfield[1].resize(this->mx + 1, this->my + 1);
	this->vfield = new glm::vec3[(this->mx + 1) * (this->my + 1)];

	this->density.resize(this->dmx, this->dmy);

	this->coef = new Eigen::VectorXd(N);
	this->coef->setZero();
	this->forces_dw = new double[N];

	// Fill lookup table
	// TODO: in these tables, what is mapped to what?
	std::cout << "Filling lookup table" << std::endl;
	this->fill_lookup_table();

	// Compute basis fields
	// This populates vel_basis
	std::cout << "Precomputing Basis Fields" << std::endl;
	this->precompute_basis_fields();

	// Compute dynamics
	std::cout << "Precomputing Dynamics" << std::endl;
	this->precompute_dynamics();

	std::cout << "Done" << std::endl;
}

glm::vec3 LaplaceEigen::getVelocity(int x, int y) {
	return this -> vfield[x * (this->my + 1) + y];
}

void LaplaceEigen::setVelocity(int x, int y, glm::vec3 value)
{
	this->vfield[x * (this->my + 1) + y] = value;
}

glm::vec3 LaplaceEigen::getVelBasis(int k, int x, int y) 
{
	return this->vel_basis[ k * (this->mx+1) * (this->my+1) + x * (this->my+1) + y];
}

void LaplaceEigen::setVelBasis(int k, int x, int y, glm::vec3 value)
{
	this->vel_basis[k * (this->mx + 1) * (this->my + 1) + x * (this->my + 1) + y] = value;
}

void LaplaceEigen::step()
{
	// Advance the simulation
	double * dw = new double[this->N];
	double prev_e = cur_energy();   // THis energy doesn't seem very important

	double* dwt = new double[4 * this->N];
	double* qn = new double[4 * this->N];

	// First order explicit Euler
	// TODO: confirm if this is a vanilla update for dw
	for (int k = 0; k < this->N; k++) {
		dw[k] = this->coef->dot((*(this->Ck[k])) * (*(this->coef)));
	}

	// TODO: higher order integrator
	// qn[0] = this.coef translated into
	// memcpy(qn, this->coef, this->N * sizeof(double));

	for (int k = 0; k < this->N; k++) {
		(*(this->coef))[k] += dw[k] * this->dt;
	}

	// renormalize energy
	// in this function, coef is also adjusted. Seems like energy is also of some use
	if (prev_e > 1e-5) {
		this->set_energy(prev_e);
	}

	for (int k = 0; k < this->N; k++) {
		(*(this->coef))[k] *= std::exp(-1.0 * this->eigs[k] * this->dt * this->visc);
		// add external forces
		(*(this->coef))[k] += this->forces_dw[k];
		forces_dw[k] = 0;
	}

	this->expand_basis();
} 

void LaplaceEigen::precompute_basis_fields()
{
	this->vel_basis = (glm::vec3*)malloc(sizeof(glm::vec3) * N * (this->mx + 1) * (this->my + 1));

	for (int i = 0; i < this->N; i++) {
		// TODO: what is this look-up doing?
		int k1 = this->basis_lookup(i, 0);
		int k2 = this->basis_lookup(i, 1);

		// Cheng: based on how lookup table is filled, we have the relation:
		// i = (k1-1) * N_sqrt + (k2-1)

		// TODO: since our vel_basis has different formats, need to rework
		// basis_field_2d_rect function
		glm::vec3* basis = this->basis_field_2d_rect(k1, k2, 1.0);
		std::memcpy(this->vel_basis + i, basis, sizeof(glm::vec3) * (this->mx + 1) * (this->my + 1));
		delete basis;
	}
}

void LaplaceEigen::precompute_dynamics()
{
	// Precomputes structure coefficients for 2-D rectangle basis functions.
	// This is performed symbolically, see paper for details

	// Cheng: wait this is Ck, as the one mentioned in page 6?
	this->Ck = (Eigen::SparseMatrix<double>**) malloc(sizeof(Eigen::SparseMatrix<double>*) * N);

	for (int i = 0; i < N; i++) {
		this->Ck[i] = new Eigen::SparseMatrix<double>(N, N);
	}

	// Calculate eigen values for each basis field
	// Cheng: eigen values are so easy to compute?
	this->eigs = new double[N];
	this->eigs_inv = new double[N];
	this->eigs_inv_root = new double[N];

	for (int d1 = 0; d1 < N; d1++) {
		int a1 = this->basis_lookup(d1, 0);
		int a2 = this->basis_lookup(d1, 1);
		int a1_2 = a1 * a1;

		// TODO: what are these?
		double lambda_a = -(a1 * a1 + a2 * a2);
		double inv_lambda_a = -1.0 / (a1 * a1 + a2 * a2);

		for (int d2 = 0; d2 < N; d2++) {
			int b1 = this->basis_lookup(d2, 0);
			int b2 = this->basis_lookup(d2, 1);

			double lambda_b = -(b1 * b1 + b2 * b2);
			double inv_lambda_b = -1.0 / (b1 * b1 + b2 * b2);

			int k1 = this->basis_rlookup(a1, a2);
			int k2 = this->basis_rlookup(b1, b2);

			int antipairs[4][2];
			antipairs[0][0] = a1 - b1; antipairs[0][1] = a2 - b2;
			antipairs[1][0] = a1 - b1; antipairs[1][1] = a2 + b2;
			antipairs[2][0] = a1 + b1; antipairs[2][1] = a2 - b2;
			antipairs[3][0] = a1 + b1; antipairs[3][1] = a2 + b2;

			for (int c = 0; c < 4; c++) {
				int i = antipairs[c][0];
				int j = antipairs[c][1];

				int index = this->basis_rlookup(i, j);

				if (index != -1) {
					double coef = -this->coefdensity(a1, a2, b1, b2, c, 0) * inv_lambda_b;
					this->Ck[index]->coeffRef(k1, k2) = coef;
					this->Ck[index]->coeffRef(k2, k1) = coef * -lambda_b / lambda_a;
				}
			}
		}
	}
}

// TODO: completely no idea what this is doing
double LaplaceEigen::coefdensity(int a1, int b1, int a2, int b2, int c, int tt)
{
	if (tt == 0) {
		//SS x SS
		if (c == 0) return 0.25 * -(a1 * b2 - a2 * b1); // --
		if (c == 1) return 0.25 * (a1 * b2 + a2 * b1);  // -+
		if (c == 2) return 0.25 * -(a1 * b2 + a2 * b1); // +-
		if (c == 3) return 0.25 * (a1 * b2 - a2 * b1);  // ++
	}
	else if (tt == 1) {
		//SC x SS
		if (c == 0) return 0.25 * -(a1 * b2 - a2 * b1); // --
		if (c == 1) return 0.25 * -(a1 * b2 + a2 * b1);  // -+
		if (c == 2) return 0.25 * (a1 * b2 + a2 * b1); // +-
		if (c == 3) return 0.25 * (a1 * b2 - a2 * b1);  // ++
	}
	else if (tt == 2) {
		//CS x SS
		if (c == 0) return 0.25 * -(a1 * b2 - a2 * b1); // --
		if (c == 1) return 0.25 * -(a1 * b2 + a2 * b1);  // -+
		if (c == 2) return 0.25 * (a1 * b2 + a2 * b1); // +-
		if (c == 3) return 0.25 * (a1 * b2 - a2 * b1);  // ++
	}
	else if (tt == 3) {
		//CS x SS
		if (c == 0) return 0.25 * -(a1 * b2 - a2 * b1); // --
		if (c == 1) return 0.25 * -(a1 * b2 + a2 * b1);  // -+
		if (c == 2) return 0.25 * (a1 * b2 + a2 * b1); // +-
		if (c == 3) return 0.25 * (a1 * b2 - a2 * b1);  // ++
	}

	return 0;
}

glm::vec3* LaplaceEigen::basis_field_2d_rect(int n, int m, double amp)
{
	// Calculate Laplacian eigenfunction for eigenvalue (k1,k2) on 2D rectangle

	// Cheng: although I don't know what laplacian eigenfunction is... but let me
	// just take this as is
	int a = n;
	int b = m;

	// TODO: what is this?
	double xfact = 1.0;
	if (n != 0) xfact = -1.0 / (a * a + b * b);
	double yfact = 1.0;
	if (m != 0) yfact = -1.0 / (a * a + b * b);

	glm::vec3* vf = (glm::vec3*) malloc(sizeof(glm::vec3) * (this->mx + 1) * (this->my + 1));

	double deltax = 3.14 / this->mx;
	double deltay = 3.14 / this->my;

	for (int i = 0; i < this->mx + 1; i++) {
		for (int j = 0; j < this->my + 1; j++) {
			double x = i * deltax;
			double y = i * deltay;

			vf[i * (this->my + 1) + j] = glm::vec3(amp * xfact * -b * std::sin(a * x) * std::cos(b * (7 + 0.5 * deltay)),
				amp * yfact * a * std::cos(a*(x+0.5*deltax)) * std::sin(b*y), 
				0.0f);
		}
	}

	return vf;
}

double LaplaceEigen::cur_energy()
{
	double energy = 0;
	for (int i = 0; i < this->N; i++) {
		energy += this->eigs_inv[i] * (*(this->coef))[i] * (*(this->coef))[i];
	}

	return energy;
}

void LaplaceEigen::set_energy(double desired_energy)
{
	double cur_e = this->cur_energy();
	double fact = std::sqrt(desired_energy) / std::sqrt(cur_e);

	(*(this->coef)) *= fact;
}

void LaplaceEigen::attract_particles()
{
}

void LaplaceEigen::add_particles(int n)
{
	this->num_particles = n;
	this->particle_index = (int*)malloc(sizeof(int) * n);

	for (int i = 0; i < this->num_particles; i++) {
		glm::vec3 rand_pos((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, 0.0f);
		this->particles.push_back(rand_pos);
	}
}

double LaplaceEigen::getInterpolatedValue(double x, double y, int index)
{
	return 0.0;
}

glm::vec3 LaplaceEigen::vel_at_bilinear(double x, double y)
{
	return glm::vec3(0.0f);
}

double* LaplaceEigen::vel_at_cubic(double x, double y)
{
	return nullptr;
}

int LaplaceEigen::clampi(int f, int a, int b)
{
	return 0;
}

double LaplaceEigen::clampd(double f, double a, double b)
{
	return 0.0;
}

void LaplaceEigen::advect_density()
{
}

double LaplaceEigen::density_at(double x, double y)
{
	return 0.0;
}

void LaplaceEigen::advect_particles()
{
	// Advects particles using RK4 and cubic velocity interpolation

	double pdt = this->dt * this->pdt_mult;
	bool RK4 = false;
	bool RK2 = false;
	bool Euler = true;

	for (int i = 0; i < this->num_particles; i++) {
		glm::vec3 position = this->particles[i];
		
		glm::vec3 new_pos(0.0f);
		if (RK4) {
			// TODO
		}
		else if (RK2) {
			// TODO
		}
		else if (Euler) {
			glm::vec3 vel = vel_at_bilinear(position[0], position[1]);
			new_pos = position + vel * (float)pdt;
		}

		//nx = clampd(nx, margin, 1.0 - margin);
		//ny = clampd(ny, margin, 1.0 - margin);
		this->particles[i] = new_pos;
	}
}

void LaplaceEigen::print_field()
{
}

int LaplaceEigen::basis_lookup(int index, int component)
{
	// TODO: is the component here related to the fact of 2-dimensional plane?
	return this->basis_lookup_table[index * 2 + component];
}

int LaplaceEigen::basis_rlookup(int k1, int k2)
{
	return 0;
}

/**
 * This matrix is building up the vector field with velocity basis
 */
void LaplaceEigen::expand_basis()
{
	// Clean matrix
	delete this->vfield;
	this->vfield = new glm::vec3[(this->mx + 1) * (this->my + 1)];
	// Calculate superposition of basis fields
	for (int k = 0; k < this->N; k++) {
		// loop through all N basis
		for (int i = 0; i < this->mx + 1; i++) {
			for (int j = 0; j < this->my + 1; j++) {
				// Loop through all positions in the grid
				this->setVelocity(i, j, (float)(*(this->coef))[k] * this->getVelBasis(k, i, j) + this->getVelocity(i, j));
			}
		}
	}
}

void LaplaceEigen::fill_lookup_table()
{
	// Assume N is perfect square number
	this->N_sqrt = (int)std::floor(std::sqrt((double)N));

	this->basis_lookup_table = (int*) malloc(sizeof(int) * this->N * 2);
	this->basis_rlookup_table = (int*) malloc(sizeof(int) * (this->N_sqrt + 1) * (this->N_sqrt + 1));

	// Initialize lookup table to -1, meaning (k1, k2) basis field does not exist
	for (int i = 0; i < (this->N_sqrt + 1) * (this->N_sqrt + 1); i++) {
		this->basis_rlookup_table[i] = -1;
	}

	int index = 0;
	for (int k1 = 0; k1 <= this->N_sqrt; k1++) {
		for (int k2 = 0; k2 <= this->N_sqrt; k2++) {
			if (k1 > this->N_sqrt || k1 < 1 || k2 > this->N_sqrt || k2 < 1) {
				//these fields do not exist
				// Cheng: what?? why don't you start the for-loop from 1 then?
				continue;
			}

			this->basis_lookup_table[index * 2 + 0] = k1;
			this->basis_lookup_table[index * 2 + 1] = k2;

			this->basis_rlookup_table[k1 * (this->N_sqrt + 1) + k2] = index;

			index++;
		} 
	}
} 

double* LaplaceEigen::project_forces(double** force_path)
{
	return nullptr;
}

void LaplaceEigen::stir(double** force_path)
{
}
