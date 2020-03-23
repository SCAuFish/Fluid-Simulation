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
	this->vfield[0].resize(this->mx + 1, this->my + 1);
	this->vfield[1].resize(this->mx + 1, this->my + 1);

	this->density.resize(this->dmx, this->dmy);

	this->coef = (double*)malloc(sizeof(double) * N);
	this->forces_dw = (double*)malloc(sizeof(double) * N);

	// Fill lookup table
	std::cout << "Filling lookup table" << std::endl;
	this->fill_lookup_table();

	// Compute basis fields
	std::cout << "Precomputing Basis Fields" << std::endl;
	this->precompute_basis_fields();

	// Compute dynamics
	std::cout << "Precomputing Dynamics" << std::endl;
	this->precompute_dynamics();

	std::cout << "Done" << std::endl;
}

void LaplaceEigen::step()
{
} 

void LaplaceEigen::precompute_basis_fields()
{
}

void LaplaceEigen::precompute_dynamics()
{
}

double LaplaceEigen::coefdensity(int a1, int b1, int a2, int b2, int c, int tt)
{
	return 0.0;
}

double*** LaplaceEigen::basis_field_2d_rect(int n, int m, double amp)
{
	return nullptr;
}

double LaplaceEigen::cur_energy()
{
	return 0.0;
}

void LaplaceEigen::attract_particles()
{
}

void LaplaceEigen::add_particles(int n)
{
	this->num_particles = n;
	this->particle_index = (int*)malloc(sizeof(int) * n);

	for (int i = 0; i < this->num_particles; i++) {
		Eigen::Vector2d rand_pos((double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
		this->particles.push_back(rand_pos);
	}
}

double LaplaceEigen::getInterpolatedValue(double x, double y, int index)
{
	return 0.0;
}

double* LaplaceEigen::vel_at_bilinear(double x, double y)
{
	return nullptr;
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
}

void LaplaceEigen::print_field()
{
}

int LaplaceEigen::basis_lookup(int index, int component)
{
	return 0;
}

int LaplaceEigen::basis_rlookup(int k1, int k2)
{
	return 0;
}

void LaplaceEigen::expand_basis()
{
}

void LaplaceEigen::fill_lookup_table()
{
}

double* LaplaceEigen::project_forces(double** force_path)
{
	return nullptr;
}

void LaplaceEigen::stir(double** force_path)
{
}
