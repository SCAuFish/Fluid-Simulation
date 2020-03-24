#pragma once
#ifndef _LAPLACEEIGEN_H_
#define _LAPLACEEIGEN_H_
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <iostream>
#include <glm/glm.hpp>

class LaplaceEigen
{
public:
	// velocity grid size (2D)
	int mx, my;

	// Velocity basis
	// TODO: why does this contain 4 dimensions?
	// Answer: N basis with 2 vector directions on x and y coordinates
	glm::vec3* vel_basis;
	// double**** vel_basis;

	// velocity basis coefficients
	double* coef;

	// current velocity field
	// TODO: why is this 3D?
	glm::vec3* vfield;
	// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vfield[2];

	// Basis field eigenvalues
	double* eigs, * eigs_inv, eigs_inv_root;

	// structure coefficient matrices
	// TODO: some coefficients but what are they?
	Eigen::SparseMatrix<double>* Ck;

	// viscosity
	double visc = 0.0;
	
	// timestep
	double dt = 0.1;

	// timestep adjustment factor for particle advection
	double pdt_mult = 1.0;

	// small numerical offsets
	// TODO: is this used to avoid the existence of strange values like whole zero or one?
	double margin = 1e-7;

	// dimension
	int N, N_sqrt;

	// tables for basis field lookup from vector eigenvalue k1, k2
	int** basis_lookup_table, basis_rlookup_table;

	// Particles
	int* particle_index;
	// a list of 2D coordinates, removed the coordinates of particle tail (we'll see how it looks)
	std::vector<glm::vec3> particles;   
	int num_particles = 0;
	int particle_tail_length = 3;

	// Density field
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> density;

	// Density field size
	int dmx, dmy;

	// External forces
	double* forces_dw;
	bool forces_pending;

	LaplaceEigen();
	LaplaceEigen(int grid_resolution, int N, int density_grid_resolution);
	glm::vec3 getVelocity(int x, int y);
	void setVelocity(int x, int y, glm::vec3 value);
	glm::vec3 getVelBasis(int k, int x, int y);
	void setVelBasis(int k, int x, int y, glm::vec3 value);

	void step();
	void precompute_basis_fields();
	void precompute_dynamics();

	double coefdensity(int a1, int b1, int a2, int b2, int c, int tt);

	double*** basis_field_2d_rect(int n, int m, double amp);

	double cur_energy(); // This seems to be a non-essential function that only outputs some useless energy value

	void attract_particles();

	void add_particles(int n);

	double getInterpolatedValue(double x, double y, int index);

	// Get the velocity at a coordinate using bilinear interpolation
	double* vel_at_bilinear(double x, double y);

	// Get velocity using cubic interpolation
	double* vel_at_cubic(double x, double y);

	int clampi(int f, int a, int b);

	double clampd(double f, double a, double b);

	void advect_density();

	double density_at(double x, double y);

	void advect_particles();

	void print_field();

	int basis_lookup(int index, int component);

	int basis_rlookup(int k1, int k2);

	void expand_basis();

	void fill_lookup_table();

	double* project_forces(double** force_path);

	void stir(double** force_path);
};

#endif