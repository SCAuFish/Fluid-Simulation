Mar. 23rd
1, Since vfield and vel_basis are both representing velocity vectors at a specific position, using glm::vec3 to represent them for better scalability to 3D.
2, Write setters and getters since when initialized, vfield, vel_basis are all in single array. In the setter and getter, it would handle the dimension transition issue.
TODO:
3, Setter/getter for vel_basis have not been implemented yet.
4, fill_lookup_table and other functions that are called in LaplaceEigen constructor. These functions might explain how vel_basis is populated.

Mar. 24th
1, basis_field_2d_rect is called to "precompute_basis_fields" while constructing. There are quite a few wierd operations that needs to be clarified. For example, when precomputing basis fields, what is basis_lookup doing?
2，In step() function, energy is used. When normalizing energy, the coefficients are also changed. So it seems we need to implement energy functions as well.
TODO:
3, implement energy-related functions
4, connect step() with Window.cpp. Try run it.

Mar. 25th
1, Many double arrays can be rewritten with Eigen::vector/matrix. For example, double * eigs can be rewritten. And all its products with this->coef could be simplified instead of using for loop.
2, Initialize some random force when intializing the window. But the particles are not moving. TODO: debug in step() function.
