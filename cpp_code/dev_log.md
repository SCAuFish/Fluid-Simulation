Mar. 24th
1, Since vfield and vel_basis are both representing velocity vectors at a specific position, using glm::vec3 to represent them for better scalability to 3D.
2, Write setters and getters since when initialized, vfield, vel_basis are all in single array. In the setter and getter, it would handle the dimension transition issue.
TODO:
3, Setter/getter for vel_basis have not been implemented yet.
4, fill_lookup_table and other functions that are called in LaplaceEigen constructor. These functions might explain how vel_basis is populated.