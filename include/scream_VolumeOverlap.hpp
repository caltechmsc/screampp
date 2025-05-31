#ifndef SCREAM_VOLUMEOVERLAP
#define SCREAM_VOLUMEOVERLAP

#include "defs.hpp"

#include "scream_vector.hpp"
#include "scream_matrix.hpp"

#include "scream_atom.hpp"

#include <set>

class VolumeOverlap;
class Box;

class Box {
public:
  Box();
  Box(double, double, double, double, double, double); //minX, minY, minZ, maxX, maxY, maxZ.
  Box(const Box&);
  ~Box();

  Box& operator=(const Box&);

  /* Max, Min get functions. */

  double getMaxX() const {return this->max_x;};
  double getMinX() const {return this->min_x;};
  double getMaxY() const {return this->max_y;};
  double getMinY() const {return this->min_y;};
  double getMaxZ() const {return this->max_z;};
  double getMinZ() const {return this->min_z;};

  double getSpacing() const {return this->spacing;};

  ScreamVector getXVector() const {return this->x_direction;};
  ScreamVector getYVector() const {return this->y_direction;};
  ScreamVector getZVector() const {return this->z_direction;};

  /* Grid generation functions. */
  void generateGrid(double); // Generates grid with grid size [double].
  
  /* Grid calculation functions. */
  std::set<ScreamVector> getEnclosedPoints(SCREAM_ATOM*); // returns the set of points that are enclosed by SCREAM_ATOM. 
  std::set<ScreamVector> getEnclosedPoints(ScreamAtomV&); // returns the set of points that are enclosed by list of SCREAM_ATOM's.

  std::set<ScreamVector> getEnclosedPoints(ScreamVector, double); // returns the set of points that are enclosed by sphere centered at specifed coordinate with radius double.


private:
  double max_x, min_x, max_y, min_y, max_z, min_z; // Multiples of x_direction, y_direction, z_direction to specify size of box.
  double spacing; // size of grid.
  ScreamVector x_direction, y_direction, z_direction; // Axis.  Default: [1 0 0], [0 1 0], [0 0 1]
  //vector<vector<vector<GridPoint> > > gridPoints;
  
  ScreamVector*** gridPoints; // Pointed to 3D array for Gridpoint, represented by ScreamVector.
  int xGrid, yGrid, zGrid; // 3D array size.
  
  double _distanceSquared(ScreamVector&, ScreamVector&); // Calculates distances squared between two coords.

};




class VolumeOverlap { 
  /* VolumeOverlap Class is an interface to routines that calculate the volume overlap between atoms. */
public:
  VolumeOverlap();
  ~VolumeOverlap();
  
  double volumeOverlapCalc(ScreamAtomV&, ScreamAtomV&); ///< General formulation.  Simpler to access than the one below which uses vector<double>.
  vector<double> volumeOverlapCalc(ScreamAtomV&, vector<ScreamAtomV>&); ///< Most general formulation, aside from specifying all xyz coords and radii.

private:
  Box box;
  double buffer; 

  /* Auxiliary functions. */

  Box getBoundingBox(ScreamAtomV&, vector<ScreamAtomV>&, double);

};




#endif
