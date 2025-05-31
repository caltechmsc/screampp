#include "defs.hpp"
#include "scream_VolumeOverlap.hpp"
#include <cmath>
#include <set>
#include <algorithm>
#include <cmath>

VolumeOverlap::VolumeOverlap() {

}

VolumeOverlap::~VolumeOverlap() {

}

double VolumeOverlap::volumeOverlapCalc(ScreamAtomV& ref, ScreamAtomV& query) {

  vector<ScreamAtomV> queryList;
  queryList.push_back(query);
  
  /* Make box. */
  Box box = this->getBoundingBox(ref, queryList, 4.5);

  /* Generates a grid with 0.5 spacing.*/
  box.generateGrid(0.5);

  /* Volume calculations. */
  std::set<ScreamVector> refSet;
  std::set<ScreamVector> querySet;
  std::set<ScreamVector> intersectionSet;;

  refSet = box.getEnclosedPoints(ref);
  querySet = box.getEnclosedPoints(query);

  std::set_intersection(refSet.begin(), refSet.end(), 
			querySet.begin(), querySet.end(), 
			std::inserter(intersectionSet, intersectionSet.end()) );

  return double(double(intersectionSet.size()) / double(refSet.size()));
  
}

vector<double> VolumeOverlap::volumeOverlapCalc(ScreamAtomV& ref, vector<ScreamAtomV>& queryList) {

  /* Make box. */
  Box box = this->getBoundingBox(ref, queryList, 4.5);

  /* Generates a grid with 0.5 spacing.*/
  box.generateGrid(0.5);

  /* Volume calculations. */
  std::set<ScreamVector> refSet;
  std::set<ScreamVector> querySet;
  std::set<ScreamVector> intersectionSet;;

  std::vector<double> volumeOverlapValues;

  refSet = box.getEnclosedPoints(ref);

  for (vector<ScreamAtomV>::iterator query = queryList.begin(); query != queryList.end(); query++) {
    querySet = box.getEnclosedPoints(*query);
    std::set_intersection(refSet.begin(), refSet.end(), 
			  querySet.begin(), querySet.end(), 
			  std::inserter(intersectionSet, intersectionSet.end()) );
    volumeOverlapValues.push_back(double(double(intersectionSet.size()) / double(refSet.size())) );
    querySet.clear();
    intersectionSet.clear();
  }

  return volumeOverlapValues;
}

/* Auxiliary functions */

Box VolumeOverlap::getBoundingBox(ScreamAtomV& ref, vector<ScreamAtomV>& queryList, double buffer) {

  /* Comment: This is the most straightforward version, where the Box is oriented in the [1 0 0], [0 1 0], [0 0 1] directions. 
     More complex would be to find the directions in which the volume of the box is minimized.  This requires:
     1) Store the points at which maxX's and minX's are recorded.
     2) Find the orientation that minimizes volume of box.  These orientation would be defined by x_direction's.
     3) 
  */

  double max_x;
  double max_y;
  double max_z;
  double min_x;
  double min_y;
  double min_z;

  /* Initialization */
  {
    ScreamAtomVItr itr = ref.begin();
    max_x = (*itr)->getX();
    max_y = (*itr)->getY();
    max_z = (*itr)->getZ();
    min_x = max_x;
    min_y = max_y;
    min_z = max_z;
  }

  /* Now find mins and maxs. */

  double x, y, z;

  for (ScreamAtomVItr itr = ref.begin(); itr != ref.end(); itr++) { // max, min for ref.
    x = (*itr)->getX();
    y = (*itr)->getY();
    z = (*itr)->getZ();
  
    if (x > max_x) {      max_x = x;    } else if (x < min_x) {      min_x = x;    }
    if (y > max_y) {      max_y = y;    } else if (y < min_y) {      min_y = y;    }
    if (z > max_z) {      max_z = z;    } else if (z < min_z) {      min_z = z;    }

  }


  for (vector<ScreamAtomV>::iterator qItr = queryList.begin(); qItr != queryList.end(); qItr++) { // max , min for query.
    for (ScreamAtomVItr itr = qItr->begin(); itr != qItr->end(); itr++) {
      x = (*itr)->getX();
      y = (*itr)->getY();
      z = (*itr)->getZ();
      
      if (x > max_x) {      max_x = x;    } else if (x < min_x) {      min_x = x;    }
      if (y > max_y) {      max_y = y;    } else if (y < min_y) {      min_y = y;    }
      if (z > max_z) {      max_z = z;    } else if (z < min_z) {      min_z = z;    }

    }
  }

  return Box(min_x - buffer, min_y - buffer, min_z - buffer, max_x+ buffer, max_y+ buffer, max_z+ buffer);

}


Box::Box() {
  
}

Box::Box(double x_min, double y_min, double z_min, double x_max, double y_max, double z_max) {
  
  this->max_x = x_max;
  this->max_y = y_max;
  this->max_z = z_max;
  this->min_x = x_min;
  this->min_y = y_min;
  this->min_z = z_min;

  this->x_direction = ScreamVector(1,0,0);
  this->y_direction = ScreamVector(0,1,0);
  this->z_direction = ScreamVector(0,0,1);

  this->gridPoints = NULL;

}

Box::Box(const Box& box) {
  this->operator=(box);
}

Box::~Box() {
  if (this->gridPoints != NULL) {    for (int i = 0; i <= this->xGrid; i++) {
      if (this->gridPoints[i] != NULL) {	for (int j = 0; j <= this->yGrid; j++) {
	delete [] this->gridPoints[i][j];
      } }
      delete [] this->gridPoints[i];
  } }
  delete[] this->gridPoints;
}

Box& Box::operator=(const Box& box) {
  this->max_x = box.getMaxX();
  this->max_y = box.getMaxY();
  this->max_z = box.getMaxZ();
  this->min_x = box.getMinX();
  this->min_y = box.getMinY();
  this->min_z = box.getMinZ();

  this->x_direction = box.getXVector();
  this->y_direction = box.getYVector();
  this->z_direction = box.getZVector();
}

void Box::generateGrid(double spacing) {

  ScreamVector Min(this->min_x, this->min_y, this->min_z);
  ScreamVector Max(this->max_x, this->max_y, this->max_z);

  this->spacing = spacing;

  this->xGrid = int(ceil((this->max_x - this->min_x)/spacing)); // counting from 0; i.e., 0, 1 ... xGrid would cover the X range.
  this->yGrid = int(ceil((this->max_y - this->min_y)/spacing));
  this->zGrid = int(ceil((this->max_z - this->min_z)/spacing));

  /*   Init 3D array.  */
  //ScreamVector* SV1 = new ScreamVector[10];
  //ScreamVector** SV2 = new ScreamVector*[10];

  this->gridPoints = new ScreamVector**[this->xGrid+1];
  if (gridPoints != NULL) {    for (int i = 0; i <= this->xGrid; i++) {
      gridPoints[i] = new ScreamVector*[this->yGrid+1];

      if (gridPoints[i] != NULL) {	for (int j = 0; j <= this->yGrid; j++) {
	  gridPoints[i][j] = new ScreamVector[this->zGrid+1];
      } }
  } }

  /* Populate 3D array. */
  for (int xx = 0; xx <= this->xGrid; xx++) {
    for (int yy = 0; yy <= this->yGrid; yy++) {
      for (int zz = 0; zz <= this->zGrid; zz++) {
	(this->gridPoints)[xx][yy][zz] = ScreamVector(this->min_x + spacing*xx, this->min_y + spacing*yy , this->min_z+spacing*zz);
      }
    }
  }
  //return &(this->gridPoints); // Pointer to 3D array?

}

std::set<ScreamVector> Box::getEnclosedPoints(SCREAM_ATOM* a) {
  std::set<ScreamVector> enclosedPoints;
  ScreamVector v(a->x[0], a->x[1], a->x[2]);
  enclosedPoints = this->getEnclosedPoints(v, a->vdw_r);
  
  return enclosedPoints;

}

std::set<ScreamVector> Box::getEnclosedPoints(ScreamAtomV& l) {
  std::set<ScreamVector> enclosedPoints;
  for (ScreamAtomVItr a = l.begin(); a != l.end(); a++) {
    std::set<ScreamVector> points = this->getEnclosedPoints(*a);
    enclosedPoints.insert(points.begin(), points.end());
  }
  return enclosedPoints;

}

std::set<ScreamVector> Box::getEnclosedPoints(ScreamVector v, double r) {
  /* This routine EFFICIENTLY finds points that are within specified radius of scream atom. */
  int rangeX_min, rangeY_min, rangeZ_min, rangeX_max, rangeY_max, rangeZ_max;
  std::set<ScreamVector> enclosedPoints;

  rangeX_min = int(floor( (v[0] - r - this->min_x)/spacing ));
  rangeY_min = int(floor( (v[1] - r - this->min_y)/spacing ));
  rangeZ_min = int(floor( (v[2] - r - this->min_z)/spacing ));

  rangeX_max = int(ceil( (v[0] + r - this->min_x)/spacing ));
  rangeY_max = int(ceil( (v[1] + r - this->min_y)/spacing ));
  rangeZ_max = int(ceil( (v[2] + r - this->min_z)/spacing ));

  for (int i = rangeX_min; i <= rangeX_max; i++ ) {
    for (int j = rangeY_min; j <= rangeY_max; j++) {
      for (int k = rangeZ_min; k <= rangeZ_max; k++) {
	if (this->_distanceSquared(v, this->gridPoints[i][j][k]) <= r*r) enclosedPoints.insert(this->gridPoints[i][j][k]);
      }
    }
  }

  return enclosedPoints;

}


double Box::_distanceSquared(ScreamVector& v1, ScreamVector& v2) {
  double x = v1[0] - v2[0];
  double y = v1[1] - v2[1];
  double z = v1[2] - v2[2];

  return x*x + y*y + z*z;
}

