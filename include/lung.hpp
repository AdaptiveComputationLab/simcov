#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "upcxx_utils.hpp"
#include "utils.hpp"
#include "options.hpp"

extern shared_ptr<Options> _options;

struct Int3D {

  int64_t x, y, z;

  Int3D() : x(0), y(0), z(0) {}

  Int3D(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}

  Int3D & operator=(const Int3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }

};

struct Double3D {

  double x, y, z;

  Double3D() : x(0.0), y(0.0), z(0.0) {}

  Double3D(double x, double y, double z) : x(x), y(y), z(z) {}

  Double3D & operator=(const Double3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }

};

struct Level {

  int64_t L = 0, d = 0, count = 0;
  double bAngle = 0.0, gAngle = 0.0;
  Int3D centroid;
  Double3D direction;

};

class Lung {

 public:

  Lung(shared_ptr<Random> rnd_gen) {
    this->rnd_gen = rnd_gen;
    gridSize.x = _options->dimensions[0];
    gridSize.y = _options->dimensions[1];
    gridSize.z = _options->dimensions[2];
    if (_options->generate_lung_model) {
      loadEstimatedParameters();
      // Draw upper right lobe 24 gen
      constructLobe(0, 15);
      // Draw middle right lobe 24 gen
      constructLobe(24, 15);
      // Draw lower right lobe 26 gen
      constructLobe(48, 17);
      // Draw upper left lobe 24 gen
      constructLobe(74, 15);
      // Draw lower left lobe 25 gen
      constructLobe(98, 16);
    } else {
      loadData();
      //TODO Draw all epitheleal cells
    }
  }

  ~Lung() {}

  const std::set<int64_t> &getAirwayEpiCellIds() { return airwayEpiCellPositions1D; }

  const std::set<int64_t> &getAlveoliEpiCellIds() { return alveoliEpiCellPositions1D; }

  const std::vector<Int3D> &getEpiLocations() { return positions; }

 private:

   double scale = 10;// 1 / 5e-3; // Convert cm to um
   std::vector<Int3D> leafs;
   std::vector<Level> levels;
   std::set<int64_t> alveoliEpiCellPositions1D;
   std::set<int64_t> airwayEpiCellPositions1D;
   std::vector<Int3D> positions;
   Int3D gridSize;
   shared_ptr<Random> rnd_gen;

   Int3D rotate(const Int3D & vec,
     const Double3D & axis,
     const double & angle) {
    /**
    * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
    */
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    double x = vec.x, y = vec.y, z = vec.z;
    //double u = 0, v = 1, w = 0;
    double u = axis.x, v = axis.y, w = axis.z;
    double dotProduct = (x * u) + (y * v) + (z * w);
    double newX = u * dotProduct * (1 - cosAngle) + x * cosAngle
      + (-w * y + v * z) * sinAngle;
    double newY = v * dotProduct * (1 - cosAngle) + y * cosAngle
      + (w * x - u * z) * sinAngle;
    double newZ = w * dotProduct * (1 - cosAngle) + z * cosAngle
      + (-v * x + u * y) * sinAngle;
    Int3D rval((int64_t) round(newX),
      (int64_t) round(newY),
      (int64_t) round(newZ));
    return rval;
   }

   void constructLobe(int startIndex, int generations) {
     // Draw root
     Level lvl = levels.at(startIndex);
     Int3D child = constructSegment({gridSize.x/2, gridSize.y/2, 0},
       lvl,
       0.0);
     // Recursively build tree
     construct(child, 1, startIndex + 1, generations - 1, lvl.bAngle);
     // Build alveolus structures
     if (scale == 1 / 5e-3) {
       for (Int3D leaf : leafs) {//TODO Int3D leaf(0, 0, 0);
         Int3D left(leaf.x + 2, leaf.y + 2, leaf.z);
         constructAlveoli(left);
         Int3D right(leaf.x + 7, leaf.y + 7, leaf.z);
         constructAlveoli(right);
       }
       SLOG("Number of alveoli ", 2 * leafs.size(), "\n");
       leafs.clear();
     }
   }

   void construct(const Int3D & root,
     int iteration,
     int index,
     int end,
     double previousBranchAngle) {
       if (iteration > end) {
         leafs.push_back({.x = root.x, .y = root.y, .z = root.z});
         return;
       }
       double rotateZ = (iteration < 7) ? 0.0 : rnd_gen->get(0,180)*M_PI/180;
       // Draw left child bronchiole
       Level lvl = levels.at(index);
       lvl.bAngle = previousBranchAngle - lvl.bAngle;
       Int3D child = constructSegment(root, lvl, rotateZ);
       construct(child, iteration + 1, index + 1, end, lvl.bAngle);
       // Draw right child bronchiole
       lvl = levels.at(index);
       lvl.bAngle = previousBranchAngle + lvl.bAngle;
       if (lvl.count > 1) {
         lvl.gAngle = -lvl.gAngle;
         child = constructSegment(root, lvl, rotateZ);
         construct(child, iteration + 1, index + 1, end, lvl.bAngle);
       }
     }

   Int3D constructSegment(const Int3D &root,
     const Level &level,
     double rotateZ) {
     std::vector<Int3D> positions;
     // Build cylinder at origin along y-axis
     int64_t radius = level.d / 2;
     double az = 0;
     double inc = M_PI / 2;
     for (int64_t z = 0; z <= level.L; z++) {
       for (az = 0; az < 2 * M_PI; az += M_PI / 180) {
         int64_t x = (int64_t) round(radius * sin(inc) * cos(az));
         int64_t y = (int64_t) round(radius * sin(inc) * sin(az));
         positions.push_back({.x = x, .y = y, .z = z});
       }
     }
     // Treat as positional vectors and apply rotations and translate
     for (int64_t i = 0; i < positions.size(); i++) {
       // Rotate y
       Int3D newPosition0 = rotate(positions.at(i),
        {.x = 0.0, .y = 1.0, .z = 0.0},
        level.bAngle);
       // Rotate z
       Int3D newPosition1 = rotate(newPosition0,
        {.x = 0.0, .y = 0.0, .z = 1.0},
        rotateZ);
       newPosition1.x += root.x;
       newPosition1.y += root.y;
       newPosition1.z += root.z;
       positions.at(i) = newPosition1;
       // Verify new location is within grid
       if ((0 <= newPosition1.x && newPosition1.x < gridSize.x)
        && (0 <= newPosition1.y && newPosition1.y < gridSize.y)
        && (0 <= newPosition1.z && newPosition1.z < gridSize.z)) {
          airwayEpiCellPositions1D.insert(newPosition1.x
            + newPosition1.y * gridSize.x
            + newPosition1.z * gridSize.x * gridSize.y);
       }
     }
     // Return root for next generation
     Int3D base0 = rotate({.x = 0, .y = 0, .z = level.L},
       {.x = 0.0, .y = 1.0, .z = 0.0},
       level.bAngle);
     Int3D base1 = rotate(base0,
       {.x = 0.0, .y = 0.0, .z = 1.0},
       rotateZ);
     Int3D rval(base1.x + root.x,
       base1.y + root.y,
       base1.z + root.z);
     return rval;
   }

   void constructAlveoli(const Int3D & pos) {
     // Cube dimensions for dim = 5 epi cells across are x,y,z = [-2 2]
     int dim = 5, idim = 2;
     for (int x = -idim; x <= idim; x++) {
       for (int y = -idim; y <= idim; y++) {
         for (int z = 0; z < dim; z++) {
           // Cells in the two x planes
           if (x == -idim || x == idim) {
             addPosition(x, y, z, pos);
           }
           // Cells in the two y planes
           if (y == -idim || y == idim) {
             addPosition(x, y, z, pos);
           }
           // Cells in the one z plane at bottom of alveoli
           if (z == dim - 1) {
             addPosition(x, y, z, pos);
           }
         }
       }
     }
     // TODO rotate
   }

   void addPosition(int x, int y, int z, const Int3D & pos) {
     positions.push_back({.x = x + pos.x,
       .y = y + pos.y,
       .z = z + pos.z});
     alveoliEpiCellPositions1D.insert((x + pos.x)
      + (y + pos.y) * gridSize.x
      + (z + pos.z) * gridSize.x * gridSize.y);
   }

   void loadEstimatedParameters() {
     // Yeh et al 1980
     std::ifstream file;
     file.open(_options->lung_parameters_filepath);
     std::string line;
     if (file.is_open()) {
       while (std::getline(file, line)) {
         std::stringstream lstream(line);
         Level e;
         int tempI;
         double tempD;
         lstream >> tempI;
         lstream >> e.count;
         lstream >> tempD;
         e.L = (int64_t) round(scale * tempD);
         lstream >> tempD;
         e.d = (int64_t) round(scale * tempD);
         lstream >> tempI;
         e.bAngle = tempI * M_PI / 180;
         lstream >> tempI;
         e.gAngle = tempI * M_PI / 180;
         if (levels.size() > 73) { // Negate for left lobe data
           e.bAngle *= -1.0;
         }
         levels.push_back(e);
       }
       SLOG("Loaded ", levels.size(), " estimated levels\n");
       file.close();
     } else {
       SDIE("Failed to open file ", _options->lung_parameters_filepath);
     }
   }

   void loadData() {
     //TODO Load cell data from file
   }

};
