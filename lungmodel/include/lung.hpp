#pragma once

#include <random>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <memory>
#include <fcntl.h>
#include <inttypes.h>

class Random {
 private:
  std::mt19937_64 generator;

  double get_prob(double max_val = 1.0) {
    return std::uniform_real_distribution<>(0, max_val)(generator);
  }

 public:
  Random(unsigned seed)
      : generator(seed) {}

  int get(int64_t begin, int64_t end) {
    return std::uniform_int_distribution<int64_t>(begin, end - 1)(generator);
  }

  bool trial_success(double thres) {
    assert(thres >= 0);
    if (thres > 1) return true;
    if (thres == 0) return false;
    return (get_prob() <= thres);
  }

  int get_normal(std::vector<int> dist_params) {
    return (int)std::normal_distribution<float>(dist_params[0], dist_params[1])(generator);
  }

  int get_poisson(int avg) { return (int)std::poisson_distribution<int>(avg)(generator); }
};

struct Int3D {
  int x, y, z;

  Int3D()
      : x(0)
      , y(0)
      , z(0) {}

  Int3D(int x, int y, int z)
      : x(x)
      , y(y)
      , z(z) {}

  Int3D &operator=(const Int3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }
};

struct Double3D {
  double x, y, z;

  Double3D()
      : x(0.0)
      , y(0.0)
      , z(0.0) {}

  Double3D(double x, double y, double z)
      : x(x)
      , y(y)
      , z(z) {}

  Double3D &operator=(const Double3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }
};

struct Level {
  int L = 0, r = 0, count = 0;
  double bAngle = 0.0, gAngle = 0.0;
  Int3D centroid;
  Double3D direction;
};

class Lung {

 public:

  Lung(bool isFullModel, int levels_to_model, double scale, Int3D gridSize)
      : isFullModel(isFullModel)
      , scale(scale)
      , levels_to_model(levels_to_model)
      , gridSize(gridSize) {
    this->rnd_gen = std::make_shared<Random>(753695190);
    if (isFullModel) {
      loadEstimatedParameters();
      // Draw upper right lobe 24 gen
      constructLobe(0, 24);
      // Draw middle right lobe 24 gen
      constructLobe(24, 24);
      // Draw lower right lobe 26 gen
      constructLobe(48, 26);
      // Draw upper left lobe 24 gen
      constructLobe(74, 24);
      // Draw lower left lobe 25 gen
      constructLobe(98, 25);
    } else { // Draw partial upper right lobe
      loadEstimatedParameters();
      constructLobe(24 - levels_to_model, levels_to_model);
    }
    std::printf("Number of alveoli " "%" PRId64 " epithileal cells " "%" PRId64 "\n", numAlveoli, numAlveoliCells);
    std::printf("Number of airways " "%" PRId64 " epithileal cells " "%" PRId64 "\n", numAirways, numAirwayCells);
    std::printf("Number of OUT OF BOUNDS epithileal cells " "%" PRId64 "\n", numOutOfBoundsCells);
  }

  ~Lung() {}

  const std::set<int64_t> &getAirwayEpiCellIds() { return airwayEpiCellPositions1D; }

  const std::set<int64_t> &getAlveoliEpiCellIds() { return alveoliEpiCellPositions1D; }

 private:

  bool isFullModel = false;
  /**
    scale = 10    => 1 unit = 1mm, 1cm = 10^-2m = 10^4 um, 10^4/10^3 = 10 units
    scale = 200   => 1 unit = 50um, 1cm = 10^-2m = 10^4 umm 10^4/50 = 200 units
    scale = 2000  => 1 unit = 5um, 1cm = 10^-2m = 10^4 um, 10^4/5 = 2000 units
    */
  double scale = 2000;
  int levels_to_model = 3;
  int64_t numAlveoli = 0, numAlveoliCells = 0;
  int64_t numAirways = 0, numAirwayCells = 0;
  int64_t numOutOfBoundsCells = 0;
  std::vector<Level> levels;
  std::set<int64_t> alveoliEpiCellPositions1D;
  std::set<int64_t> airwayEpiCellPositions1D;
  Int3D gridSize = {300, 300, 300};
  std::shared_ptr<Random> rnd_gen;
  double cylinderRadialIncrement = M_PI / 2;

  Int3D rotate(const Int3D &vec, const Double3D &axis, const double &angle) {
    /**
     * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
     */
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    double x = vec.x, y = vec.y, z = vec.z;
    // double u = 0, v = 1, w = 0;
    double u = axis.x, v = axis.y, w = axis.z;
    double dotProduct = (x * u) + (y * v) + (z * w);
    double newX = u * dotProduct * (1 - cosAngle) + x * cosAngle + (-w * y + v * z) * sinAngle;
    double newY = v * dotProduct * (1 - cosAngle) + y * cosAngle + (w * x - u * z) * sinAngle;
    double newZ = w * dotProduct * (1 - cosAngle) + z * cosAngle + (-v * x + u * y) * sinAngle;
    Int3D rval((int)round(newX), (int)round(newY), (int)round(newZ));
    return rval;
  }

  void constructLobe(int startIndex, int generations) {
    // Draw root and recursively build tree
    Level lvl = levels.at(startIndex);
    if (isFullModel) {
      Int3D child = constructSegment({gridSize.x / 2, gridSize.y / 2, 0}, lvl, 0.0, false);
      construct(child, 1, startIndex + 1, generations - 1, lvl.bAngle, 0.0);
    } else {
      Int3D child = constructSegment({0, gridSize.y / 4, gridSize.z / 2}, lvl, 0.0, false);
      construct(child, 1, startIndex + 1, generations - 1, lvl.bAngle, 0.0);
    }
  }

  void construct(const Int3D &root, int iteration, int index, int end, double previousBranchAngle,
                 double previousRotAngle) {
    if (iteration > end) {
      return;
    }
    bool isTerminal = (iteration == end) ? true : false;
    // Uniform randomly rotate each branch
    double rotateZ = getUniformRandomRotation(iteration);
    rotateZ = previousRotAngle - rotateZ;
    // Draw left child bronchiole
    Level lvl = levels.at(index);
    lvl.bAngle = previousBranchAngle - lvl.bAngle;
    Int3D child = constructSegment(root, lvl, rotateZ, isTerminal);
    construct(child, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
    // Uniform randomly rotate each branch
    rotateZ = getUniformRandomRotation(iteration);
    rotateZ = previousRotAngle + rotateZ;
    // Draw right child bronchiole
    lvl = levels.at(index);
    lvl.bAngle = previousBranchAngle + lvl.bAngle;
    if (lvl.count > 1) {
      lvl.gAngle = -lvl.gAngle;
      child = constructSegment(root, lvl, rotateZ, isTerminal);
      construct(child, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
    }
  }

  double getUniformRandomRotation(int iteration) {
    //TODO pull from distribution
    double rval = 0.0;
    if ((isFullModel && iteration >= 7) || (!isFullModel && iteration >= 2)) {
      rval = rnd_gen->get(0, 180) * M_PI / 180;
    }
    return rval;
  }

  Int3D constructSegment(const Int3D &root, const Level &level, double rotateZ, bool isTerminal) {
    // Build cylinder at origin along y-axis
    for (int z = 0; z <= level.L; z++) {
      for (double az = 0; az < 2 * M_PI; az += M_PI / 180) {
#ifndef COMPUTE_ONLY
        int x = (int)round(level.r * sin(cylinderRadialIncrement) * cos(az));
        int y = (int)round(level.r * sin(cylinderRadialIncrement) * sin(az));
        Int3D position({.x = x, .y = y, .z = z});
        // Treat as positional vectors and apply rotations and translate y
        Int3D newPosition0 = rotate(position, {.x = 0.0, .y = 1.0, .z = 0.0}, level.bAngle);
        // Rotate z
        Int3D newPosition1 = rotate(newPosition0, {.x = 0.0, .y = 0.0, .z = 1.0}, rotateZ);
        newPosition1.x += root.x;
        newPosition1.y += root.y;
        newPosition1.z += root.z;
        // Set and verify new location is within grid
        if ((0 <= newPosition1.x && newPosition1.x < gridSize.x) &&
            (0 <= newPosition1.y && newPosition1.y < gridSize.y) &&
            (0 <= newPosition1.z && newPosition1.z < gridSize.z)) {
              airwayEpiCellPositions1D.insert(newPosition1.x + newPosition1.y
                * gridSize.x + newPosition1.z * gridSize.x * gridSize.y);
              numAirwayCells++;
        } else {
          numOutOfBoundsCells++;
        }
#else
        numAirwayCells++;
#endif
      }
    }
    // Create root for next generation
    Int3D base0 =
        rotate({.x = 0, .y = 0, .z = level.L}, {.x = 0.0, .y = 1.0, .z = 0.0}, level.bAngle);
    Int3D base1 = rotate(base0, {.x = 0.0, .y = 0.0, .z = 1.0}, rotateZ);
    Int3D rval(base1.x + root.x, base1.y + root.y, base1.z + root.z);
    numAirways++;
    // Draw alveolus
    if (isTerminal && scale == 2000) {
      constructAlveoli(rval, level.bAngle, rotateZ);
    }
    return rval;
  }

  void constructAlveoli(const Int3D &pos, double bAngle, double rotateZ) {
    // Single alveolar volume 200um x 200um x 200um, ? et al ?
    // 40 x 40 x 40 units, [-20, 20]
    int idim = 20;
    for (int x = -idim; x <= idim; x++) {
      for (int y = -idim; y <= idim; y++) {
        for (int z = 0; z < (2 * idim); z++) {
          // Cells in the two x planes
          if (x == -idim || x == idim) {
            addPosition(x, y, z, pos, bAngle, rotateZ);
          }
          // Cells in the two y planes
          if (y == -idim || y == idim) {
            addPosition(x, y, z, pos, bAngle, rotateZ);
          }
          // Cells in the one z plane at bottom of alveoli
          if (z == (2 * idim) - 1) {
            addPosition(x, y, z, pos, bAngle, rotateZ);
          }
        }
      }
    }
    numAlveoli++;
  }

  void addPosition(int x, int y, int z, const Int3D &pos, double bAngle, double rotateZ) {
#ifndef COMPUTE_ONLY
    // Rotate y
    Int3D newPosition0 = rotate({.x = x, .y = y, .z = z}, {.x = 0.0, .y = 1.0, .z = 0.0}, bAngle);
    // Rotate z
    Int3D newPosition1 = rotate(newPosition0, {.x = 0.0, .y = 0.0, .z = 1.0}, rotateZ);
    // Verify new location is within grid
    newPosition1.x += pos.x;
    newPosition1.y += pos.y;
    newPosition1.z += pos.z;
    if ((0 <= newPosition1.x && newPosition1.x < gridSize.x) &&
        (0 <= newPosition1.y && newPosition1.y < gridSize.y) &&
        (0 <= newPosition1.z && newPosition1.z < gridSize.z)) {
          alveoliEpiCellPositions1D.insert(newPosition1.x + newPosition1.y
            * gridSize.x + newPosition1.z * gridSize.x * gridSize.y);
          numAlveoliCells++;
    } else {
      numOutOfBoundsCells++;
    }
#else
    numAlveoliCells++;
#endif
  }

  void loadEstimatedParameters() {
    // Yeh et al 1980
    std::ifstream file;
    file.open("table.txt");
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
        e.L = (int)round(scale * tempD);
        lstream >> tempD;
        e.r = (int)(round(scale * tempD) / 2);
        lstream >> tempI;
        e.bAngle = tempI * M_PI / 180;
        lstream >> tempI;
        e.gAngle = tempI * M_PI / 180;
        if (levels.size() > 73) {  // Negate for left lobe data
          e.bAngle *= -1.0;
        }
        levels.push_back(e);
      }
      std::printf("Loaded estimated levels %ld\n", levels.size());
      file.close();
    } else {
      std::printf("Failed to open table.txt\n");
    }
  }
};
