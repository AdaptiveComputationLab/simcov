// LungModel
// Author: Akil Andrews

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
using namespace std;
namespace upcxx {
int rank_me() { return 0; }
}  // namespace upcxx
#include "CLI11.hpp"
#include "lung.hpp"

int main(int argc, char** argv) {
  CLI::App app("SIMCoV lung model generator");
  string output_dir = "lung_model_data";
  bool full_lung = false;
  int levels = 3;
  double scale = 2000;
  vector<int> dimensions{300, 300, 300};
  app.add_option("-d,--dim", dimensions, "Dimensions: x y z")
      ->delimiter(',')
      ->expected(3)
      ->capture_default_str();
  app.add_option("--levels", levels, "Number of levels to model a single lobe (0 for full lung)")
      ->check(CLI::Range(0, 24))
      ->capture_default_str();
  app.add_option("--scale", scale, "Scale factor for model")->capture_default_str();
  app.add_option("-o,--output", output_dir, "Output directory")->capture_default_str();
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    app.exit(e);
  }
  if (levels == 0) full_lung = true;
  Lung lung(full_lung,
    levels,
    scale,
    {dimensions[0], dimensions[1], dimensions[2]});
#ifndef COMPUTE_ONLY
  // Create output directory
  if (mkdir(output_dir.c_str(), S_IRWXU) == -1) {
    // could not create the directory
    if (errno == EEXIST) {
      cerr << "WARNING: Output directory " << output_dir
           << " already exists. May overwrite existing files\n";
    } else {
      cerr << "ERROR: Could not create output directory " << output_dir << ": " << strerror(errno)
           << endl;
      exit(1);
    }
  }
  // Write alveolus epithileal cells
  string fname = output_dir + "/alveolus.dat";
  auto fileno =
      open(fname.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) {
    cerr << "ERROR: Cannot create file " << fname << " " << strerror(errno) << "\n";
    exit(1);
  }
  int64_t bytes_written =
      write(fileno, reinterpret_cast<const char*>(&dimensions[0]), sizeof(int64_t) * dimensions.size());
  //TODO change simcov to input int64_t
  for (const int64_t& a : lung.getAlveoliEpiCellIds()) {
    bytes_written += write(fileno, reinterpret_cast<const char*>(&a), sizeof(int64_t));
  }
  if (bytes_written !=
      sizeof(int64_t) * (lung.getAlveoliEpiCellIds()).size() + sizeof(int64_t) * dimensions.size()) {
    cerr << "ERROR: Could not write all " << lung.getAlveoliEpiCellIds().size() << " bytes, wrote "
         << bytes_written << endl;
    exit(1);
  }
  close(fileno);
  // Write bronchiole epithileal cells
  fname = output_dir + "/bronchiole.dat";
  fileno = open(fname.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) {
    cerr << "ERROR: Cannot create file " << fname << " " << strerror(errno) << "\n";
    exit(1);
  }
  bytes_written =
      write(fileno, reinterpret_cast<const char*>(&dimensions[0]), sizeof(int64_t) * dimensions.size());
  for (const int64_t& a : lung.getAirwayEpiCellIds()) {
    bytes_written += write(fileno, reinterpret_cast<const char*>(&a), sizeof(int64_t));
  }
  if (bytes_written !=
      sizeof(int64_t) * (lung.getAirwayEpiCellIds()).size() + sizeof(int64_t) * dimensions.size()) {
    cerr << "ERROR: Could not write all " << lung.getAirwayEpiCellIds().size() << " bytes, wrote "
         << bytes_written << endl;
    exit(1);
  }
  close(fileno);
#endif

  return 0;
}
