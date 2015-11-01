#include "Viewer2d.h"
#include "Viewer3d.h"
//#include "Test.h"

#include <cxxopts.h>

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Melting Ice Simulator");

    bool run2d = false;
    bool run3d = false;

    options.add_options()
    ("h,help", "Print help")
    ("2d", "Run 2D mode (default)", cxxopts::value<bool>(run2d))
    ("3d", "Run 3D mode", cxxopts::value<bool>(run3d))
    ("input", "Input files", cxxopts::value<std::vector<std::string>>())
    ;

    options.parse_positional("input");

    // Parse command line arguments
    try {
        options.parse(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << "Error during command line parsing: " << e.what() << std::endl;
        return -1;
    }

    // Show help if requested
    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    //pbs::test::obj();
    //pbs::test::marchingCubes();

    nanogui::init();

    int mode = 0;
    if (run3d) mode = 1;

    switch (mode) {
    case 0: { // 2d mode
        std::unique_ptr<Viewer2d> screen(new Viewer2d());
        nanogui::mainloop();
        break;
    }
    case 1: { // 3d mode
        std::unique_ptr<Viewer3d> screen(new Viewer3d());
        nanogui::mainloop();
        break;
    }
    }

    nanogui::shutdown();

    return 0;
}