#include "Viewer.h"

#include <cxxopts.h>

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Fluid Simulator Viewer");

    options.add_options()
    ("h,help", "Print help")
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

    nanogui::init();
    std::unique_ptr<pbs::Viewer> screen(new pbs::Viewer());
    nanogui::mainloop();
    nanogui::shutdown();

    return 0;
}
