#include "Viewer.h"

#include <cxxopts.h>

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Fluid Simulator Viewer");

    pbs::ViewerSettings settings;

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

    if (options.count("input") > 1) {
        std::cerr << "Provide at most one scene file!" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        return 0;
    } else if (options.count("input") == 0) {
        std::cerr << "Provide a scene file!" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
    }
    settings.filename = options["input"].as<std::vector<std::string>>().front();

    try {
        nanogui::init();
        std::unique_ptr<pbs::Viewer> screen(new pbs::Viewer(settings));
        nanogui::mainloop();
        nanogui::shutdown();
    } catch (const std::exception &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}
