#include "Renderer.h"

#include <tinyformat.h>
#include <cxxopts.h>
#include <json11.h>

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Fluid Simulator Renderer");

    pbs::RendererSettings settings;
    std::string renderMode = pbs::RendererSettings::renderModeToString(settings.renderMode);

    options.add_options()
    ("h,help", "Print help")
    ("width", tfm::format("Display Width (default: %d)", settings.width), cxxopts::value<int>(settings.width), "N")
    ("height", tfm::format("Display Height (default: %d)", settings.height), cxxopts::value<int>(settings.height), "N")
    ("framerate", tfm::format("Frame Rate (default: %d)", settings.framerate), cxxopts::value<int>(settings.framerate), "")
    ("skipToFrame", tfm::format("Skip to frame (default: %d)", settings.skipToFrame), cxxopts::value<int>(settings.skipToFrame), "")
    ("renderMode", tfm::format("Render Mode (default: %s)", renderMode), cxxopts::value<std::string>(renderMode), "[particles|mesh]")
    ("animation", "Animation Scene", cxxopts::value<std::string>(settings.filenameAnimation), "scene file")
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

    if (options.count("input") != 1) {
        std::cerr << "Provide a scene file!" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (options.count("animation") != 1) {
        std::cerr << "Provide an animation file!" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        return 0;
    }

    settings.filenameScene = options["input"].as<std::vector<std::string>>().front();
    settings.renderMode = pbs::RendererSettings::stringToRenderMode(renderMode);

    try {
        nanogui::init();
        std::unique_ptr<pbs::Renderer> screen(new pbs::Renderer(settings));
        nanogui::mainloop();
        nanogui::shutdown();
    } catch (const std::exception &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}