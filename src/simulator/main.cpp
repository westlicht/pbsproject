#include "Simulator.h"

#include <tinyformat.h>
#include <cxxopts.h>
#include <json11.h>

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Fluid Simulator");

    pbs::SimulatorSettings settings;
    std::string renderMode = pbs::SimulatorSettings::renderModeToString(settings.renderMode);

    options.add_options()
    ("h,help", "Print help")
    ("width", tfm::format("Display Width (default: %d)", settings.width), cxxopts::value<int>(settings.width), "N")
    ("height", tfm::format("Display Height (default: %d)", settings.height), cxxopts::value<int>(settings.height), "N")
    ("duration", tfm::format("Duration (default: %.1f s)", settings.duration), cxxopts::value<float>(settings.duration), "s")
    ("timescale", tfm::format("Time Scaling (default: %.1f)", settings.timescale), cxxopts::value<float>(settings.timescale), "")
    ("framerate", tfm::format("Frame Rate (default: %d)", settings.framerate), cxxopts::value<int>(settings.framerate), "")
    ("renderMode", tfm::format("Render Mode (default: %s)", renderMode), cxxopts::value<std::string>(renderMode), "[particles|mesh]")
    ("cacheParticles", tfm::format("Write particle cache (default: %d)", settings.cacheParticles), cxxopts::value<bool>(settings.cacheParticles), "")
    ("cacheMesh", tfm::format("Write mesh cache (default: %d)", settings.cacheMesh), cxxopts::value<bool>(settings.cacheMesh), "")
    ("D,define", "Parameter definition", cxxopts::value<std::vector<std::string>>(), "name=value")
    ("tag", "Tag to append to output files", cxxopts::value<std::string>(settings.tag), "")
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

    // Parse definitions
    json11::Json::object definitionValues;
    if (options.count("D")) {
        const auto &definitions = options["D"].as<std::vector<std::string>>();
        for (const auto &definition : definitions) {
            auto tokens = pbs::tokenize(definition, "=");
            if (tokens.size() == 2) {
                std::string error;
                auto json = json11::Json::parse(tfm::format("{\"%s\":%s}", tokens[0], tokens[1]), error);
                if (json == json11::Json()) {
                    json = json11::Json::parse(tfm::format("{\"%s\":\"%s\"}", tokens[0], tokens[1]), error);
                    if (json == json11::Json()) {
                        std::cerr << tfm::format("Invalid definition %s (JSON error: %d)", definition, error) << std::endl;
                        return -1;
                    }
                }
                for (auto kv : json.object_items()) {
                    definitionValues.emplace(kv.first, kv.second);
                }
            } else {
                std::cerr << tfm::format("Invalid definition '%s' (needs to be in the form of 'name=value')", definition) << std::endl;
                return -1;
            }
        }
    }
    settings.sceneSettings = json11::Json(definitionValues);

    if (options.count("input") != 1) {
        std::cerr << "Provide a scene file!" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        return 0;
    }
    settings.filename = options["input"].as<std::vector<std::string>>().front();
    settings.renderMode = pbs::SimulatorSettings::stringToRenderMode(renderMode);

    try {
        nanogui::init();
        std::unique_ptr<pbs::Simulator> screen(new pbs::Simulator(settings));
        nanogui::mainloop();
        nanogui::shutdown();
    } catch (const std::exception &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}