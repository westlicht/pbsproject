#include <tinyformat.h>
#include <cxxopts.h>

#include "core/FileUtils.h"
#include "core/Timer.h"

#include "sim/Scene.h"
#include "sim/SPH.h"
#include "sim/Cache.h"

#include "geometry/ParticleMesher.h"

using namespace pbs;

static inline MatrixXf toMatrix(const std::vector<Vector3f> &data) {
    MatrixXf result;
    result.resize(3, data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        result.col(i) = data[i];
    }
    return std::move(result);
}

struct Settings {
    std::string task = "mesh";
    int startFrame = -1;
    int endFrame = -1;
};

void run(const std::string &filename, const Settings &settings) {
    DBG("Loading scene from '%s' ...", filename);
    Scene scene = Scene::load(filename);
    DBG("%s", scene.toString());
    SPH sph(scene);
    std::string cachePath = FileUtils::splitExtension(filename).first + ".cache";
    Cache cache(cachePath);

    if (settings.task == "mesh") {
        if (!cache.valid()) {
            throw Exception("Cache is invalid!");
        }

        Timer timer;

        int startFrame = std::max(settings.startFrame, 0);
        int endFrame = std::min(settings.endFrame + 1, cache.frameCount());

        for (int frame = startFrame; frame < endFrame; ++frame) {
            DBG("processing frame %d ...", frame);

            cache.setFrame(frame);

            std::vector<Vector3f> particles;
            cache.readParticles(particles);

            ParticleMesher::Parameters params;
            params.particleRadius = sph.parameters().particleRadius;
            params.particleDiameter = sph.parameters().particleDiameter;
            params.kernelRadius = sph.parameters().kernelRadius;
            params.kernelSupportParticles = sph.parameters().kernelSupportParticles;
            params.particleMass = sph.parameters().particleMass;
            params.restDensity = sph.parameters().restDensity;

            MatrixXf positions = toMatrix(particles);
            Box3f bounds = sph.bounds().expanded(sph.bounds().extents() * 0.05f); // expand bounds by 5% of diagonal
            Vector3f extents = bounds.extents();

            Vector3i cells(
                int(std::ceil(extents.x() / params.particleRadius)),
                int(std::ceil(extents.y() / params.particleRadius)),
                int(std::ceil(extents.z() / params.particleRadius))
            );

            bool anisotropic = false;
            params.isoLevel = anisotropic ? 0.2f : 0.75f;
            Mesh mesh = anisotropic ? ParticleMesher::createMeshAnisotropic(positions, bounds, cells, params)
                                    : ParticleMesher::createMeshIsotropic(positions, bounds, cells, params);
            cache.writeMesh(mesh);

            // Show time estimate
            float elapsed = timer.elapsed();
            float progress = float(frame - startFrame + 1) / (endFrame - startFrame);
            float eta = progress != 0.f ? elapsed / progress - elapsed : 0.f;
            DBG("Progress %.1f%% (Elapsed: %s ETA: %s)", progress * 100.f, timeString(elapsed), timeString(eta));

        }
        cache.commit();

    } else {
        throw Exception("Unknown task '%s'", settings.task);
    }
}

int main(int argc, char *argv[]) {

    cxxopts::Options options(argv[0], " - Fluid Simulator Processor");

    Settings settings;

    options.add_options()
    ("h,help", "Print help")
    ("task", tfm::format("Processing task (default: %s)", settings.task), cxxopts::value<std::string>(settings.task), "")
    ("startFrame", tfm::format("First frame to process (default: %d)", settings.startFrame), cxxopts::value<int>(settings.startFrame), "")
    ("endFrame", tfm::format("Last frame to process (default: %d)", settings.endFrame), cxxopts::value<int>(settings.endFrame), "")
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
    std::string filename = options["input"].as<std::vector<std::string>>().front();

    try {
        run(filename, settings);
    } catch (const std::exception &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}