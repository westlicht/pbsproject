#include "Renderer.h"

#include "core/FileUtils.h"
#include "core/StringUtils.h"
#include "sim/Scene.h"
#include "sim/SPH.h"
#include "geometry/PlyWriter.h"

#include <exec-stream.h>
#include <tinydir.h>

namespace pbs {

Renderer::Renderer(const RendererSettings &settings) :
    nanogui::Screen(nanogui::Vector2i(1280, 720), "Fluid Simulator Renderer"),
    _settings(settings),
    _engine(mNVGContext, mSize, Vector2i(settings.width, settings.height))
{
    setVisible(true);

    _engine.loadScene(settings.filenameScene);
    _engine.viewOptions().showDebug = false;

    _animationScene = Scene::load(settings.filenameAnimation);
    if (_animationScene.cameraKeyframes.size() != 2) {
        throw Exception("Provide exactly two camera keyframes!");
    }
    _startCamera = _animationScene.cameraKeyframes[0];
    _endCamera = _animationScene.cameraKeyframes[1];
    _frameFirst = _startCamera.frame;
    _frameCount = _endCamera.frame - _startCamera.frame;
    _frameIndex = 0;
    if (settings.skipToFrame >= 0) {
        _frameIndex = settings.skipToFrame;
    }

    _animationName = FileUtils::splitExtension(_settings.filenameAnimation).first;
    _basename = FileUtils::splitExtension(_settings.filenameScene).first + "-" + _animationName;

    setupEmptyDirectory("images");

    if (settings.renderMode == RendererSettings::SLG) {
        FileUtils::createDir(_animationName + "-images");
    }

    _timer.reset();
}

Renderer::~Renderer() {
}

bool Renderer::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (!Screen::keyboardEvent(key, scancode, action, modifiers)) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_ESCAPE:
                terminate();
                break;
            }
        }
    }
    return true;
}

void Renderer::drawContents() {

    _engine.viewOptions().showDomain = true;
    _engine.viewOptions().showBoundaryParticles = false;
    _engine.viewOptions().showBoundaryMeshes = true;
    _engine.viewOptions().showDebug = false;
    _engine.viewOptions().showCache = true;

    switch (_settings.renderMode) {
    case RendererSettings::Particles:
        _engine.viewOptions().showFluidParticles = true;
        _engine.viewOptions().showFluidMesh = false;
        break;
    case RendererSettings::Mesh:
    case RendererSettings::SLG:
        _engine.viewOptions().showFluidParticles = false;
        _engine.viewOptions().showFluidMesh = true;
        break;
    }

    // Interpolate camera
    {
        float t = float(_frameIndex) / _frameCount;
        _engine.camera().setPosition(lerp(t, _startCamera.position, _endCamera.position));
        _engine.camera().setTarget(lerp(t, _startCamera.target, _endCamera.target));
        _engine.camera().setUp(lerp(t, _startCamera.up, _endCamera.up).normalized());
        _engine.camera().setFov(lerp(t, _startCamera.fov, _endCamera.fov));
        _engine.camera().setNear(lerp(t, _startCamera.near, _endCamera.near));
        _engine.camera().setFar(lerp(t, _startCamera.far, _endCamera.far));
    }

    _engine.setCacheFrame(_frameFirst + _frameIndex);

    _engine.render();

    switch (_settings.renderMode) {
    case RendererSettings::Particles:
    case RendererSettings::Mesh:
        _engine.renderToPNG(tfm::format("images/frame%04d.png", _frameIndex));
        break;
    case RendererSettings::SLG:
        renderSLG();
        break;
    }

    ++_frameIndex;

    if (_frameIndex >= _frameCount) {
        terminate();
    }

    // Show time estimate
    float elapsed = _timer.elapsed();
    float progress = float(_frameIndex) / _frameCount;
    float eta = progress != 0.f ? elapsed / progress - elapsed : 0.f;
    std::cout << tfm::format("Progress %.1f%% (Elapsed: %s ETA: %s)", progress * 100.f, timeString(elapsed), timeString(eta)) << std::flush;

    glfwPostEmptyEvent();
}

void Renderer::initialize() {

}

void Renderer::terminate() {
    createVideo();
    setVisible(false);
}

void Renderer::createVideo() {
    // FFMPEG candidate locations
    std::vector<std::string> candidates = {
        "/opt/local/bin/ffmpeg",
        "/usr/bin/avconv"
    };

    std::string ffmpeg;
    for (const auto &candidate : candidates) {
        if (filesystem::path(candidate).exists()) {
            ffmpeg = candidate;
            break;
        }
    }
    if (ffmpeg.empty()) {
        std::cerr << "Cannot find FFMPEG to encode video!" << std::endl;
        return;
    }

    std::string filename = _basename + ".mp4";

    std::cout << std::endl;
    std::cout << tfm::format("Encoding video '%s' ...", filename) << std::endl;

    try {
        std::string arguments = tfm::format("-y -framerate %d -i images/frame%%04d.png -pix_fmt yuv420p -vcodec libx264 -r %d -preset slow -crf 10 %s", _settings.framerate, _settings.framerate, filename);
        if (_settings.renderMode == RendererSettings::SLG) {
            arguments = tfm::format("-y -framerate %d -i %s-images/frame-%%04d.png -pix_fmt yuv420p -vcodec libx264 -r %d -preset slow -crf 10 %s", _settings.framerate, _animationName, _settings.framerate, filename);
        }
        exec_stream_t es;
        // Set 5 minutes timeout
        es.set_wait_timeout(exec_stream_t::s_out | exec_stream_t::s_err | exec_stream_t::s_child, 300*1000);
        es.start(ffmpeg, arguments);

        std::string s;
        while (std::getline(es.err(), s).good()) {
            std::cout << s << std::endl;
        }

        if (!es.close()) {
            std::cerr << "FFMPEG timeout!" << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << "exec-stream error: " << e.what() << "\n";
    }
}

void Renderer::renderSLG() {

    FileUtils::createDir("slg");
    std::string currentDir = FileUtils::getCurrentDir();
    FileUtils::changeCurrentDir("slg");

    std::string imageFilename = tfm::format("../%s-images/frame-%04d.png", _animationName, _frameIndex);

    int width = 1280;
    int height = 720;
    int spp = 64;

    {
        std::ofstream os("render.cfg");
        os << "opencl.platform.index = -1" << std::endl;
        os << "opencl.cpu.use = 0" << std::endl;
        os << "opencl.gpu.use = 1" << std::endl;
        os << "opencl.gpu.workgroup.size = 64" << std::endl;
        os << "scene.epsilon.min = 9.99999971718068536574719e-10" << std::endl;
        os << "scene.epsilon.max = 0.100000001490116119384766" << std::endl;
        os << "renderengine.type = PATHOCL" << std::endl;
        os << "accelerator.instances.enable = 0" << std::endl;
        os << "path.maxdepth = 12" << std::endl;
        os << tfm::format("film.width = %d", width) << std::endl;
        os << tfm::format("film.height = %d", height) << std::endl;
        os << tfm::format("image.filename = %s", imageFilename) << std::endl;
        os << "sampler.type = RANDOM" << std::endl;
        os << "film.filter.type = GAUSSIAN" << std::endl;
        os << "film.filter.xwidth = 1" << std::endl;
        os << "film.filter.ywidth = 1" << std::endl;
        os << "scene.file = scene.scn" << std::endl;
        os << "film.imagepipeline.0.type = TONEMAP_LINEAR" << std::endl;
        os << "film.imagepipeline.1.type = GAMMA_CORRECTION" << std::endl;
        os << "film.imagepipeline.1.value = 2.2" << std::endl;
        os << tfm::format("batch.haltspp = %d", spp) << std::endl;
    }

    {
        std::ofstream os("scene.scn");

        Vector3f origin = _engine.camera().position();
        Vector3f target = _engine.camera().target();
        Vector3f up = _engine.camera().up();
        float fov = _engine.camera().fov() * 1.6f;

        os << tfm::format("scene.camera.lookat.orig = %f %f %f", origin.x(), origin.y(), origin.z()) << std::endl;
        os << tfm::format("scene.camera.lookat.target = %f %f %f", target.x(), target.y(), target.z()) << std::endl;
        os << tfm::format("scene.camera.up = %f %f %f", up.x(), up.y(), up.z()) << std::endl;
        os << tfm::format("scene.camera.fieldofview = %f", fov) << std::endl;

#if 0
#scene.camera.screenwindow = -1 1 -1 1
scene.camera.cliphither = 0.0010000000474974513053894
scene.camera.clipyon = 1.00000001504746621987669e+30
scene.camera.lensradius = 0.00218750000931322574615479
scene.camera.focaldistance = 0.28457677364349365234375
scene.camera.horizontalstereo.enable = 0
scene.camera.horizontalstereo.oculusrift.barrelpostpro.enable = 0
scene.camera.autofocus.enable = 1
#endif

        Vector3f dir = Vector3f(1.f).normalized();
#if 0
        os << tfm::format("scene.sunlight.dir = %f %f %f", dir.x(), dir.y(), dir.z()) << std::endl;
        os << "scene.sunlight.turbidity = 2.2" << std::endl;
        os << "scene.sunlight.relsize = 1" << std::endl;
        os << "scene.sunlight.gain = 0.0001 0.0001 0.0001" << std::endl;
#endif
        os << tfm::format("scene.skylight.dir = %f %f %f", dir.x(), dir.y(), dir.z()) << std::endl;
        os << "scene.skylight.turbidity = 2.2" << std::endl;
        os << "scene.skylight.gain = 0.00007 0.00007 0.00007" << std::endl;

        os << "scene.materials.fluid.type = matte" << std::endl;
        os << "scene.materials.fluid.kd = 0.3 0.3 1.0" << std::endl;
        os << "scene.materials.boundary.type = matte" << std::endl;
        os << "scene.materials.boundary.kd = 0.5 0.5 0.5" << std::endl;
        os << "scene.materials.floor.type = matte" << std::endl;
        os << "scene.materials.floor.kd = 0.2 0.2 0.2" << std::endl;

        // export fluid
        PlyWriter::save(_engine.fluidMesh(), "fluid.ply");
        os << "scene.objects.fluid.material = fluid" << std::endl;
        os << "scene.objects.fluid.ply = fluid.ply" << std::endl;

        // export boundaries
        const auto &boundaryMeshes = _engine.sph().boundaryMeshes();
        for (size_t i = 0; i < boundaryMeshes.size(); ++i) {
            std::string name = tfm::format("boundary-%03d", i);
            std::string filename = name + ".ply";
            PlyWriter::save(boundaryMeshes[i], filename);
            os << tfm::format("scene.objects.%s.material = boundary", name) << std::endl;
            os << tfm::format("scene.objects.%s.ply = %s", name, filename) << std::endl;
        }

        // export floor
        Box3f box = _engine.scene().world.bounds;
        box.max.y() = box.min.y();
        box.min.y() -= 0.01f * box.extents().norm();
        Mesh floorMesh = Mesh::createBox(box);
        PlyWriter::save(floorMesh, "floor.ply");
        os << "scene.objects.floor.material = floor" << std::endl;
        os << "scene.objects.floor.ply = floor.ply" << std::endl;
    }


    std::cout << std::endl;
    std::cout << tfm::format("Rendering frame %d with SLG4 ...", _frameIndex) << std::endl;

    try {
        std::string slg = "/home/freak/luxmark-v3.1/slg4";
        std::string arguments = "-o render.cfg";
        exec_stream_t es;
        // Set 1 minutes timeout
        es.set_wait_timeout(exec_stream_t::s_out | exec_stream_t::s_err | exec_stream_t::s_child, 60*1000);
        es.start(slg, arguments);

        std::string s;
        while (std::getline(es.err(), s).good()) {
            std::cout << s << std::endl;
            if (StringUtils::endsWith(StringUtils::trim(s), "Done.")) {
                DBG("detected finish");
                break;
            }
        }

        DBG("closing slg");
        if (!es.close()) {
            std::cerr << "SLG4 timeout!" << std::endl;
        }
        DBG("closed");
    } catch (const std::exception &e) {
        std::cerr << "exec-stream error: " << e.what() << "\n";
    }

    FileUtils::changeCurrentDir(currentDir);
}

void Renderer::setupEmptyDirectory(const filesystem::path &path) {
    FileUtils::createDir(path.str());
    tinydir_dir dir;
    tinydir_open_sorted(&dir, path.str().c_str());
    for (size_t i = 0; i < dir.n_files; ++i) {
        tinydir_file file;
        tinydir_readfile_n(&dir, &file, i);
        if (file.is_reg) {
            FileUtils::deleteFile((path / file.name).str());
        }
    }
    tinydir_close(&dir);
}

} // namespace pbs
