#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2

#include <iostream>
#include <vector>
#include <fstream>

#include "nlohmann_json.hpp"
// for convenience
using json = nlohmann::json;
using namespace std;
#if __EMSCRIPTEN__
#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/emscripten.h>

EM_JS(void, MountFilesystem, (const char *repo_root_path), {
    // inside of emscripten
    console.log("MountFilesystem mkdir c++");
    FS.mkdir("/working");

    console.log("MountFilesystem c++");
    FS.mount(NODEFS, {root : Module.UTF8ToString(repo_root_path)}, "/working");
});

#endif //__EMSCRIPTEN__
struct Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()
struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    {                     // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

json read_json_file(std::string file_path)
{
    std::ifstream i(file_path);
    json j;
    i >> j;
    return j;
}

Refl_t extract_reflection_type(std::string r_type){
    Refl_t reflection_type;
    if (r_type == "DIFF")
    {
        reflection_type = DIFF;
    }
    else if (r_type == "SPEC")
    {
        reflection_type = SPEC;
    }
    else
    {
        reflection_type = REFR;
    }
    return reflection_type;
}

std::vector<Sphere> spheres;
void load_scene(std::string scene_data)
{
    std::cout << "load_scene(std::string scene_data)" << std::endl;
    std::cout << scene_data << std::endl;
    json data = json::parse(scene_data); //read_json_file(scene_file);

    std::cout << "dataSpheres" << std::endl;
    for(auto &sphere : data["Spheres"]){
        std::cout << "------------------------------------" << std::endl;

        int radius = sphere["radius"];
        std::cout << "radius: " << radius << std::endl;

        Vec position{sphere["poistion"][0], sphere["poistion"][1], sphere["poistion"][2]};
        std::cout << "position: " << sphere["poistion"] << std::endl;

        Vec emission{sphere["emission"][0], sphere["emission"][1], sphere["emission"][2]};
        std::cout << "emission: " << sphere["emission"] << std::endl;

        Vec color{sphere["color"][0], sphere["color"][1], sphere["color"][2]};
        std::cout << "color: " << sphere["color"] << std::endl;

        Refl_t reflection_type = extract_reflection_type(sphere["reflection_type"]);
        std::cout << "reflection_type: " << sphere["reflection_type"] << std::endl;

        spheres.push_back(Sphere(radius, position, emission, color, reflection_type));
    }
}

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id)
{
    int n = spheres.size();
    double d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id))
        return Vec();                // if miss, return black
    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z; // max refl
    if (++depth > 5){
        if (erand48(Xi) < p)
        {
            f = f * (1 / p);
        }
        else
        {
            return obj.e; //R.R.
        }
    }

    if (obj.refl == DIFF)
    { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    else if (obj.refl == SPEC) // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.e + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                           radiance(reflRay, depth, Xi) * RP
                                                       : radiance(Ray(x, tdir), depth, Xi) * TP)
                                    : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}
int render(int spp)
{
    int w = 1024, h = 768;
    int samps = spp / 4;                                                  // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        auto spp_info = samps * 4;
        auto percent_complete = 100.0 * y / (h - 1);
        std::cout << "spp_info:" << spp_info << " status: " << percent_complete << std::endl;
        //fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)}; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)       // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }

    std::cout << "write render to file" << std::endl;
#if __EMSCRIPTEN__
    std::cout << "/working/image.ppm" << std::endl;
    FILE *f = fopen("/working/image.ppm", "w"); // Write image to PPM file.
#else
    std::cout << "image.ppm" << std::endl;
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
#endif //__EMSCRIPTEN__
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++){
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
    fclose(f);
    std::cout << "DOne" << std::endl;
    return 0;
}

#if __EMSCRIPTEN__

int renderJS(int spp, std::string scene_data)
{
    std::cout << "hello from renderJS" << std::endl;
    std::cout << " argv[1]: " << spp << std::endl;

    std::string cur_dir = "./";
    std::cout << "MountFilesystem: " << cur_dir << std::endl;
    MountFilesystem(cur_dir.c_str());

    std::cout << "load_scene" << std::endl;
    load_scene(scene_data);

    std::cout << "render" << std::endl;
    int ret = render(spp);
    return ret;
}

int main(int argc, char *argv[]){
    std::cout << "hello from wasm main function" << std::endl;
    /*****
    //first argv is spp (sample per pixel)
    //second argv is scene data
    if(argc < 3){
        std::cout << "Usage error: ./smallpt <SPP> <SCENE_DATA>" << std::endl;
        return 1;
    }
    ****/
    // std::string ssp = argv[1];
    // std::cout << " argv[1]: " << ssp << std::endl;

    // std::string scene_data = argv[2];
    // std::cout << " argv[2]: " << scene_data << std::endl;

    // std::string cur_dir = "./";
    // std::cout << "MountFilesystem: " << cur_dir << std::endl;
    // MountFilesystem(cur_dir.c_str());

    // std::cout << "load_scene: " << scene_data << std::endl;
    // load_scene(scene_data);

    // std::cout << "render" << std::endl;
    // return render(spp);
}

EMSCRIPTEN_BINDINGS(renderJS)
{
    emscripten::function("renderJS", &renderJS);
}

#else

json spheresToJson(std::vector<Sphere> spheres)
{
    json j;
    j["Spheres"] = json::array();
    for(auto& sphere : spheres){
        json json_object = json::object();
        json_object["radius"] = sphere.rad;
        json_object["poistion"] = json::array({sphere.p.x, sphere.p.y, sphere.p.z});
        json_object["emission"] = json::array({sphere.e.x, sphere.e.y, sphere.e.z});
        json_object["color"] = json::array({sphere.c.x, sphere.c.y, sphere.c.z});
        if (sphere.refl == DIFF){
            json_object["reflection_type"] = "DIFF";
        }
        else if (sphere.refl == SPEC)
        {
            json_object["reflection_type"] = "SPEC";
        }else{
            json_object["reflection_type"] = "REFR";
        }
        j["Spheres"].push_back(json_object);
    }
    return j;
}

int main(int argc, char *argv[])
{
    std::cout << "hello from wasm" << std::endl;

    /*****
    //first argv is spp (sample per pixel)
    //second argv is scene data
    if(argc < 3){
        std::cout << "Usage error: ./smallpt <spp> <SCENE_DATA>" << std::endl;
        return 1;
    }
    ****/
    std::string spp = argv[1];
    std::cout << " argv[1]: " << spp << std::endl;

    std::string scene_data = argv[2];
    std::cout << " argv[2]: " << scene_data << std::endl;

    std::cout << "load_scene: " << scene_data << std::endl;
    load_scene(scene_data);

    std::cout << "render" << std::endl;
    return render(atoi(spp.c_str()));
}
#endif //__EMSCRIPTEN__