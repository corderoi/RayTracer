//
//  raytrace.cpp
//  RayTracer
//
//  Modified by Ian Cordero on 12/1/15.
//  CS 174A | Prof. Diana Ford | UCLA Fall 2015 | Dis 1B | TA: Donyang
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept> // runtime_error, logic_error
#include <algorithm> // min()
#include <cassert> // assert()
using namespace std;

//////////////////////////////////////////////////////////////////////

// Configuration constants

const int MaximumReflectionDepth = 3;
const string DefaultOutputFileName = "output.ppm";

//////////////////////////////////////////////////////////////////////

// General utilities

/**
 * Converts any object with the << operator to string
 * @param obj T the object to convert
 * @return string representation of obj
 */
template <typename T>
string toString(T obj) {
    ostringstream oss;
    oss << obj;
    return oss.str();
}

/**
 * Exception representing programming error
 */
class LogicException: public logic_error {
public:
    LogicException() : logic_error("LogicException") {}
    LogicException(string message) : logic_error(message) {}
};

/**
 * Exception representing I/O or other runtime errors
 */
class RuntimeException: public runtime_error {
public:
    RuntimeException() : runtime_error("RuntimeException") {}
    RuntimeException(string message) : runtime_error(message) {}
};

/**
 * Checks if a string is only white space
 * @param text The string to check
 * @return true if 'text' is empty or only white space, false otherwise
 */
bool isWhiteSpace(const string& text) {
    for (int i = 0; i < text.length(); i++) {
        const char c = text[i];
        if (!isspace(c)) {
            return false;
        }
    }
    return true;
}

/**
 * Assistance with debugging
 */
class Debug {
public:
    #define Debug_ENABLED__ false

    static void assertC(bool condition) {
        #if Debug_ENABLED__
        assert(condition);
        #endif
    }
    
    static void log(string message) {
        #if Debug_ENABLED__
        cerr << message << endl;
        #endif
    }
    
    #undef Debug_ENABLED__
};

//////////////////////////////////////////////////////////////////////

// Core classes/structs

struct Ray
{
    vec4 origin;
    vec4 dir;
};

/**
 * Represents a sphere object in the scene.
 */
struct Sphere {
    string name;
    vec4 pos;
    vec3 scale;
    vec4 color;
    float Ka,
        Kd,
        Ks,
        Kr;
    float n;
};

/**
 * Represents a point light source in the scene.
 */
struct Light {
    string name;
    vec4 pos;
    vec4 color;
};

//////////////////////////////////////////////////////////////////////

// Scene parameters

int g_width;
int g_height;

vector<vec4> g_colors;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

vec4 bgColor;
vec4 ambientLightIntensity;

string outputFileName = DefaultOutputFileName;

vector<Sphere> spheres;
vector<Light> lights;


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

/**
 * Represents all valid descriptors to be expected in the input file.
 */
struct Descriptor {
    enum ID {
        NEAR, LEFT, RIGHT, BOTTOM, TOP, RES, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT
    };
    
    /**
     * Parses a string to determine its corresponding descriptor.
     * @arg text: string The string to examine
     * @return ID The corresponding descriptor ID
     * @throws RuntimeException if the string is unrecognized
     */
    static ID parse(const string& text) {
        #define dm_(n) if (text == #n) { return n; } else
        dm_(NEAR) dm_(LEFT) dm_(RIGHT) dm_(BOTTOM) dm_(TOP) dm_(RES) dm_(SPHERE) dm_(LIGHT) dm_(BACK) dm_(AMBIENT) dm_(OUTPUT)
        #undef dm_
        throw RuntimeException("Invalid descriptor field received: " + text);
    }
};


void parseLine(const vector<string>& vs)
{
    // Ignore line if empty string or white space
    if (isWhiteSpace(vs[0])) {
        return;
    }
    
    const Descriptor::ID d = Descriptor::parse(vs[0]);
    
    // Handle each possible descriptor case
    // Set the corresponding variables to parsed values. It is assumed that all
    // required values will be specified (besides the optional and multiple potential
    // number of lights/spheres), otherwise behavior is undefined
    switch (d) {
        case Descriptor::NEAR:
            g_near = toFloat(vs[1]);
            break;
        case Descriptor::LEFT:
            g_left = toFloat(vs[1]);
            break;
        case Descriptor::RIGHT:
            g_right = toFloat(vs[1]);
            break;
        case Descriptor::BOTTOM:
            g_bottom = toFloat(vs[1]);
            break;
        case Descriptor::TOP:
            g_top = toFloat(vs[1]);
            break;
        case Descriptor::RES:
            g_width = (int)toFloat(vs[1]);
            g_height = (int)toFloat(vs[2]);
            g_colors.resize(g_width * g_height);
            break;
        case Descriptor::SPHERE:
        {
            // Push a new sphere onto 'spheres'
            const string name = vs[1];
            const vec4 pos = toVec4(vs[2], vs[3], vs[4]);
            const vec3 scale = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
            const vec4 color = toVec4(vs[8], vs[9], vs[10]);
            const float Ka = toFloat(vs[11]),
                Kd = toFloat(vs[12]),
                Ks = toFloat(vs[13]),
                Kr = toFloat(vs[14]);
            const float n = toFloat(vs[15]);
            
            const Sphere s = { name, pos, scale, color, Ka, Kd, Ks, Kr, n };
            
            spheres.push_back(s);
        }
            break;
        case Descriptor::LIGHT:
        {
            // Push a new light onto 'lights'
            const string name = vs[1];
            const vec4 pos = toVec4(vs[2], vs[3], vs[4]);
            const vec4 color = toVec4(vs[5], vs[6], vs[7]);
            
            const Light l = { name, pos, color };
            
            lights.push_back(l);
        }
            break;
        case Descriptor::BACK:
            bgColor = toVec4(vs[1], vs[2], vs[3]);
            break;
        case Descriptor::AMBIENT:
            ambientLightIntensity = toVec4(vs[1], vs[2], vs[3]);
            break;
        case Descriptor::OUTPUT:
            outputFileName = vs[1];
            break;
        default:
            throw LogicException("Unknown descriptor encountered: " + toString(d));
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

/**
 * Calculates the discriminant of a quadratic equation given its coefficients
 */
static float calculateDiscriminant(float A, float B, float C) {
    return pow(B, 2.0) - A * C;
}

/**
 * Solves the quadratic equation.
 * @return vector containing the root(s), if any
 */
static vector<float> quadratic(float A, float B, float C) {
    vector<float> roots;
    const float discriminant = calculateDiscriminant(A, B, C);
    if (discriminant < 0.0) {
        return roots;
    }
    const float D = B / A;
    if (discriminant == 0.0) {
        roots.push_back(D);
        return roots;
    }
    const float E = sqrt(discriminant) / A;
    roots.push_back(D + E);
    roots.push_back(D - E);
    return roots;
}

/**
 * Purges roots not corresponding to valid intersections for ray tracing, such as roots of 0 or less.
 * @param roots The roots to check
 * @return vector containing roots that correspond to valid intersections
 */
static vector<float> purgeRoots(const vector<float>& roots) {
    vector<float> newRoots;
    const float IntersectionThreshold = 0.0001; // don't let objects cast a shadow on themselves
    for (int i = 0; i < roots.size(); i++) {
        const float root = roots[i];
        if (root > IntersectionThreshold) {
            newRoots.push_back(root);
        }
    }
    return newRoots;
}

/**
 * Checks for intersection of 'sphere' and 'ray'.
 * @param sphere The sphere of interest with regards to intersection
 * @param ray The ray we are checking intersection for
 * @param intersections Rays to potentially set the point and normal of intersection 
 * to, if and only if there is an intersection
 * @return true if ray intersects sphere and sets intersections to the point(s) and
 * normal(s) of intersection, false otherwise.
 */
bool intersect(const Sphere& sphere, const Ray& ray, vector<Ray>& intersections) {
    const vec4 offsetPos = sphere.pos - ray.origin;
    const mat4 scale = Scale(sphere.scale);
    mat4 invertedScale;
    InvertMatrix(scale, invertedScale);
    
    Ray scaledRay = { invertedScale * offsetPos, invertedScale * ray.dir };
    
    const float A = dot(scaledRay.dir, scaledRay.dir),
        B = dot(scaledRay.origin, scaledRay.dir),
        C = dot(scaledRay.origin, scaledRay.origin) - 1.0;
    
    vector<float> roots = purgeRoots(quadratic(A, B, C));
    
    if (roots.size() == 0) { // no real roots of interest
        return false;
    }
    sort(roots.begin(), roots.end()); // guarantee to caller that intersections will be ordered by distance
    
    for (int i = 0; i < roots.size(); i++) {
        const float root = roots[i];
        
        // Calculate point of intersection
        const vec4 pointOfIntersection = ray.origin + root * ray.dir;
        
        // Calculate normal at point of intersection
        vec4 scaledNormal = pointOfIntersection - sphere.pos;
        InvertMatrix(Scale(sphere.scale), invertedScale);
        const vec4 normalAtIntersection = normalize(2.0 * invertedScale * invertedScale * scaledNormal);

        const Ray intersection = { pointOfIntersection, normalAtIntersection };
        
        intersections.push_back(intersection);
    }
    
    return true;
}

bool intersect(const Sphere& sphere, const Ray& ray) { // if we only care about whether an intersection happens or not (i.e. shadow ray)
    vector<Ray> dontCare;
    return intersect(sphere, ray, dontCare);
}

// -------------------------------------------------------------------
// Ray tracing

/**
 * Calculates the distance between two points.
 * @param p1 The first point
 * @param p2 The second point
 * @return float representing the Euclidean distance between the points p1 and p2.
 */
static float calculateDistance(const vec4& p1, const vec4& p2) {
    const float distance = sqrt(pow(p1.x - p2.x, 2.0) + pow(p1.y - p2.y, 2.0) + pow(p1.z - p2.z, 2.0));
    return distance;
}

/**
 * Calculates the reflection ray for a given incident ray and reflection point/normal.
 * @param incident Ray representing the ray to reflect
 * @param normal Ray representing the surface to reflect against
 * @return Ray representing the reflection ray from the surface
 */
static Ray reflect(const Ray& incident, const Ray& normal) {
    const vec4 Ri = incident.dir,
        Rn = normal.dir;
    const vec4 Rr = Ri - 2.0 * Rn * dot(Ri, Rn);
    
    const Ray reflection = { normal.origin, normalize(Rr) };
    
    return reflection;
}

/**
 * Determines if the shadow ray is blocked by any spheres with regards to the light.
 * @param shadowRay The shadow ray
 * @param light The light of interest
 * @param spheres The spheres of interest
 * @return true if any spheres block the shadowRay, false otherwise
 */
static bool shadowRayBlocked(const Ray& shadowRay, const Light& light, const vector<Sphere>& spheres) {
    for (int i = 0; i < spheres.size(); i++) {
        const Sphere* sphere = &spheres[i];
        if (intersect(*sphere, shadowRay)) {
            return true;
        }
    }
    return false;
}

/**
 * Traces a ray through the scene and computes the color corresponding to the ray.
 * @param ray The ray to trace
 * @param reflectionDepth How many times this ray has been reflected. 0 represents
 * initial camera rays (rays starting from the camera)
 * @return vec4 representing the color that this ray maps to
 */
static vec4 trace(const Ray& ray, int reflectionDepth)
{
    // Only recurse 'MaximumReflectionDepth' levels
    if (reflectionDepth > MaximumReflectionDepth) {
        return bgColor;
    }
    
    // Find closest intersection P of ray r with objects
    Ray closestIntersection;
    float minimumDistance = INFINITY;
    Sphere* intersectionSphere = nullptr;
    
    for (int i = 0; i < spheres.size(); i++) {
        const Sphere* sphere = &spheres[i];
        
        vector<Ray> intersections;
        if (intersect(*sphere, ray, intersections)) {
            for (int j = 0; j < intersections.size(); j++) {
                const Ray* intersection = &intersections[j];
                
                // Calculate distance of intersection point
                const float distance = calculateDistance(ray.origin, intersection->origin);
                
                if (distance < minimumDistance) {
                    // If this is an initial camera ray, cut off any intersection points before the image plane
                    const bool cutOffByImagePlane = (reflectionDepth == 0 && fabs(intersection->origin.z) < g_near);
                    
                    if (!cutOffByImagePlane) {
                        closestIntersection = *intersection;
                        minimumDistance = distance;
                        intersectionSphere = const_cast<Sphere*>(sphere);
                    }
                }
            }
        }
    }
    
    if (intersectionSphere == nullptr) { // no intersection at all
        // Return the background color
        return bgColor;
    }
    
    // Calculate ambient illumination
    vec4 ambientIllumination = intersectionSphere->Ka * (ambientLightIntensity * intersectionSphere->color);
    
    // Calculate illumination of point P with respect to the scene's lights
    // clocal = sum(shadowRays(P, Light))
    vec4 lightIllumination;
    
    const bool insideSphere = (dot(ray.dir, closestIntersection.dir) > 0.0);

    const vec4 V = normalize(ray.origin - closestIntersection.origin);
    
    for (int i = 0; i < lights.size(); i++) {
        const Light* light = &lights[i];
        const Ray shadowRay = { closestIntersection.origin, normalize(light->pos - closestIntersection.origin) };
        const vec4 R = 2.0 * dot(closestIntersection.dir, shadowRay.dir) * closestIntersection.dir - shadowRay.dir;
        
        if (!shadowRayBlocked(shadowRay, *light, spheres)) {
            const float diffuseContribution = dot(closestIntersection.dir, shadowRay.dir);
            
            if (diffuseContribution < 0.0) {
                continue;
            }
            const vec4 diffuseIllumination = intersectionSphere->Kd * light->color * diffuseContribution * intersectionSphere->color;
            lightIllumination += diffuseIllumination;
            
            if (!insideSphere) { // only consider specular illumination for outer surface
                const vec4 specularIllumination =  intersectionSphere->Ks * light->color * pow(dot(R, V), intersectionSphere->n);
                lightIllumination += specularIllumination;
            }
        }
    }
    
    // Calculate contribution of reflections to point P
    // c_rfl = raytrace(r_rfl)
    vec4 reflectionIllumination;
    const Ray reflection = reflect(ray, closestIntersection);
    
    const vec4 reflectionContribution = trace(reflection, reflectionDepth + 1);
    
    if (!(reflectionContribution.x == bgColor.x && reflectionContribution.y == bgColor.y && reflectionContribution.z == bgColor.z)) {
            // omit contribution of reflection if it only produces the background color
        reflectionIllumination = intersectionSphere->Kr * reflectionContribution;
    }
    
    // Reduce calculated data to single color using the following illumination model:
    // return c = clocal + k_rfl * c_rfl + k_rfa * c_rfa
    const vec4 pixelColor = ambientIllumination + lightIllumination + reflectionIllumination;
    
    return pixelColor;
}

vec4 trace(const Ray& ray) {
    return trace(ray, 0); // initial ray (camera ray)
}

vec4 getDir(int ix, int iy)
{
    // Return the direction from the origin to pixel (ix, iy), normalized.
    const float X = g_left + (g_right - g_left) * (static_cast<float>(ix) / (static_cast<float>(g_width) - 1.0)),
        Y = g_bottom + (g_top - g_bottom) * (static_cast<float>(iy) / (static_cast<float>(g_height) - 1.0)),
        Z = -fabs(g_near);
    
    return normalize(vec4(X, Y, Z, 0.0f));
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++) {
        for (int x = 0; x < g_width; x++) {
            for (int i = 0; i < 3; i++) {
                const int BufIndex = y * g_width * 3 + x * 3 + i;
                const int ColorsIndex = y * g_width + x;
                const float PixelColor = min(1.0f, max(0.0f, ((float*) g_colors[ColorsIndex])[i])); // bind to interval [0.0, 1.0]
                buf[BufIndex] = (unsigned char) (PixelColor * 255.9f);
            }
        }
    }
    savePPM(g_width, g_height, (char*) outputFileName.c_str(), buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: " << "raytrace" << " <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

