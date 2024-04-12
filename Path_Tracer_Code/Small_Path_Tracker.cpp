#include<stdio.h>
#include<math.h>
#include<stdlib.h>

// It's convention to name classes and structures in UpperCamelCase
struct Vec{
    double x, y, z;

    // Constructor with default values

    Vec (double x_ = 0, double y_ = 0, double z_ = 0){
        x = x_;
        y = y_;
        z = z_;
    }

    // Defines the way two Vec objects are added. operator+ is a way to overload the standard C++ addition operator for stuff in a structure.
    // const used at the end of the member function basically says that the function doesn't change the object it operates on.
    Vec operator+ (const Vec &b) const{
        return Vec(x + b.x, y + b.y, z + b.z);
    }

    Vec operator- (const Vec&b) const{
        return Vec(x - b.x, y - b.y, z - b.z);
    }

    // For scalar multiplication
    Vec operator* (double b) const{ 
        return Vec(x * b, y * b, z * b);
    }

    // Component wise multiplication
    Vec mult (const Vec &b) const{
        return Vec(x * b.x, y * b.y, z * b.z);
    }

    // Unit vector or normalized vector
    Vec& norm(){
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    // Dot product
    double dot (const Vec &b) const{
        return x * b.x + y * b.y + z * b.z;
    }

    // Cross product
    // % operator is being overloaded to represent the cross product between the Vec object
    Vec operator% (Vec &b){
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

// It's convention to name classes and structures in UpperCamelCase
struct Ray{
    Vec o, d; //origin and direction Vectors

    // Constructor with an initializer list. A more efficient way of initializing.
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

// Diffusive reflection, specular reflection and refractive reflection
enum Refl_t {DIFF, SPEC, REFR};

// It's convention to name classes and structures in UpperCamelCase
struct Sphere{
    double rad; //radius
    Vec p, e, c; // position, emission, color. Emission vector specifies the direction and intensity of emitted light. Color is an RGB vector.
    Refl_t refl; // reflection type

    //Constructor with and initializer list.
    Sphere (double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    // Function to determine the intersection of a ray with the sphere specified

    double intersect (const Ray &r) const{
        Vec op = p - r.o; // Vector joining ray origin to sphere centre

        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;

        if (det < 0)
            return 0; // returns 0 if no interesction
        else
            det = sqrt(det);

        // returns distance along the ray to the intersection point. ? Is a ternanry operators. Look up the syntax to understand.
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[] = {// Array of spheres to set the scene
    Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left 
    Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Right 
    Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back 
    Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Front 
    Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Bottomm 
    Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top 
    Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirror 
    Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glass 
    Sphere(600, Vec(50,681.6 -0.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Light
};

inline double clamp(double x) {return x < 0 ? 0 : x > 1 ? 1 : x;} // Function essentially bounds values into [0,1]

inline int toInt(double x) {return int(pow(clamp(x), 1/2.2) * 255 + 0.5);} // 2.2 is gamma correction factor, 255 is used to scale the color value in [0,256]

inline bool intersect(const Ray &r, double &t, int &id){
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20; // number of spheres, distance and maximum distance.
    
    for (int i = int(n); i--;){ //iterates through every sphere in the scene
        if ((d = spheres[i].intersect(r)) && d < t){ // Checks if the distance along ray for ith sphere is the minimum with amongst all the spheres so far.
            t = d; // Reassigns the minimum distance along ray
            id = i; // Notes down the id of the closest sphere
        }
    }
    return t < inf; // Returns bool
}

// This function returns the color of a pixel corresponding to a ray.
Vec radiance(const Ray &r, int depth, unsigned short *Xi){ // *Xi is the pointer to an array of possible seeds.
    double t; // Distance to intersection
    int id = 0; // ID of intersected object

    // Returning black color if no intersection

    if (!intersect(r, t, id))
        return Vec();

    const Sphere &obj = spheres[id]; // The object hit by the ray
    Vec x = r.o + r.d * t; // The intersection point
    Vec n = (x - obj.p).norm(); // Normal vector at the intersection
    Vec n1 = n.dot(r.d) < 0 ? n : n * -1; // Ensuring the normal vector makes outwards from the object
    Vec f = obj.c; // Surface color of the object
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // Max reflection coefficient

    // Russian roulette : randomly terminating rays to avoid infinite recursion
    if (++depth > 5) // If depth of recursion (number of reflections or refractions) exceeds a number, terminate the ray early
        if (erand48(Xi) < p)
            f = f * (1 / p); // If the randomg number is < p, the ray's intensity is scaled 1/p to conserve energy by compensating for early termination.
        else
            return obj.e;

    // Ideal diffuse reflection
    if (obj.refl == DIFF){
        double r1 = 2 * M_PI * erand48(Xi); // Generates a random angle in [0, 2pi]
        double r2 = erand48(Xi); // cosine-weighted distribution of driections used to simulate more likely reflections along direction closer to normal.
        double r2s = sqrt(r2); // Square root of r2

        Vec w = n1; // Corrected surface normal
        Vec u = (fabs(w.x) > 0.1 ? Vec(0, 1) : Vec (1))%w.norm(); // u is a 2D tangent vector to n1. The condition checks if w is parallel to the y-axis.
        Vec v = w%u; // Bitangent vector perpendicular to both u and w
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); // Random reflected direction following a cosine-weighted distribution
        // Note that the direction of the reflected ray is basically a linear combination of u, v and w.

        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }

    // Ideal specular reflection
    else if (obj.refl == SPEC){
        return obj.e + f.mult(radiance(Ray(x, r.d - n1 * 2 * n1.dot(r.d)), depth, Xi));
    }

    // Ideal Dielectric Refraction
    Ray reflRay(x, r.d - n1 * 2 * n1.dot(r.d));
    bool into = n.dot(n1)> 0; // Checking if the normal and corrected normal are pointing in the same way
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(n1), cos2t;
    // nc is the mu of medium 1. nt is the mu of medium 2. nnt is the ratio of refractive indices depending on whether the light is going in or out.
    // ddn is the cosine of the angle between the ray and the correct normal
    
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Checks is the angle of refraction is obtuse. If so TIR occurs
        return obj.e + f.mult(radiance(reflRay, depth, Xi)); // In case of TIR, only the reflected component of the ray is considered

    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt  + sqrt(cos2t)))).norm(); // Direction of refracted ray

    double a = nt - nc, b =  nt + nc, RO = a * a / (b * b), c = 1 - (into ? - ddn : tdir.dot(n));
    double Re = RO + (1 - RO) * c * c * c * c * c, Tr = 1 - Re, P = 0.25 + 0.5 * Re, RP = Re / P, TP = Tr / (1 - P); 
    // Fresnel equation to determine coefficients of reflection and transmission
    // RP and TP are reflection and transmission probabilities
    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) + Tr); 
    // For more than two recursions, there is a Russian roulettet terminate either the reflected or transmitted ray.
    // For less than two recursions, the new ray is simply a linear combination of the reflected and transmitted rays.
}

int main(int argc, char *argv[]){
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // Number of samples per pixel and dimensions of the image.
    // atoi("12345") returns 12345. Converts string numbers to numbers. 
    // samps has a /4 as each pixel is divided into a 2x2 grid of subpixels to enable supersampling, an anti-aliasing technique.

    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // Camera position and direction
    Vec cx = Vec(w * 0.5135 / h), cy = (cx % cam.d).norm() * 0.5315, r, *c = new Vec[w * h];
    // cx is the x-axis of the camera. cy is the y-axis of the camera which is the cross product of cx and direction of the camera.
    // r is used to obtain radiance value of a particular pixel. c is a pointer to an array of vectors. There is 1 color vector for each pixel.
    // "new" is used the mention that the memory to c is allocated dynamically.

#pragma omp parallel for schedule(dynamic, 1) private(r)
// a compiler directive to use OpenMP that will parallelise the following loops by dynamically allocating one interation to a thread at any time.
// private(r) means that each thread has it's own private copy of r to avoid data conflicts when multiple threads are operating on r.

    for (unsigned short y = 0; y < h; y++){
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.0 * y / (h - 1));
        // %5.2f%% says that each number has at least 5 characters (with spaces padded to the left) and 2 digits after the decimal point.

        for (unsigned short x = 0, Xi[3] = {0, 0, y}; x < w; x++)
            for (int sy = 0, i = (h - y -1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()){ // r is reset to 0 at the beginning of each subpixel to obtain radiance of each subpizel.
                    for (int s = 0; s < samps; s++){
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        // Introduce random offsets dx and dy within the subpixel. Called Jittered sampling, anti-aliasing

                        Vec d = cx * (((sx + 0.5 + dx) / 2 + x) / w - 0.5) + cy * (((sy + 0.5 + dy) / 2 + y) / h - 0.5) + cam.d;
                        // Adjusts the direction of the camera by considering the camera to look exactly at the direction at the point from within the subpixel
                        // from which the sampling is done from. / w and / h are done to normalize to [0,1]. - 0.5 is then performed to make sure centre of 
                        // each subpixel is (0,0).

                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1.0 / samps);
                        // Radiance accumulation. Dividing by samps to normalize the color obtained.
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
    }

    FILE *f = fopen("image.ppm", "w"); // Writing image to ppm file
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}