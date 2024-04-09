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

    // Norm
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
    Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Light
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