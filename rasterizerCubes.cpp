#include <iostream>
#include <fstream>
#include <cmath>

#define WIDTH 500
#define HEIGHT 500
#define INF 1e9

using namespace std;

class Vec3 {
public:
    float x, y, z;
    Vec3() {}
    Vec3(float a, float b, float c): x(a), y(b), z(c) {}

    Vec3 operator+(Vec3 v) { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(Vec3 v) { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator*(float s) { return Vec3(x*s, y*s, z*s); }
};

class Plane {
public:
    Vec3 normal;
    float D;
};

class Triangle {
public:
    Vec3 v0, v1, v2;
    unsigned char color[3];
};

class BoundingSphere {
public:
    Vec3 center;
    float radius;
};

class Instance {
public:
    Triangle triangles[12];
    BoundingSphere sphere;
};

class Scene {
public:
    Instance instances[3];
    int count = 0;
};

float Dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

float SignedDistance(Plane p, Vec3 pt) {
    return Dot(p.normal, pt) + p.D;
}

Vec3 Intersect(Vec3 A, Vec3 B, Plane p) {
    Vec3 AB = B - A;
    float t = -(Dot(p.normal, A) + p.D) / Dot(p.normal, AB);
    return A + AB * t;
}

bool Inside(Plane p, Vec3 v) {
    return SignedDistance(p, v) >= 0;
}

int ClipTriangle(Triangle in, Plane p, Triangle out[2]) {
    Vec3 v[3] = {in.v0, in.v1, in.v2};
    float d[3];
    bool inside[3];
    int n_in = 0;

    for (int i = 0; i < 3; i++) {
        d[i] = SignedDistance(p, v[i]);
        inside[i] = d[i] >= 0;
        if (inside[i]) n_in++;
    }

    if (n_in == 3) {
        out[0] = in;
        return 1;
    }
    if (n_in == 0) return 0;

    Vec3 in_v[2], out_v[2];
    int i_in = 0, i_out = 0;

    for (int i = 0; i < 3; i++) {
        if (inside[i]) in_v[i_in++] = v[i];
        else out_v[i_out++] = v[i];
    }

    if (n_in == 1) {
        Vec3 A = in_v[0];
        Vec3 Bp = Intersect(A, out_v[0], p);
        Vec3 Cp = Intersect(A, out_v[1], p);
        out[0] = {A, Bp, Cp, {in.color[0], in.color[1], in.color[2]}};
        return 1;
    }

    Vec3 A = in_v[0], B = in_v[1];
    Vec3 Ap = Intersect(A, out_v[0], p);
    Vec3 Bp = Intersect(B, out_v[0], p);
    out[0] = {A, B, Ap, {in.color[0], in.color[1], in.color[2]}};
    out[1] = {Ap, B, Bp, {in.color[0], in.color[1], in.color[2]}};
    return 2;
}

int ClipInstance(Instance in, Plane p, Triangle out_tris[24]) {
    int count = 0;
    for (int i = 0; i < 12; i++) {
        Triangle clipped[2];
        int n = ClipTriangle(in.triangles[i], p, clipped);
        for (int j = 0; j < n; j++) {
            out_tris[count++] = clipped[j];
        }
    }
    return count;
}

unsigned char framebuffer[HEIGHT][WIDTH][3];
float zbuffer[HEIGHT][WIDTH];

int Project(Vec3 v, int& x, int& y) {
    if (v.z == 0) return 0;
    float fx = v.x / v.z;
    float fy = v.y / v.z;
    x = (int)((fx + 1) * WIDTH / 2);
    y = (int)((1 - fy) * HEIGHT / 2);
    return 1;
}

void DrawTriangle(Triangle t) {
    int x0, y0, x1, y1, x2, y2;
    if (!Project(t.v0, x0, y0)) return;
    if (!Project(t.v1, x1, y1)) return;
    if (!Project(t.v2, x2, y2)) return;

    int minX = min(x0, min(x1, x2));
    int maxX = max(x0, max(x1, x2));
    int minY = min(y0, min(y1, y2));
    int maxY = max(y0, max(y1, y2));

    float denom = (float)((y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2));
    if (denom == 0) return;

    for (int y = minY; y <= maxY; y++) {
        if (y < 0 || y >= HEIGHT) continue;
        for (int x = minX; x <= maxX; x++) {
            if (x < 0 || x >= WIDTH) continue;
            float l1 = ((y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)) / denom;
            float l2 = ((y2 - y0)*(x - x2) + (x0 - x2)*(y - y2)) / denom;
            float l3 = 1 - l1 - l2;
            if (l1 >= 0 && l2 >= 0 && l3 >= 0) {
                float z = l1*t.v0.z + l2*t.v1.z + l3*t.v2.z;
                if (z < zbuffer[y][x]) {
                    zbuffer[y][x] = z;
                    framebuffer[y][x][0] = t.color[0];
                    framebuffer[y][x][1] = t.color[1];
                    framebuffer[y][x][2] = t.color[2];
                }
            }
        }
    }
}

Triangle MakeTri(Vec3 a, Vec3 b, Vec3 c, unsigned char col[3]) {
    Triangle t;
    t.v0 = a; t.v1 = b; t.v2 = c;
    t.color[0] = col[0];
    t.color[1] = col[1];
    t.color[2] = col[2];
    return t;
}

Instance MakeCube(Vec3 pos, float scale) {
    Vec3 v[8] = {
        {1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1},
        {1,1,-1},{-1,1,-1},{-1,-1,-1},{1,-1,-1}
    };
    for (int i = 0; i < 8; i++) v[i] = v[i]*scale + pos;

    unsigned char red[3]     = {255, 0, 0};
    unsigned char green[3]   = {0, 255, 0};
    unsigned char blue[3]    = {0, 0, 255};
    unsigned char yellow[3]  = {255, 255, 0};
    unsigned char magenta[3] = {255, 0, 255};
    unsigned char cyan[3]    = {0, 255, 255};

    Triangle tris[12] = {
        MakeTri(v[0],v[1],v[2], red),     MakeTri(v[0],v[2],v[3], red),      // front
        MakeTri(v[4],v[0],v[3], green),   MakeTri(v[4],v[3],v[7], green),    // left
        MakeTri(v[5],v[4],v[7], blue),    MakeTri(v[5],v[7],v[6], blue),     // back
        MakeTri(v[1],v[5],v[6], yellow),  MakeTri(v[1],v[6],v[2], yellow),   // right
        MakeTri(v[4],v[5],v[1], magenta), MakeTri(v[4],v[1],v[0], magenta),  // bottom
        MakeTri(v[2],v[6],v[7], cyan),    MakeTri(v[2],v[7],v[3], cyan)      // top 
    };

    Instance inst;
    for (int i = 0; i < 12; i++) inst.triangles[i] = tris[i];
    inst.sphere.center = pos;
    inst.sphere.radius = scale * 1.73;
    return inst;
}

void ClearBuffers() {
    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++) {
            framebuffer[y][x][0] = framebuffer[y][x][1] = framebuffer[y][x][2] = 0;
            zbuffer[y][x] = INF;
        }
}

void SavePPM(const char* fname) {
    ofstream out(fname);
    out << "P3\n" << WIDTH << " " << HEIGHT << "\n255\n";
    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++)
            out << (int)framebuffer[y][x][0] << " " << (int)framebuffer[y][x][1] << " " << (int)framebuffer[y][x][2] << "\n";
    out.close();
}

int main() {
    Scene scene;
    scene.instances[0] = MakeCube({0, 0, 6}, 1);
    scene.instances[1] = MakeCube({-4, 2, 9}, 1);
    scene.instances[2] = MakeCube({4, -2, 12}, 0.75);
    scene.count = 3;

    Plane clip = {{0,0,1}, -1};
    ClearBuffers();

    for (int i = 0; i < scene.count; i++) {
        Triangle clipped[24];
        int n = ClipInstance(scene.instances[i], clip, clipped);
        for (int j = 0; j < n; j++) {
            DrawTriangle(clipped[j]);
        }
    }

    SavePPM("rasterizerCubes.ppm");
    cout << "Saved to rasterizerCubes.ppm";
    return 0;
}
