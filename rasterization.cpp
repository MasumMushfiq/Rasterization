#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <random>
#include "bitmap_image.hpp"

using namespace std;

using pdd = pair<double, double>;
using ppdd = pair<pdd, pdd>;

const double PI = 3.141592653589793;
const double EPSILON = 1.0e-6;

class Vector {
public:
    double x, y, z;

    Vector() { x = y = z = 0; }

    // constructs a vector with given components
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize() {
        double r = sqrt(x * x + y * y + z * z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vector operator+(const Vector &v) const {
        Vector v1(x + v.x, y + v.y, z + v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector &v) const {
        Vector v1(x - v.x, y - v.y, z - v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vector operator*(double m) const {
        Vector v(x * m, y * m, z * m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b) {
        Vector v(a.y * b.z - a.z * b.y, b.x * a.z - b.z * a.x, a.x * b.y - a.y * b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print() const {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};

typedef class homogeneous_point {
public:
    double x, y, z, w;

    // set the three coordinates, set w to 1
    homogeneous_point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    /*
    default constructor. does nothing. allows declarations like below:
        matrix m;
    therefore, usage is dangerous
    */
    homogeneous_point() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
        this->w = 1;
    }

    // constructs a homogeneous point with given coordinates. forces w to be 1.0
    // if w is zero, raises error
    homogeneous_point(double x, double y, double z, double w) {
        assert(w != 0);
        this->x = x / w;
        this->y = y / w;
        this->z = z / w;
        this->w = 1;
    }

    // adds two points. returns a point forcing w to be 1.0
    homogeneous_point operator+(const homogeneous_point &point) const {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    homogeneous_point operator+(const Vector &v) const {
        return homogeneous_point(x + v.x, y + v.y, z + v.z);
    }

    // subtracts one point from another. returns a vector
    Vector operator-(const homogeneous_point &point) const {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        Vector v(x, y, z);
        return v;
    }

    bool operator==(const homogeneous_point &p) {
        return abs(x - p.x) < EPSILON and abs(y - p.y) < EPSILON and (z - p.z) < EPSILON;
    }

    // Print the coordinates of a point. exists for testing purpose.
    void print() const {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }
} point;

/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix {
public:
    double values[4][4]{};
    int num_rows, num_cols;

    // only set the number of rows and cols
    matrix(int rows, int cols) {
        assert(rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }

    // prepare an nxn square matrix
    explicit matrix(int n) {
        assert(n <= 4);
        num_rows = num_cols = n;
    }

    // prepare and return an identity matrix of size nxn
    static matrix make_identity(int n) {
        assert(n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    // print the matrix. exists for testing purposes
    void print() {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // add the two matrices. Raise error if dimension mismatches
    matrix operator+(const matrix &m) {
        assert(this->num_rows == m.num_rows);
        assert(this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    // subtract a matrix from another. raise error if dimension mismatches
    matrix operator-(const matrix &m) {
        assert(this->num_rows == m.num_rows);
        assert(this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    // multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
    matrix operator*(const matrix &m) {
        assert(this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    // multiply a matrix with a constant
    matrix operator*(double m) {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    // multiply a 4x4 matrix with a homogeneous point and return the resulting point.
    // usage: homogeneous_point p = m * p1;
    // here, m is a 4x4 matrix, intended to be the transformation matrix
    // p1 is the point on which the transformation is being made
    // p is the resulting homogeneous point
    homogeneous_point operator*(const homogeneous_point &p) {
        assert(this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this) * m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    double *operator[](int index) {
        return values[index];
    }

    // return the transpose of a matrix
    matrix transpose() {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }
};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
    double r, g, b;

    color(double r, double g, double b) : r(r), g(g), b(b) {
    }

    color()
    = default;

    void print() const {
        cout << "Color: " << endl;
        cout << r << " " << g << " " << b << endl;
    }
};

// here plane will be parallel to xy plane
struct plane {
    point p;  // point on the plane
    Vector n; // normal to plane

    plane(point p, Vector n) : p(p), n(n) {}

    void print() const {
        cout << "Plane: \n";
        p.print();
        n.print();
    }
};

typedef struct line_segment {
    point p;
    point q;
    Vector v;

    line_segment() = default;

    line_segment(point p, point q) {
        this->p = p;
        this->q = q;
        v = q - p;
    }

    void print() const {
        cout << "Line: \n";
        p.print();
        q.print();
        v.print();
    }
} line;

struct Triangle {
    point p, q, r;

    line a, b, c;

    color colour{};

    Triangle() = default;

    Triangle(point p, point q, point r) {
        this->p = p;
        this->q = q;
        this->r = r;
        a = line(p, q);
        b = line(q, r);
        c = line(r, p);
    }

    void set_sides() {
        a = line(p, q);
        b = line(q, r);
        c = line(r, p);
    }

    void set_color(color colour) {
        this->colour = colour;
    }

    void sort_points_wrt_y() {
        vector<point> ps = {p, q, r};
        sort(ps.begin(), ps.end(), [](point p1, point p2) -> bool { return p1.y > p2.y; }); //descending order
        p = ps[0], q = ps[1], r = ps[2];

    }

    void print() const {
        cout << p.x << " " << p.y << " " << p.z << endl;
        cout << q.x << " " << q.y << " " << q.z << endl;
        cout << r.x << " " << r.y << " " << r.z << endl;
        cout << endl;
    }

    void print_detail() {
        cout << "Triangle: \n";
        p.print();
        q.print();
        r.print();
        a.print();
        b.print();
        c.print();
        cout << endl;
    }
};

double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspect_ratio, near, far;
color background;
int screen_x, screen_y;

vector<color> colors;
vector<Triangle> triangles;

enum CASE {
    INTERSECTING = 0,
    NON_INTERSECTING = 1,
    TOUCHING = 2,
    LYING = 3
};

pair<CASE, point> get_intersecting_point(const line &l, const plane &p) {
    // l.print();
    // p.print();

    double z = p.p.z;
    if ((l.p.z > z and l.q.z > z) or (l.p.z < z and l.q.z < z)) {
//        cout << "Points on the same side" << endl;
        return {NON_INTERSECTING, point()};
    }

    if (Vector::dot(p.n, l.v) == 0) {
//        cout << "Line parallel to plane" << endl;
        return {LYING, l.p};
    }

    double t = Vector::dot((p.p - l.p), p.n) / Vector::dot(l.v, p.n);
    point pnt = l.p + l.v * t;
    if (pnt == l.p or pnt == l.q) {
//        cout << "One point on the plane" << endl;
        return {TOUCHING, pnt};
    }
//    cout << "Intersecting" << endl;
    return {INTERSECTING, pnt};
}

void clip_along_z_near(Triangle tg) {
    // clip by near plane
    plane near_plane(point(0, 0, -near), Vector(0, 0, 1));

    auto inside = [](point p) -> bool { return p.z < -near; };
    auto outside = [](point p) -> bool { return p.z > -near; };
    auto on_the_plane = [](point p) -> bool { return p.z == -near; };

    auto res1 = get_intersecting_point(tg.a, near_plane);
    auto res2 = get_intersecting_point(tg.b, near_plane);
    auto res3 = get_intersecting_point(tg.c, near_plane);

    int cut_values[] = {0, 0, 0, 0};
    cut_values[res1.first]++;
    cut_values[res2.first]++;
    cut_values[res3.first]++;

    if (cut_values[LYING] == 3) { // tg totally on the plane
        triangles.push_back(tg);

    } else if (cut_values[INTERSECTING] == 1 and cut_values[TOUCHING] == 2) { // tg 1 point out, 1 in, 1 on the plane
        point new_point;
        if (inside(tg.p)) {
            if (outside(tg.q)) {
                new_point = res1.second;
                triangles.emplace_back(tg.p, tg.r, new_point);
            } else {
                new_point = res3.second;
                triangles.emplace_back(tg.p, tg.q, new_point);
            }
        } else if (inside(tg.q)) {
            if (outside(tg.p)) {
                new_point = res1.second;
                triangles.emplace_back(tg.q, tg.r, new_point);
            } else {
                new_point = res2.second;
                triangles.emplace_back(tg.q, tg.p, new_point);
            }
        } else if (inside(tg.r)) {
            if (outside(tg.q)) {
                new_point = res2.second;
                triangles.emplace_back(tg.r, tg.p, new_point);
            } else {
                new_point = res3.second;
                triangles.emplace_back(tg.r, tg.q, new_point);
            }
        } else {
            cout << "Should not be here" << endl;
        }

        triangles.back().set_color(tg.colour);
    } else if (cut_values[NON_INTERSECTING] == 3) { // tg totally on 1 side
        if (inside(tg.p)) { // tg totally inside
            triangles.push_back(tg);
        }
    } else if (cut_values[NON_INTERSECTING] == 1 and cut_values[TOUCHING] == 2) { // 1 point on the plane 2 out/in plane
        if (on_the_plane(tg.p)) { // p on the plane
            if (inside(tg.q)) { // tg inside
                triangles.push_back(tg);
            }
        } else if (on_the_plane(tg.q)) { // q on the plane
            if (inside(tg.r)) { // tg inside
                triangles.push_back(tg);
            }
        }
        if (on_the_plane(tg.r)) { // r on the plane
            if (inside(tg.p)) { // tg inside
                triangles.push_back(tg);
            }
        }
    } else if (cut_values[TOUCHING] == 2 and cut_values[LYING] == 1) { // 2 on the plane 1 out/in
        if (inside(tg.p) or inside(tg.q) or inside(tg.r)) {
            triangles.push_back(tg);
        }
    } else if (cut_values[INTERSECTING] == 2 and cut_values[NON_INTERSECTING] == 1) { // 2 one side 1 different
        if (res1.first == NON_INTERSECTING) {
            if (inside(tg.r)) {
                Triangle new_tr = Triangle(tg.r, res2.second, res3.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.r)) {
                Triangle new_tr = Triangle(tg.p, tg.q, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.p, res2.second, res3.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        } else if (res2.first == NON_INTERSECTING) {
            if (inside(tg.p)) {
                Triangle new_tr = Triangle(tg.p, res3.second, res1.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.p)) {
                Triangle new_tr = Triangle(res3.second, tg.q, tg.r);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.q, res3.second, res1.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        } else if (res3.first == NON_INTERSECTING) {
            if (inside(tg.q)) {
                Triangle new_tr = Triangle(tg.q, res1.second, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.q)) {
                Triangle new_tr = Triangle(tg.p, tg.r, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.p, res1.second, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        }
    } else {
        cout << "Should not be here 1" << endl;
    }
}

void clip_along_z_far(Triangle tg) {
    // clip by far plane
    plane far_plane(point(0, 0, -far), Vector(0, 0, 1));

    auto inside = [](point p) -> bool { return p.z > -far; };
    auto outside = [](point p) -> bool { return p.z < -far; };
    auto on_the_plane = [](point p) -> bool { return p.z == -far; };

    auto res1 = get_intersecting_point(tg.a, far_plane);
    auto res2 = get_intersecting_point(tg.b, far_plane);
    auto res3 = get_intersecting_point(tg.c, far_plane);

    int cut_values[] = {0, 0, 0, 0};
    cut_values[res1.first]++;
    cut_values[res2.first]++;
    cut_values[res3.first]++;

    if (cut_values[LYING] == 3) { // tg totally on the plane
        triangles.push_back(tg);

    } else if (cut_values[INTERSECTING] == 1 and cut_values[TOUCHING] == 2) { // tg 1 point out, 1 in, 1 on the plane
        point new_point;
        if (inside(tg.p)) {
            if (outside(tg.q)) {
                new_point = res1.second;
                triangles.emplace_back(tg.p, tg.r, new_point);
            } else {
                new_point = res3.second;
                triangles.emplace_back(tg.p, tg.q, new_point);
            }
        } else if (inside(tg.q)) {
            if (outside(tg.p)) {
                new_point = res1.second;
                triangles.emplace_back(tg.q, tg.r, new_point);
            } else {
                new_point = res2.second;
                triangles.emplace_back(tg.q, tg.p, new_point);
            }
        } else if (inside(tg.r)) {
            if (outside(tg.q)) {
                new_point = res2.second;
                triangles.emplace_back(tg.r, tg.p, new_point);
            } else {
                new_point = res3.second;
                triangles.emplace_back(tg.r, tg.q, new_point);
            }
        } else {
            cout << "Should not be here" << endl;
        }

        triangles.back().set_color(tg.colour);
    } else if (cut_values[NON_INTERSECTING] == 3) { // tg totally on 1 side
        if (inside(tg.p)) { // tg totally inside
            triangles.push_back(tg);
        }
    } else if (cut_values[NON_INTERSECTING] == 1 and cut_values[TOUCHING] == 2) { // 1 point on the plane 2 out/in plane
        if (on_the_plane(tg.p)) { // p on the plane
            if (inside(tg.q)) { // tg inside
                triangles.push_back(tg);
            }
        } else if (on_the_plane(tg.q)) { // q on the plane
            if (inside(tg.r)) { // tg inside
                triangles.push_back(tg);
            }
        }
        if (on_the_plane(tg.r)) { // r on the plane
            if (inside(tg.p)) { // tg inside
                triangles.push_back(tg);
            }
        }
    } else if (cut_values[TOUCHING] == 2 and cut_values[LYING] == 1) { // 2 on the plane 1 out/in
        if (inside(tg.p) or inside(tg.q) or inside(tg.r)) {
            triangles.push_back(tg);
        }
    } else if (cut_values[INTERSECTING] == 2 and cut_values[NON_INTERSECTING] == 1) { // 2 one side 1 different
        if (res1.first == NON_INTERSECTING) {
            if (inside(tg.r)) {
                Triangle new_tr = Triangle(tg.r, res2.second, res3.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.r)) {
                Triangle new_tr = Triangle(tg.p, tg.q, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.p, res2.second, res3.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        } else if (res2.first == NON_INTERSECTING) {
            if (inside(tg.p)) {
                Triangle new_tr = Triangle(tg.p, res3.second, res1.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.p)) {
                Triangle new_tr = Triangle(res3.second, tg.q, tg.r);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.q, res3.second, res1.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        } else if (res3.first == NON_INTERSECTING) {
            if (inside(tg.q)) {
                Triangle new_tr = Triangle(tg.q, res1.second, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            } else if (outside(tg.q)) {
                Triangle new_tr = Triangle(tg.p, tg.r, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);

                new_tr = Triangle(tg.p, res1.second, res2.second);
                new_tr.set_color(tg.colour);
                triangles.push_back(new_tr);
            }
        }
    } else {
        cout << "Should not be here 1" << endl;
    }
}

double interpolate_z_from_y(line l, double ys) {
    double z1 = l.p.z, z2 = l.q.z;
    double y1 = l.p.y, y2 = l.q.y;

    assert(y1 != y2);
    return z1 - (z1 - z2) / (y1 - y2) * (y1 - ys);
}

double interpolate_x_from_y(line l, double ys) {
    double x1 = l.p.x, x2 = l.q.x;
    double y1 = l.p.y, y2 = l.q.y;

    assert(y1 != y2);
    return x1 - (x1 - x2) / (y1 - y2) * (y1 - ys);
}

double interpolate_z_from_x(ppdd xz, double xp) {
    double xa = xz.first.first;
    double xb = xz.first.second;

    double za = xz.second.first;
    double zb = xz.second.second;

    assert(xa != xb);
    return zb - (zb - za) / (xb - xa) * (xb - xp);
}

ppdd get_x_and_z_from_y(point p, point q) {
    auto is_greater = [](double d1, double d2) -> bool { return (d1 - d2) > EPSILON; }; //if d1 greater than d2

    if (is_greater(p.x, q.x)) {
        // q.x < p.x
        return {{q.x, p.x},
                {q.z, p.z}};
    }
    return {{p.x, q.x},
            {p.z, q.z}};
}

ppdd get_x_and_z_from_y(line a, line b, double ys) {
    auto is_greater = [](double d1, double d2) -> bool { return (d1 - d2) > EPSILON; }; //if d1 greater than d2

    double xa, xb, za, zb;
    xa = interpolate_x_from_y(a, ys);
    xb = interpolate_x_from_y(b, ys);

    za = interpolate_z_from_y(a, ys);
    zb = interpolate_z_from_y(b, ys);

    if (is_greater(xa, xb)) {
        swap(xa, xb);
        swap(za, zb);
    }
    return {{xa, xb},
            {za, zb}};
}

void z_buffer_algorithm(vector<vector<color>> &pixels, vector<vector<double>> &zs) {
    double dy = 1.0 / screen_y;
    double dx = 1.0 / screen_x;

    auto get_index_from_y = [&dy](double y) -> int {
        return static_cast<int>(floor(((1.0 - y) / (2.0 * dy)) + EPSILON));
    };
    auto get_y_from_index = [&dy](int index) -> double { return 1.0 - (2.0 * index + 1) * dy; };

    auto get_index_from_x = [&dx](double x) -> int {
        return static_cast<int>(floor(((1.0 + x) / (2.0 * dx)) + EPSILON));
    };
    auto get_x_from_index = [&dx](int index) -> double { return 1.0 - (2.0 * index + 1) * dx; };

    auto is_same_value = [](double d1, double d2) -> bool { return abs(d1 - d2) < EPSILON; };
    auto is_greater = [](double d1, double d2) -> bool { return (d1 - d2) > EPSILON; }; //if d1 greater than d2
    auto is_less = [](double d1, double d2) -> bool { return (d2 - d1) > EPSILON; };

//    shuffle(begin(triangles), end(triangles), default_random_engine(20));
    for (auto &triangle: triangles) {
        double high_y = triangle.p.y;
        double mid_y = triangle.q.y;
        double low_y = triangle.r.y;
//        bool flag1 = true, flag2 = true;

        int low_y_index = get_index_from_y(high_y); // high y value means low index
        int high_y_index = get_index_from_y(low_y);

        int ly = max(0, low_y_index);
        int hy = min(high_y_index, screen_y - 1);
        for (int y_index = ly; y_index <= hy; y_index++) {   // TODO: check boundary
            double ys = get_y_from_index(y_index);
            if (ys > 1 or ys < -1) {
                cout << "Should never be printed Y\n";
                continue;
            }
            ppdd xz;
            auto za = 2.0, zb = 2.0;
            auto xa = 2.0, xb = 2.0;  // xa < xb must
            pair<line, line> lines_to_intersect;

            if (is_same_value(high_y, mid_y)) {
                lines_to_intersect = {triangle.b, triangle.c};
            } else if (is_same_value(mid_y, low_y)) {
                lines_to_intersect = {triangle.a, triangle.c};
            } else if (is_greater(ys, mid_y)) {
                lines_to_intersect = {triangle.a, triangle.c};
            } else {
                lines_to_intersect = {triangle.b, triangle.c};
            };
            /*if (ys >= mid_y) {
                // here we are considering triangles line (a and c)
                if (is_same_value(high_y, mid_y) and flag1) {
                    xz = get_x_and_z_from_y(triangle.p, triangle.q);
                    flag1 = false;
                } else {
                    xz = get_x_and_z_from_y(triangle.a, triangle.c, ys);
                }
            } else {
                // here we are considering triangles line (b and c)
                if (is_same_value(mid_y, low_y) and flag2) {
                    // q.x < r.x
                    xz = get_x_and_z_from_y(triangle.q, triangle.r);
                    flag2 = false;
                } else {
                    xz = get_x_and_z_from_y(triangle.b, triangle.c, ys);
                }
            }*/
            xz = get_x_and_z_from_y(lines_to_intersect.first, lines_to_intersect.second, ys);
            xa = xz.first.first;
            xb = xz.first.second;

            za = xz.second.first;
            zb = xz.second.second;

            if (is_same_value(xa, xb)) {
                cout << "One point to draw\n";
            }

            int low_x_index = get_index_from_x(xa); // x value low means index high
            int high_x_index = get_index_from_x(xb);

            int lx = max(0, low_x_index);
            int hx = min(high_x_index, screen_x - 1);

            for (int x_index = lx; x_index <= hx; x_index++) { // TODO: check boundary
                double xp = get_x_from_index(x_index);
                if (xp > 1 or xp < -1) {
                    cout << "Should never be printed X\n";
                    continue;
                }
                double zp = interpolate_z_from_x(xz, xp);

                if (is_less(zp, zs[y_index][x_index])) {
                    zs[y_index][x_index] = zp;
                    pixels[y_index][x_index] = triangle.colour;
                }
            }
        }
    }
}

void scan_convert() {
    /*ifstream stage3;
    stage3.open("stage3.txt");*/

    vector<vector<color>> pixels(screen_y, vector<color>(screen_x, background));
    vector<vector<double>> zs(screen_y, vector<double>(screen_x, +20));

    // perform scan conversion, populate the 2D array pixels
    // the array zs is the z-buffer.
    z_buffer_algorithm(pixels, zs);

    // the following code generates a bmp image. do not change this.
    bitmap_image image(screen_x, screen_y);
    for (unsigned y = 0; y < screen_y; y++) {
        for (unsigned x = 0; x < screen_x; x++) {
            image.set_pixel(x, y, pixels[y][x].r, pixels[y][x].g, pixels[y][x].b);
        }
    }
    image.save_image("image.bmp");
}

void stage3() {
    if (near == far)
        return;
    ifstream stage2;
    ofstream stage3;
    stage2.open("stage2.txt");
    stage3.open("stage3.txt");
    if (!stage2.is_open() or !stage3.is_open()) {
        cerr << "Stage 3 file(s) not found\n";
        exit(EXIT_FAILURE);
    }
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    // process input from stage2 and write to stage3
    fov_x = fov_y * aspect_ratio;
    double t = near * tan(fov_y / 2 * PI / 180);
    double r = near * tan(fov_x / 2 * PI / 180);

    matrix P = matrix::make_identity(4);
    P[0][0] = near / r;
    P[1][1] = near / t;
    P[2][2] = -(far + near) / (far - near);
    P[2][3] = -(2 * far * near) / (far - near);
    P[3][2] = -1;
    P[3][3] = 0;
    // P.print();
    vector<Triangle> temp_triangles;
    for (const auto &color : colors) {
        Triangle T;
        stage2 >> T.p.x >> T.p.y >> T.p.z;
        stage2 >> T.q.x >> T.q.y >> T.q.z;
        stage2 >> T.r.x >> T.r.y >> T.r.z;
        T.set_sides();
        T.set_color(color);
        temp_triangles.push_back(T);
        // T.print();
    }
    for_each(temp_triangles.begin(), temp_triangles.end(), clip_along_z_near);
    temp_triangles.clear();

    temp_triangles = move(triangles);
    for_each(temp_triangles.begin(), temp_triangles.end(), clip_along_z_far);

    for (auto &triangle : triangles) {
        point _p = P * triangle.p;
        stage3 << _p.x << " " << _p.y << " " << _p.z << " " /* << _p.w  */ << "\n";
        triangle.p = _p;

        _p = P * triangle.q;
        stage3 << _p.x << " " << _p.y << " " << _p.z << " " /* << _p.w  */ << "\n";
        triangle.q = _p;


        _p = P * triangle.r;
        stage3 << _p.x << " " << _p.y << " " << _p.z << " " /* << _p.w  */ << "\n";
        triangle.r = _p;

        stage3 << endl;
        triangle.sort_points_wrt_y();
        triangle.set_sides();
    }

    stage3.close();
    stage2.close();
}

void stage2() {
    ifstream stage1;
    ofstream stage2;
    stage1.open("stage1.txt");
    stage2.open("stage2.txt");
    if (!stage1.is_open() or !stage2.is_open()) {
        cerr << "Stage 2 file(s) not found\n";
        exit(EXIT_FAILURE);
    }
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    // collect input from stage1 and process, write output to stage2
    Vector l(look_x - eye_x, look_y - eye_y, look_z - eye_z);
    l.normalize();
    Vector r = Vector::cross(l, Vector(up_x, up_y, up_z));
    r.normalize();
    Vector u = Vector::cross(r, l);

    matrix T = matrix::make_identity(4);
    T[0][3] = -eye_x;
    T[1][3] = -eye_y;
    T[2][3] = -eye_z;
    // T.print();
    matrix R = matrix::make_identity(4);
    R[0][0] = r.x, R[0][1] = r.y, R[0][2] = r.z;
    R[1][0] = u.x, R[1][1] = u.y, R[1][2] = u.z;
    R[2][0] = -l.x, R[2][1] = -l.y, R[2][2] = -l.z;
    // R.print();
    matrix V = R * T;
    // V.print();

    for (int i = 0; i < colors.size(); i++) {
        for (int j = 0; j < 3; j++) {
            point p;
            stage1 >> p.x >> p.y >> p.z;
            point _p = V * p;
            stage2 << _p.x << " " << _p.y << " " << _p.z << " " /* << _p.w  */ << "\n";
        }
        stage2 << endl;
    }

    stage1.close();
    stage2.close();
}

void stage1() {
    ifstream scene;
    ofstream stage1;
    scene.open("scene.txt");
    stage1.open("stage1.txt");
    if (!scene.is_open() or !stage1.is_open()) {
        cerr << "Stage 1 file(s) not found\n";
        exit(EXIT_FAILURE);
    }
    stage1 << std::fixed;
    stage1 << std::setprecision(7);

    string command;

    scene >> eye_x >> eye_y >> eye_z;
    scene >> look_x >> look_y >> look_z;
    scene >> up_x >> up_y >> up_z;
    scene >> fov_y >> aspect_ratio >> near >> far;
    scene >> screen_x >> screen_y;
    scene >> background.r >> background.g >> background.b;

    stack<matrix> Stack;
    Stack.push(matrix::make_identity(4));
    // take other commands as input from scene in a loop
    // process accordingly
    // write to stage1
    while (true) {
        scene >> command;
        if (command == "end")
            break;
        else if (command == "triangle") {
            // Stack.top().print();
            for (int i = 0; i < 3; i++) {
                point p;
                scene >> p.x >> p.y >> p.z;
                point _p = Stack.top() * p;
                stage1 << _p.x << " " << _p.y << " " << _p.z << " " /* << _p.w  */ << "\n";
            }
            stage1 << endl;

            color triangle_color{};
            scene >> triangle_color.r >> triangle_color.g >> triangle_color.b;
            colors.push_back(triangle_color);
        } else if (command == "translate") {
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            matrix translation_matrix = matrix::make_identity(4);
            translation_matrix[0][3] = tx;
            translation_matrix[1][3] = ty;
            translation_matrix[2][3] = tz;
            // translation_matrix.print();
            matrix new_transformation_mat = translation_matrix * Stack.top();
            Stack.pop();
            Stack.push(new_transformation_mat);
        } else if (command == "scale") {
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            matrix scaling_matrix = matrix::make_identity(4);
            scaling_matrix[0][0] = sx;
            scaling_matrix[1][1] = sy;
            scaling_matrix[2][2] = sz;
            // scaling_matrix.print();
            matrix new_transformation_mat = scaling_matrix * Stack.top();
            Stack.pop();
            Stack.push(new_transformation_mat);
        } else if (command == "rotate") {
            double deg_angle, ax, ay, az;
            scene >> deg_angle >> ax >> ay >> az;
            double angle = deg_angle * PI / 180.0;
            Vector a(ax, ay, az);
            a.normalize();
            Vector i(1, 0, 0), j(0, 1, 0), k(0, 0, 1);
            Vector C1 =
                    i * cos(angle) + (a * (Vector::dot(a, i))) * (1 - cos(angle)) + Vector::cross(a, i) * sin(angle);
            Vector C2 =
                    j * cos(angle) + (a * (Vector::dot(a, j))) * (1 - cos(angle)) + Vector::cross(a, j) * sin(angle);
            Vector C3 =
                    k * cos(angle) + (a * (Vector::dot(a, k))) * (1 - cos(angle)) + Vector::cross(a, k) * sin(angle);
            matrix rotation_matrix(4);
            rotation_matrix[0][0] = C1.x, rotation_matrix[1][0] = C1.y, rotation_matrix[2][0] = C1.z, rotation_matrix[3][0] = 0;
            rotation_matrix[0][1] = C2.x, rotation_matrix[1][1] = C2.y, rotation_matrix[2][1] = C2.z, rotation_matrix[3][1] = 0;
            rotation_matrix[0][2] = C3.x, rotation_matrix[1][2] = C3.y, rotation_matrix[2][2] = C3.z, rotation_matrix[3][2] = 0;
            rotation_matrix[0][3] = 0, rotation_matrix[1][3] = 0, rotation_matrix[2][3] = 0, rotation_matrix[3][3] = 1;
            // cout << "Rotation Matrix: \n";
            // rotation_matrix.print();
            matrix new_transformation_mat = rotation_matrix * Stack.top();
            Stack.pop();
            Stack.push(new_transformation_mat);
        } else if (command == "push") {
            Stack.push(Stack.top());
        } else if (command == "pop") {
            if (Stack.size() > 1) {
                Stack.pop();
            } else {
                cerr << "Push failed" << endl;
            }
        } else {
            cerr << "Unrecognized command" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // cout << colors.size() << endl;

    scene.close();
    stage1.close();
}

int main() {
    cout << std::fixed;
    cout << std::setprecision(4);

    stage1();
    stage2();
    stage3();
    scan_convert();

    /*point p(1, 2, 4), q(5, 20, 40), r(-1, -22, 4);
    line l(p, q), m(q, r);
    double ys = 2.0;
    while (true) {
        cin >> ys;
        cout << interpolate_z_from_y(l, ys) << endl;
        cout << interpolate_x_from_y(l, ys) << endl;
        cout << interpolate_z_from_y(m, ys) << endl;
        cout << interpolate_x_from_y(m, ys) << endl;
    }*/

    /*screen_x = 500;
    screen_y = 500;

    double del_y = 1.0 / screen_y;
    double del_x = 1.0 / screen_x;

    auto get_index_from_y = [&del_y](double y) -> int {
        return static_cast<int>(floor(((1.0 - y) / (2.0 * del_y)) + EPSILON));
    };
    auto get_y_from_index = [&del_y](int index) -> double { return 1.0 - (2.0 * index + 1) * del_y; };

    auto get_index_from_x = [&del_x](double x) -> int {
        return static_cast<int>(floor(((1.0 + x) / (2.0 * del_x)) + EPSILON));
    };
    auto get_x_from_index = [&del_x](int index) -> double { return 1.0 - (2.0 * index + 1) * del_x; };

    auto is_same_value = [](double d1, double d2) -> bool { return abs(d1 - d2) < EPSILON; };
    auto is_greater = [](double d1, double d2) -> bool { return (d1 - d2) > EPSILON; }; //if d1 greater than d2
    auto is_less = [](double d1, double d2) -> bool { return (d2 - d1) > EPSILON; };    // if d1 less than d2

    double y = -2.0;
    double x = 2.0;
    int count = 0;
    for (int i = 0; i < 8000; ++i) {
        int i1 = get_index_from_y(y), i2 = get_index_from_x(x);
        cout << i1 << " " << i2 << endl;
//        if (i1 + i2 != 499)
//            count++;
        y += 2.0 / 2000;
        x -= 2.0 / 2000;
    }
    cout << count << endl; */

    /*auto x = 1.0, y = 2.0;

    for (int i = 0; i < 10; ++i) {
        cout << x << " " << y << endl;
        if (is_greater(y, x))
            cout << "greater" << endl;
        else
            cout << "Not" << endl;
        x /= 10;
        y /= 10;
    }*/

    /* point p(0, 0, 67), q(6, 6, 67);
    point near_p(0, 0, 67), far_p(0, 0, 6);
    plane pl(near_p, Vector(0, 0, 1));
    auto x = get_intersecting_point(line(p, q), pl);
    cout << x.first << "\n"; 
    x.second.print(); */
    return 0;
}
