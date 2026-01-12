#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

const double pi = 3.14159265358979323846;

// Compute shader source code
const char* computeShaderSource = R"(
#version 430 core
#extension GL_ARB_gpu_shader_fp64 : enable

// Constants for double precision
const double PI_D = 3.14159265358979323846LF;
const double TWO_PI_D = 6.28318530717958647692LF;
const double HALF_PI_D = 1.57079632679489661923LF;

// Range reduction: bring angle to [-PI, PI]
double reduceAngle(double x) {
    x = mod(x + PI_D, TWO_PI_D);
    if (x < 0.0LF)
        x += TWO_PI_D;
    return x - PI_D;
}

// Double-precision sine using Taylor series
double sin(double x) {
    x = reduceAngle(x);
    
    double x2 = x * x;
    double result = x;
    double term = x;
    
    // Taylor series
    term *= -x2 / (2.0LF * 3.0LF);
    result += term;
    
    term *= -x2 / (4.0LF * 5.0LF);
    result += term;
    
    term *= -x2 / (6.0LF * 7.0LF);
    result += term;
    
    term *= -x2 / (8.0LF * 9.0LF);
    result += term;
    
    term *= -x2 / (10.0LF * 11.0LF);
    result += term;
    
    term *= -x2 / (12.0LF * 13.0LF);
    result += term;
    
    term *= -x2 / (14.0LF * 15.0LF);
    result += term;
    
    term *= -x2 / (16.0LF * 17.0LF);
    result += term;
    
    return result;
}

// Double-precision cosine using Taylor series
double cos(double x) {
    x = reduceAngle(x);
    
    double x2 = x * x;
    double result = 1.0LF;
    double term = 1.0LF;
    
    // Taylor series.
    term *= -x2 / (1.0LF * 2.0LF);
    result += term;
    
    term *= -x2 / (3.0LF * 4.0LF);
    result += term;
    
    term *= -x2 / (5.0LF * 6.0LF);
    result += term;
    
    term *= -x2 / (7.0LF * 8.0LF);
    result += term;
    
    term *= -x2 / (9.0LF * 10.0LF);
    result += term;
    
    term *= -x2 / (11.0LF * 12.0LF);
    result += term;
    
    term *= -x2 / (13.0LF * 14.0LF);
    result += term;
    
    term *= -x2 / (15.0LF * 16.0LF);
    result += term;
    
    return result;
}

// Convenience: compute both at once (more efficient)
void sincos_d(double x, out double s, out double c) {
    x = reduceAngle(x);
    
    double x2 = x * x;
    
    // Sine
    double s_result = x;
    double s_term = x;
    
    s_term *= -x2 / 6.0LF;    s_result += s_term;
    s_term *= -x2 / 20.0LF;   s_result += s_term;
    s_term *= -x2 / 42.0LF;   s_result += s_term;
    s_term *= -x2 / 72.0LF;   s_result += s_term;
    s_term *= -x2 / 110.0LF;  s_result += s_term;
    s_term *= -x2 / 156.0LF;  s_result += s_term;
    s_term *= -x2 / 210.0LF;  s_result += s_term;
    s_term *= -x2 / 272.0LF;  s_result += s_term;
    
    // Cosine
    double c_result = 1.0LF;
    double c_term = 1.0LF;
    
    c_term *= -x2 / 2.0LF;    c_result += c_term;
    c_term *= -x2 / 12.0LF;   c_result += c_term;
    c_term *= -x2 / 30.0LF;   c_result += c_term;
    c_term *= -x2 / 56.0LF;   c_result += c_term;
    c_term *= -x2 / 90.0LF;   c_result += c_term;
    c_term *= -x2 / 132.0LF;  c_result += c_term;
    c_term *= -x2 / 182.0LF;  c_result += c_term;
    c_term *= -x2 / 240.0LF;  c_result += c_term;
    
    s = s_result;
    c = c_result;
}



layout(local_size_x = 256) in;

// Input uniforms
uniform double emitter_radius;
uniform double receiver_distance;
uniform double receiver_distance_plus;
uniform double receiver_radius;
uniform uint total_samples;
uniform uint seed_offset;

// Output buffer for partial sums
layout(std430, binding = 0) buffer ResultBuffer {
    double partial_sums[];
};

layout(std430, binding = 1) buffer ResultBufferPlus {
    double partial_sums_plus[];
};

const double PI = 3.14159265358979323846LF;

// PCG random number generator state
uint pcg_state;

void pcg_init(uint seed) {
    pcg_state = seed;
}

uint pcg_next() {
    uint oldstate = pcg_state;
    pcg_state = oldstate * 747796405u + 2891336453u;
    uint word = ((oldstate >> ((oldstate >> 28u) + 4u)) ^ oldstate) * 277803737u;
    return (word >> 22u) ^ word;
}

double random_double() {
    return double(pcg_next()) / double(4294967295.0);
}

dvec3 normalize_dvec3(dvec3 v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len > 0.0LF) {
        return v / len;
    }
    return v;
}

double length_dvec3(dvec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double dot_dvec3(dvec3 a, dvec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

dvec3 cross_dvec3(dvec3 a, dvec3 b) {
    return dvec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

dvec3 random_unit_vector() {
    double z = random_double() * 2.0LF - 1.0LF;
    double a = random_double() * 2.0LF * PI;
    
    double r = sqrt(1.0LF - z * z);
    double x = r * cos(a);
    double y = r * sin(a);
    
    return normalize_dvec3(dvec3(x, y, z));
}

dvec3 random_cosine_weighted_hemisphere(dvec3 normal) {
    double u1 = random_double();
    double u2 = random_double();
    
    double r = sqrt(u1);
    double theta = 2.0LF * PI * u2;
    
    double x = r * cos(theta);
    double y = r * sin(theta);
    double z = sqrt(1.0LF - u1);
    
    dvec3 n = normalize_dvec3(normal);
    
    dvec3 arbitrary;
    if (abs(n.x) > 0.9LF)
        arbitrary = dvec3(0.0LF, 1.0LF, 0.0LF);
    else
        arbitrary = dvec3(1.0LF, 0.0LF, 0.0LF);
    
    dvec3 tangent = normalize_dvec3(cross_dvec3(n, arbitrary));
    dvec3 bitangent = normalize_dvec3(cross_dvec3(n, tangent));
    
    dvec3 result;
    result.x = tangent.x * x + bitangent.x * y + n.x * z;
    result.y = tangent.y * x + bitangent.y * y + n.y * z;
    result.z = tangent.z * x + bitangent.z * y + n.z * z;
    
    return normalize_dvec3(result);
}

double intersect_AABB(dvec3 min_location, dvec3 max_location, dvec3 ray_origin, dvec3 ray_dir) {
    double tmin = (min_location.x - ray_origin.x) / ray_dir.x;
    double tmax = (max_location.x - ray_origin.x) / ray_dir.x;
    
    if (tmin > tmax) {
        double temp = tmin;
        tmin = tmax;
        tmax = temp;
    }
    
    double tymin = (min_location.y - ray_origin.y) / ray_dir.y;
    double tymax = (max_location.y - ray_origin.y) / ray_dir.y;
    
    if (tymin > tymax) {
        double temp = tymin;
        tymin = tymax;
        tymax = temp;
    }
    
    if ((tmin > tymax) || (tymin > tmax))
        return 0.0LF;
    
    if (tymin > tmin)
        tmin = tymin;
    
    if (tymax < tmax)
        tmax = tymax;
    
    double tzmin = (min_location.z - ray_origin.z) / ray_dir.z;
    double tzmax = (max_location.z - ray_origin.z) / ray_dir.z;
    
    if (tzmin > tzmax) {
        double temp = tzmin;
        tzmin = tzmax;
        tzmax = temp;
    }
    
    if ((tmin > tzmax) || (tzmin > tmax))
        return 0.0LF;
    
    if (tzmin > tmin)
        tmin = tzmin;
    
    if (tzmax < tmax)
        tmax = tzmax;
    
    if (tmin < 0.0LF || tmax < 0.0LF)
        return 0.0LF;
    
    dvec3 ray_hit_start = ray_origin + ray_dir * tmin;
    dvec3 ray_hit_end = ray_origin + ray_dir * tmax;
    
    return length_dvec3(ray_hit_end - ray_hit_start);
}

double intersect(dvec3 location, dvec3 normal, double recv_distance, double recv_radius) {
    dvec3 circle_origin = dvec3(recv_distance, 0.0LF, 0.0LF);
    
    if (dot_dvec3(normal, circle_origin) <= 0.0LF)
        return 0.0LF;
    
    dvec3 min_location = dvec3(-recv_radius + recv_distance, -recv_radius, -recv_radius);
    dvec3 max_location = dvec3(recv_radius + recv_distance, recv_radius, recv_radius);
    
    return intersect_AABB(min_location, max_location, location, normal);
}

void main() {
    uint gid = gl_GlobalInvocationID.x;
    
    if (gid >= total_samples)
        return;
    
    // Initialize RNG with unique seed per thread
    pcg_init(gid + seed_offset * 1000000u);
    
    // Generate random point on sphere
    dvec3 location = random_unit_vector();
    location *= emitter_radius;
    
    // Surface normal at that point
    dvec3 surface_normal = normalize_dvec3(location);
    
    // Random direction in hemisphere
    dvec3 normal = random_cosine_weighted_hemisphere(surface_normal);
    
    // Compute intersections
    double count = intersect(location, normal, receiver_distance, receiver_radius);
    double count_plus = intersect(location, normal, receiver_distance_plus, receiver_radius);
    
    // Store results
    partial_sums[gid] = count;
    partial_sums_plus[gid] = count_plus;
}
)";

// Reduction compute shader
const char* reductionShaderSource = R"(
#version 430 core
#extension GL_ARB_gpu_shader_fp64 : enable

layout(local_size_x = 256) in;

uniform uint array_size;
uniform uint stride;

layout(std430, binding = 0) buffer DataBuffer {
    double data[];
};

shared double shared_data[256];

void main() {
    uint gid = gl_GlobalInvocationID.x;
    uint lid = gl_LocalInvocationID.x;
    
    // Load data into shared memory
    if (gid < array_size) {
        shared_data[lid] = data[gid];
    } else {
        shared_data[lid] = 0.0LF;
    }
    
    barrier();
    
    // Parallel reduction in shared memory
    for (uint s = gl_WorkGroupSize.x / 2u; s > 0u; s >>= 1u) {
        if (lid < s && gid + s < array_size) {
            shared_data[lid] += shared_data[lid + s];
        }
        barrier();
    }
    
    // Write result
    if (lid == 0u) {
        data[gl_WorkGroupID.x] = shared_data[0];
    }
}
)";

GLuint computeProgram;
GLuint reductionProgram;
GLuint resultBuffer;
GLuint resultBufferPlus;

void checkShaderCompilation(GLuint shader, const string& name) {
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar infoLog[1024];
        glGetShaderInfoLog(shader, 1024, NULL, infoLog);
        cerr << "Shader compilation error (" << name << "):\n" << infoLog << endl;
        exit(1);
    }
}

void checkProgramLinking(GLuint program, const string& name) {
    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        GLchar infoLog[1024];
        glGetProgramInfoLog(program, 1024, NULL, infoLog);
        cerr << "Program linking error (" << name << "):\n" << infoLog << endl;
        exit(1);
    }
}

GLuint createComputeProgram(const char* source, const string& name) {
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);
    checkShaderCompilation(shader, name);

    GLuint program = glCreateProgram();
    glAttachShader(program, shader);
    glLinkProgram(program);
    checkProgramLinking(program, name);

    glDeleteShader(shader);
    return program;
}

double reduceBuffer(GLuint buffer, GLuint size) {
    glUseProgram(reductionProgram);

    GLuint currentSize = size;

    while (currentSize > 1) {
        GLuint numGroups = (currentSize + 255) / 256;

        glUniform1ui(glGetUniformLocation(reductionProgram, "array_size"), currentSize);
        glUniform1ui(glGetUniformLocation(reductionProgram, "stride"), 1);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, buffer);
        glDispatchCompute(numGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        currentSize = numGroups;
    }

    // Read back single result
    double result;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(double), &result);

    return result;
}

double get_intersecting_line_density_gpu(
    unsigned long long n,
    double emitter_radius,
    double receiver_distance,
    double receiver_distance_plus,
    double receiver_radius)
{
    GLuint numSamples = static_cast<GLuint>(n);
    GLuint numGroups = (numSamples + 255) / 256;

    // Create/resize buffers
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, resultBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, numSamples * sizeof(double), NULL, GL_DYNAMIC_COPY);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, resultBufferPlus);
    glBufferData(GL_SHADER_STORAGE_BUFFER, numSamples * sizeof(double), NULL, GL_DYNAMIC_COPY);

    // Set uniforms and dispatch compute shader
    glUseProgram(computeProgram);

    glUniform1d(glGetUniformLocation(computeProgram, "emitter_radius"), emitter_radius);
    glUniform1d(glGetUniformLocation(computeProgram, "receiver_distance"), receiver_distance);
    glUniform1d(glGetUniformLocation(computeProgram, "receiver_distance_plus"), receiver_distance_plus);
    glUniform1d(glGetUniformLocation(computeProgram, "receiver_radius"), receiver_radius);
    glUniform1ui(glGetUniformLocation(computeProgram, "total_samples"), numSamples);
    glUniform1ui(glGetUniformLocation(computeProgram, "seed_offset"), 0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, resultBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, resultBufferPlus);

    glDispatchCompute(numGroups, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    // Reduce results
    double count = reduceBuffer(resultBuffer, numSamples);

    // Re-upload data for plus buffer reduction (since we modified the buffer)
    // Actually we need to copy the plus buffer data first
    vector<double> plusData(numSamples);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, resultBufferPlus);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, numSamples * sizeof(double), plusData.data());

    // Use resultBuffer for reduction of plus data
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, resultBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, numSamples * sizeof(double), plusData.data(), GL_DYNAMIC_COPY);

    double count_plus = reduceBuffer(resultBuffer, numSamples);

    return count_plus - count;
}

void initGL() {
    // Initialize GLEW
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        cerr << "GLEW initialization failed: " << glewGetErrorString(err) << endl;
        exit(1);
    }

    // Check for compute shader support
    if (!GLEW_ARB_compute_shader) {
        cerr << "Compute shaders not supported!" << endl;
        exit(1);
    }

    // Check for double precision support
    if (!GLEW_ARB_gpu_shader_fp64) {
        cerr << "Double precision not supported!" << endl;
        exit(1);
    }

    cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
    cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

    // Create compute programs
    computeProgram = createComputeProgram(computeShaderSource, "main compute");
    reductionProgram = createComputeProgram(reductionShaderSource, "reduction");

    // Create buffers
    glGenBuffers(1, &resultBuffer);
    glGenBuffers(1, &resultBufferPlus);
}

void runSimulation() {
    ofstream outfile("ratio");

    const double emitter_radius_geometrized = sqrt(1e11 * log(2.0) / pi);
    const double receiver_radius_geometrized = emitter_radius_geometrized * 0.01;
    const double emitter_area_geometrized = 4.0 * pi * emitter_radius_geometrized * emitter_radius_geometrized;
    const double n_geometrized = emitter_area_geometrized / (log(2.0) * 4.0);
    const double emitter_mass_geometrized = emitter_radius_geometrized / 2.0;

    double start_pos = emitter_radius_geometrized + receiver_radius_geometrized;
    double end_pos = start_pos * 10;

    const size_t pos_res = 10;
    const double pos_step_size = (end_pos - start_pos) / (pos_res - 1);
    const double epsilon = receiver_radius_geometrized;

    cout << "Emitter radius: " << emitter_radius_geometrized << endl;
    cout << "Receiver radius: " << receiver_radius_geometrized << endl;
    cout << "Number of samples: " << static_cast<unsigned long long>(n_geometrized) << endl;
    cout << endl;

    for (size_t i = 0; i < pos_res; i++) {
        const double receiver_distance_geometrized = start_pos + i * pos_step_size;
        const double receiver_distance_plus_geometrized = receiver_distance_geometrized + epsilon;

        cout << "Processing step " << (i + 1) << "/" << pos_res << endl;
        cout << "Receiver distance: " << receiver_distance_geometrized << endl;

        const double collision_count_plus_minus_collision_count =
            get_intersecting_line_density_gpu(
                static_cast<unsigned long long>(n_geometrized),
                emitter_radius_geometrized,
                receiver_distance_geometrized,
                receiver_distance_plus_geometrized,
                receiver_radius_geometrized);

        const double gradient_integer = collision_count_plus_minus_collision_count / epsilon;

        double gradient_strength = -gradient_integer /
            (2.0 * receiver_radius_geometrized *
                receiver_radius_geometrized *
                receiver_radius_geometrized);

        const double a_Newton_geometrized = sqrt(
            n_geometrized * log(2.0) /
            (4.0 * pi * pow(receiver_distance_geometrized, 4.0)));

        const double a_flat_geometrized =
            gradient_strength * receiver_distance_geometrized * log(2.0) /
            (8.0 * emitter_mass_geometrized);

        const double dt_Schwarzschild = sqrt(1 - emitter_radius_geometrized / receiver_distance_geometrized);

        const double a_Schwarzschild_geometrized =
            emitter_radius_geometrized /
            (pi * pow(receiver_distance_geometrized, 2.0) * dt_Schwarzschild);

        cout << "a_Schwarzschild_geometrized " << a_Schwarzschild_geometrized << endl;
        cout << "a_Newton_geometrized " << a_Newton_geometrized << endl;
        cout << "a_flat_geometrized " << a_flat_geometrized << endl;
        cout << a_Schwarzschild_geometrized / a_flat_geometrized << endl;
        cout << endl;
        cout << a_Newton_geometrized / a_flat_geometrized << endl;
        cout << endl << endl;

        outfile << receiver_distance_geometrized << " "
            << (a_Schwarzschild_geometrized / a_flat_geometrized) << endl;
    }

    outfile.close();
    cout << "Results written to 'ratio' file." << endl;
}

void cleanup() {
    glDeleteProgram(computeProgram);
    glDeleteProgram(reductionProgram);
    glDeleteBuffers(1, &resultBuffer);
    glDeleteBuffers(1, &resultBufferPlus);
}

void display() {
    // Empty display callback required by GLUT
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
}

void idle() {
    static bool hasRun = false;
    if (!hasRun) {
        hasRun = true;
        runSimulation();
        cleanup();
        exit(0);
    }
}


#pragma comment(lib, "freeglut")
#pragma comment(lib, "glew32")



int main(int argc, char** argv) {
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(100, 100);
    glutCreateWindow("Compute Shader Simulation");

    // Initialize OpenGL
    initGL();

    // Set callbacks
    glutDisplayFunc(display);
    glutIdleFunc(idle);

    // Run main loop
    glutMainLoop();

    return 0;
}