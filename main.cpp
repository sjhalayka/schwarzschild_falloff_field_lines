#include <glm/glm.hpp>

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;

#include <utility>
using std::swap;

#include <cmath>

#include <random>


std::mt19937 generator(0);
std::uniform_real_distribution<double> dis(0.0, 1.0);

const double pi = 4.0 * atan(1.0);


double intersect_AABB(const glm::dvec3 min_location, const glm::dvec3 max_location, const glm::dvec3 ray_origin, const glm::dvec3 ray_dir, double& tmin, double& tmax)
{
	tmin = (min_location.x - ray_origin.x) / ray_dir.x;
	tmax = (max_location.x - ray_origin.x) / ray_dir.x;

	if (tmin > tmax)
		swap(tmin, tmax);

	double tymin = (min_location.y - ray_origin.y) / ray_dir.y;
	double tymax = (max_location.y - ray_origin.y) / ray_dir.y;

	if (tymin > tymax)
		swap(tymin, tymax);

	if ((tmin > tymax) || (tymin > tmax))
		return 0;

	if (tymin > tmin)
		tmin = tymin;

	if (tymax < tmax)
		tmax = tymax;

	double tzmin = (min_location.z - ray_origin.z) / ray_dir.z;
	double tzmax = (max_location.z - ray_origin.z) / ray_dir.z;

	if (tzmin > tzmax)
		swap(tzmin, tzmax);

	if ((tmin > tzmax) || (tzmin > tmax))
		return 0;

	if (tzmin > tmin)
		tmin = tzmin;

	if (tzmax < tmax)
		tmax = tzmax;

	if (tmin < 0 || tmax < 0)
		return 0;

	glm::dvec3 ray_hit_start = ray_origin;
	ray_hit_start.x += ray_dir.x * tmin;
	ray_hit_start.y += ray_dir.y * tmin;
	ray_hit_start.z += ray_dir.z * tmin;

	glm::dvec3 ray_hit_end = ray_origin;
	ray_hit_end.x += ray_dir.x * tmax;
	ray_hit_end.y += ray_dir.y * tmax;
	ray_hit_end.z += ray_dir.z * tmax;

	double l = glm::length(ray_hit_end - ray_hit_start);

	return l;
}

double intersect(
	const glm::dvec3 location,
	const glm::dvec3 normal,
	const double receiver_distance,
	const double receiver_radius)
{
	const glm::dvec3 circle_origin(receiver_distance, 0, 0);

	if (glm::dot(normal, circle_origin) <= 0)
		return 0.0;

	glm::dvec3 min_location(-receiver_radius + receiver_distance, -receiver_radius, -receiver_radius);
	glm::dvec3 max_location(receiver_radius + receiver_distance, receiver_radius, receiver_radius);

	double tmin = 0, tmax = 0;

	return intersect_AABB(min_location, max_location, location, normal, tmin, tmax);
}

glm::dvec3 random_cosine_weighted_hemisphere(const glm::dvec3& normal)
{
	// Method 1:
	glm::dvec2 r = glm::vec2(dis(generator), dis(generator));
	glm::dvec3 uu = glm::normalize(glm::cross(normal, glm::dvec3(0.0, 1.0, 1.0)));
	glm::dvec3 vv = glm::cross(uu, normal);

	double ra = sqrt(r.y);
	double rx = ra * cos(2.0 * pi * r.x);
	double ry = ra * sin(2.0 * pi * r.x);
	double rz = sqrt(1.0 - r.y);
	glm::dvec3 rr = glm::dvec3(rx * uu + ry * vv + rz * normal);

	return normalize(rr);


	// Method 2:
	//double u1 = dis(generator);
	//double u2 = dis(generator);

	//double r = sqrt(u1);
	//double theta = 2.0 * pi * u2;

	//double x = r * cos(theta);
	//double y = r * sin(theta);
	//double z = sqrt(1.0 - u1);

	//glm::dvec3 n = normalize(normal);

	//glm::dvec3 arbitrary;
	//if (fabs(n.x) > 0.9)
	//	arbitrary = glm::dvec3(0, 1, 0);
	//else
	//	arbitrary = glm::dvec3(1, 0, 0);

	//glm::dvec3 tangent = glm::normalize(glm::cross(n, arbitrary));
	////tangent.normalize();

	//glm::dvec3 bitangent = glm::normalize(glm::cross(n, tangent));
	////bitangent.normalize();

	//glm::dvec3 result;
	//result.x = 
	//	tangent.x * x +
	//	bitangent.x * y +
	//	n.x * z;

	//result.y = 
	//	tangent.y * x +
	//	bitangent.y * y +
	//	n.y * z;

	//result.z = 
	//	tangent.z * x +
	//	bitangent.z * y +
	//	n.z * z;

	//return glm::normalize(result);
}

glm::dvec3 random_unit_vector(void)
{
	const double z = dis(generator) * 2.0 - 1.0;
	const double a = dis(generator) * 2.0 * pi;

	const double r = sqrt(1.0f - z * z);
	const double x = r * cos(a);
	const double y = r * sin(a);

	return glm::normalize(glm::dvec3(x, y, z));
}

double get_intersecting_line_density(
	const long long unsigned int n,
	const double emitter_radius,
	const double receiver_distance,
	const double receiver_distance_plus,
	const double receiver_radius)
{
	double count = 0;
	double count_plus = 0;

	generator.seed(static_cast<unsigned>(0));

	for (long long unsigned int i = 0; i < n; i++)
	{
		if (i % 10000000 == 0)
			cout << double(i) / double(n) << endl;

		glm::dvec3 location = random_unit_vector();

		location.x *= emitter_radius;
		location.y *= emitter_radius;
		location.z *= emitter_radius;

		glm::dvec3 surface_normal = glm::normalize(location);


		// Newtonian gravitation
		// glm::dvec3 normal = surface_normal;


		// Schwarzschild gravitation
		glm::dvec3 normal =
			random_cosine_weighted_hemisphere(
				surface_normal);


		// Schwarzschild gravitation using a useful trick
		// https://pema.dev/obsidian/math/light-transport/cosine-weighted-sampling.html
		//glm::dvec3 normal = glm::normalize(surface_normal + random_unit_vector());


		// Emulate Quantum Graphity
		//glm::dvec3 normal = glm::normalize(random_unit_vector() * emitter_radius - random_unit_vector() * emitter_radius);


		if (dot(normal, surface_normal) < 0)
			surface_normal = -surface_normal;






		count += intersect(
			location, normal,
			receiver_distance, receiver_radius);

		count_plus += intersect(
			location, normal,
			receiver_distance_plus, receiver_radius);
	}

	return count_plus - count;
}

int main(int argc, char** argv)
{
	ofstream outfile("ratio");

	const double emitter_radius_geometrized =
		sqrt(1e8 * log(2.0) / pi);

	const double receiver_radius_geometrized =
		emitter_radius_geometrized * 0.01; // Minimum one Planck unit

	const double emitter_area_geometrized =
		4.0 * pi * pow(emitter_radius_geometrized, 2.0);

	// Field line count
	const double n_geometrized =
		emitter_area_geometrized
		/ (log(2.0) * 4.0);

	const double emitter_mass_geometrized =
		emitter_radius_geometrized
		/ 2.0;

	double start_pos =
		emitter_radius_geometrized
		+ receiver_radius_geometrized;

	double end_pos = start_pos * 10;

	const size_t pos_res = 10; // Minimum 2 steps

	const double pos_step_size =
		(end_pos - start_pos)
		/ (pos_res - 1);

	const double epsilon =
		receiver_radius_geometrized;

	for (size_t i = 0; i < pos_res; i++)
	{
		const double receiver_distance_geometrized =
			start_pos + i * pos_step_size;

		const double receiver_distance_plus_geometrized =
			receiver_distance_geometrized + epsilon;

		const double collision_count_plus_minus_collision_count =
			get_intersecting_line_density(
				static_cast<long long unsigned int>(n_geometrized),
				emitter_radius_geometrized,
				receiver_distance_geometrized,
				receiver_distance_plus_geometrized,
				receiver_radius_geometrized);

		// alpha variable
		const double gradient_integer =
			collision_count_plus_minus_collision_count
			/ epsilon;

		// g variable
		double gradient_strength =
			-gradient_integer /
			(2.0 * pow(receiver_radius_geometrized, 3.0));

		const double a_Newton_geometrized =
			sqrt(n_geometrized * log(2.0) /
				(4.0 * pi * pow(receiver_distance_geometrized, 4.0)));

		const double a_flat_geometrized =
			gradient_strength * receiver_distance_geometrized * log(2)
			/ (8.0 * emitter_mass_geometrized);


		const double dt_Schwarzschild =
			sqrt(1 - emitter_radius_geometrized /
				receiver_distance_geometrized);

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

		outfile << receiver_distance_geometrized <<
			" " <<
			(a_Schwarzschild_geometrized / a_flat_geometrized) <<
			endl;
	}
}




