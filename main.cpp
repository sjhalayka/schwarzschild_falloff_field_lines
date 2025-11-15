#include "main.h"

real_type intersect_AABB(const vector_3 min_location, const vector_3 max_location, const vector_3 ray_origin, const vector_3 ray_dir, real_type& tmin, real_type& tmax)
{
	tmin = (min_location.x - ray_origin.x) / ray_dir.x;
	tmax = (max_location.x - ray_origin.x) / ray_dir.x;

	if (tmin > tmax)
		swap(tmin, tmax);

	real_type tymin = (min_location.y - ray_origin.y) / ray_dir.y;
	real_type tymax = (max_location.y - ray_origin.y) / ray_dir.y;

	if (tymin > tymax)
		swap(tymin, tymax);

	if ((tmin > tymax) || (tymin > tmax))
		return 0;

	if (tymin > tmin)
		tmin = tymin;

	if (tymax < tmax)
		tmax = tymax;

	real_type tzmin = (min_location.z - ray_origin.z) / ray_dir.z;
	real_type tzmax = (max_location.z - ray_origin.z) / ray_dir.z;

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

	vector_3 ray_hit_start = ray_origin;
	ray_hit_start.x += ray_dir.x * tmin;
	ray_hit_start.y += ray_dir.y * tmin;
	ray_hit_start.z += ray_dir.z * tmin;

	vector_3 ray_hit_end = ray_origin;
	ray_hit_end.x += ray_dir.x * tmax;
	ray_hit_end.y += ray_dir.y * tmax;
	ray_hit_end.z += ray_dir.z * tmax;

	real_type l = (ray_hit_end - ray_hit_start).length();

	return l;
}

vector_3 random_unit_vector(void)
{
	const real_type z = dis(generator) * 2.0 - 1.0;
	const real_type a = dis(generator) * 2.0 * pi;

	const real_type r = sqrt(1.0f - z * z);
	const real_type x = r * cos(a);
	const real_type y = r * sin(a);

	return vector_3(x, y, z).normalize();
}

std::optional<real_type> intersect(
	const vector_3 location,
	const vector_3 normal,
	const real_type receiver_distance,
	const real_type receiver_radius)
{
	const vector_3 circle_origin(receiver_distance, 0, 0);

	if (normal.dot(circle_origin) <= 0)
		return std::nullopt;

	vector_3 min_location(-receiver_radius + receiver_distance, -receiver_radius, -receiver_radius);
	vector_3 max_location(receiver_radius + receiver_distance, receiver_radius, receiver_radius);

	real_type tmin = 0, tmax = 0;

	real_type AABB_hit = intersect_AABB(min_location, max_location, location, normal, tmin, tmax);

	if (AABB_hit > 0)
		return AABB_hit;

	return std::nullopt;
}

vector_3 random_tangent_vector(const vector_3& point_on_sphere)
{
	// Normalize to ensure it's on unit sphere
	vector_3 normal = point_on_sphere;
	normal.normalize();

	// Choose an arbitrary vector that's not parallel to normal
	vector_3 arbitrary;
	if (fabs(normal.x) > 0.9)
		arbitrary = vector_3(0, 1, 0);  // If normal is mostly along x, use y
	else
		arbitrary = vector_3(1, 0, 0);  // Otherwise use x

	// Get first basis vector perpendicular to normal
	vector_3 tangent1 = normal.cross(arbitrary);
	tangent1.normalize();

	// Get second basis vector perpendicular to both
	vector_3 tangent2 = normal.cross(tangent1);
	tangent2.normalize();

	// Generate random angle for rotation in tangent plane
	real_type angle = dis(generator) * 2.0 * pi;

	// Combine the two tangent vectors with random weights
	vector_3 result = tangent1 * cos(angle) + tangent2 * sin(angle);

	return result.normalize();
}


vector_3 random_cosine_weighted_hemisphere(const vector_3& normal)
{
	// Generate two random numbers
	real_type u1 = dis(generator);
	real_type u2 = dis(generator);

	// Malley's method
	// (cosine-weighted hemisphere sampling)
	// Sample uniformly on a disk, 
	// then project up to hemisphere
	real_type r = sqrt(u1);
	real_type theta = 2.0 * pi * u2;

	// Point on unit disk
	real_type x = r * cos(theta);
	real_type y = r * sin(theta);
	real_type z = sqrt(1.0 - u1); // Height above disk

	// Create orthonormal basis around normal
	vector_3 n = normal;
	n.normalize();

	// Choose an arbitrary vector not parallel to normal
	vector_3 arbitrary;
	if (fabs(n.x) > 0.9)
		arbitrary = vector_3(0, 1, 0);
	else
		arbitrary = vector_3(1, 0, 0);

	// Create tangent and bitangent
	vector_3 tangent = n.cross(arbitrary);
	tangent.normalize();

	vector_3 bitangent = n.cross(tangent);
	bitangent.normalize();

	// Transform from local coordinates
	// to world coordinates
	vector_3 result;
	result.x = tangent.x * x + 
		bitangent.x * y + 
		n.x * z;

	result.y = tangent.y * x + 
		bitangent.y * y + 
		n.y * z;

	result.z = tangent.z * x + 
		bitangent.z * y + 
		n.z * z;

	return result.normalize();
}


real_type get_intersecting_line_density(
	const long long unsigned int n,
	const real_type emitter_radius,
	const real_type receiver_distance,
	const real_type receiver_distance_plus,
	const real_type receiver_radius)
{
	real_type count = 0;
	real_type count_plus = 0;

	generator.seed(static_cast<unsigned>(0));

	for (long long unsigned int i = 0; i < n; i++)
	{
		if (i % 100000000 == 0)
			cout << float(i) / float(n) << endl;

		// Random hemisphere outward
		vector_3 location = random_unit_vector();

		location.x *= emitter_radius;
		location.y *= emitter_radius;
		location.z *= emitter_radius;

		vector_3 surface_normal = location;
		surface_normal.normalize();

		vector_3 normal = 
			random_cosine_weighted_hemisphere(
				surface_normal);

		std::optional<real_type> i_hit = intersect(
			location, normal, 
			receiver_distance, receiver_radius);

		if (i_hit)
			count += *i_hit / (2.0 * receiver_radius);
	
		i_hit = intersect(
			location, normal,
			receiver_distance_plus, receiver_radius);

		if (i_hit)
			count_plus += *i_hit / (2.0 * receiver_radius);
	}

	return count_plus - count;
}

int main(int argc, char** argv)
{
	ofstream outfile("ratio");

	const real_type emitter_radius_geometrized =
		sqrt(1e8 * log(2.0) / pi);

	const real_type receiver_radius_geometrized =
		emitter_radius_geometrized * 0.01; // Minimum one Planck unit

	const real_type emitter_area_geometrized =
		4.0 * pi
		* emitter_radius_geometrized
		* emitter_radius_geometrized;

	// Field line count
	const real_type n_geometrized =
		emitter_area_geometrized
		/ (log(2.0) * 4.0);

	const real_type emitter_mass_geometrized =
		emitter_radius_geometrized
		/ 2.0;

	real_type start_pos =
		emitter_radius_geometrized
		+ receiver_radius_geometrized;

	real_type end_pos = start_pos * 10;

	//swap(end_pos, start_pos);

	const size_t pos_res = 10; // Minimum 2 steps

	const real_type pos_step_size =
		(end_pos - start_pos)
		/ (pos_res - 1);

	const real_type epsilon =
		receiver_radius_geometrized;


	for (size_t i = 0; i < pos_res; i++)
	{
		const real_type receiver_distance_geometrized =
			start_pos + i * pos_step_size;

		const real_type receiver_distance_plus_geometrized =
			receiver_distance_geometrized + epsilon;

		// beta function
		const real_type collision_count_plus_minus_collision_count =
			get_intersecting_line_density(
				static_cast<long long unsigned int>(n_geometrized),
				emitter_radius_geometrized,
				receiver_distance_geometrized,
				receiver_distance_plus_geometrized,
				receiver_radius_geometrized);

		// alpha variable
		const real_type gradient_integer =
			collision_count_plus_minus_collision_count
			/ epsilon;

		// g variable
		real_type gradient_strength =
			-gradient_integer
			/
			(receiver_radius_geometrized
				* receiver_radius_geometrized
				);

		//cout << gradient_strength << " " << n_geometrized / (2 * pow(receiver_distance_geometrized, 3.0)) << endl;
		//cout << gradient_strength / (n_geometrized / (2 * pow(receiver_distance_geometrized, 3.0))) << endl;


		const real_type a_Newton_geometrized =
			sqrt(
				n_geometrized * log(2.0)
				/
				(4.0 * pi *
					pow(receiver_distance_geometrized, 4.0))
			);

		const real_type a_flat_geometrized =
			gradient_strength * receiver_distance_geometrized * log(2)
			/ (8.0 * emitter_mass_geometrized);


		//const real_type g_approx = n_geometrized / (2 * pow(receiver_distance_geometrized, 3.0));
		//const real_type a_approx_geometrized =
		//	g_approx * receiver_distance_geometrized * log(2)
		//	/ (8.0 * emitter_mass_geometrized);


		const real_type dt_Schwarzschild = sqrt(1 - emitter_radius_geometrized / receiver_distance_geometrized);

		const real_type a_Schwarzschild_geometrized =
			emitter_radius_geometrized / (pi * pow(receiver_distance_geometrized, 2.0) * dt_Schwarzschild);

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




