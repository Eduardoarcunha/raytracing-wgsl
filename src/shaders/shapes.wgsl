// Tests if a ray intersects with a sphere
// Parameters: center and radius of sphere, the ray, hit record to store results, and max distance
fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
    // Calculate vector from ray origin to sphere center
    let oc = r.origin - center;
    
    // Calculate quadratic equation coefficients
    // For a sphere: (P(t) - C)² = R², where P(t) = A + tb
    // Expanding gives us: (A + tb - C)² = R²
    // This becomes: t²(b·b) + 2tb·(A-C) + ((A-C)·(A-C)) - R² = 0
    // Which is a quadratic equation: at² + bt + c = 0
    let a = dot(r.direction, r.direction);    // b·b
    let half_b = dot(oc, r.direction);        // b·(A-C)
    let c = dot(oc, oc) - radius * radius;    // (A-C)·(A-C) - R²
    
    // Calculate discriminant to determine number of intersections
    // discriminant = b² - 4ac
    // We use half_b instead of b, so formula becomes: discriminant = 4(half_b² - ac)
    let discriminant = half_b * half_b - a * c;

    // If discriminant < 0, ray misses sphere (no real solutions)
    if (discriminant < 0.0) {
        record.hit_anything = false;
        return;
    }

    // Calculate the nearest intersection distance
    // t = (-b ± √discriminant) / 2a
    // Since we used half_b = b/2, formula becomes: t = (-half_b ± √discriminant) / a
    let sqrtd = sqrt(discriminant);
    
    // Try smaller t first (nearest intersection)
    var root = (-half_b - sqrtd) / a;
    
    // If nearest intersection is outside valid range, try farther intersection
    if (root < RAY_TMIN || root > max) {
        root = (-half_b + sqrtd) / a;
        if (root < RAY_TMIN || root > max) {
            record.hit_anything = false;
            return;
        }
    }

    record.t = root;
    record.p = ray_at(r, root);
    let outward_normal = (record.p - center) / radius;
    
    // Set frontface and adjust normal direction
    record.frontface = dot(r.direction, outward_normal) < 0.0;
    record.normal = select(-outward_normal, outward_normal, record.frontface);
    
    record.hit_anything = true;
}

// Tests if a ray intersects with a quadrilateral (a flat 4-sided shape)
// Parameters: ray, Q (corner point), u and v (vectors defining quad edges), hit record, and max distance
fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
    // Calculate the normal of the quad using cross product of its edges
    var n = cross(u.xyz, v.xyz);
    var normal = normalize(n);
    // D is the distance from origin to the plane containing the quad
    var D = dot(normal, Q.xyz);
    // w is used for calculating barycentric coordinates
    var w = n / dot(n.xyz, n.xyz);

    // Check if ray is parallel to quad's plane (would cause division by zero)
    var denom = dot(normal, r.direction);
    if (abs(denom) < 0.0001)
    {
        record.hit_anything = false;
        return;
    }

    // Calculate distance along ray to plane intersection
    var t = (D - dot(normal, r.origin)) / denom;
    // Check if intersection is within valid distance range
    if (t < RAY_TMIN || t > max)
    {
        record.hit_anything = false;
        return;
    }

    // Calculate the actual intersection point
    var intersection = ray_at(r, t);
    // Vector from quad corner to intersection point
    var planar_hitpt_vector = intersection - Q.xyz;
    // Calculate barycentric coordinates (alpha, beta)
    // These tell us if the point lies within the quad
    var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
    var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

    // Check if point lies outside the quad's boundaries
    if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
    {
        record.hit_anything = false;
        return;
    }

    // Check if we're hitting the back face of the quad
    if (dot(normal, r.direction) > 0.0)
    {
        record.hit_anything = false;
        return;
    }

    // Record successful hit information
    record.t = t;
    record.p = intersection;
    record.normal = normal;
    record.hit_anything = true;
}

// Tests if a ray intersects with a triangle
// Parameters: ray, three vertices (v0, v1, v2), hit record, and max distance
// Uses Möller–Trumbore algorithm
fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
    // Calculate two edges of the triangle
    var v1v0 = v1 - v0;
    var v2v0 = v2 - v0;
    var rov0 = r.origin - v0;

    // Calculate triangle normal and helper vector
    var n = cross(v1v0, v2v0);
    var q = cross(rov0, r.direction);

    // Calculate denominator for barycentric coordinates
    var d = 1.0 / dot(r.direction, n);

    // Calculate barycentric coordinates (u,v)
    var u = d * dot(-q, v2v0);
    var v = d * dot(q, v1v0);
    var t = d * dot(-n, rov0);

    // Check if point lies outside triangle (using barycentric coordinates)
    if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
    {
        record.hit_anything = false;
        return;
    }

    // Check if intersection is within valid distance range
    if (t < RAY_TMIN || t > max)
    {
        record.hit_anything = false;
        return;
    }

    // Record successful hit information
    record.t = t;
    record.p = ray_at(r, t);
    record.normal = normalize(n);
    record.hit_anything = true;
}

// Tests if a ray intersects with an axis-aligned box
// Parameters: ray, center of box, radius (half-dimensions), hit record, and max distance
// Uses "slab" method for efficient intersection
fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
    // Calculate inverse of ray direction and related values
    var m = 1.0 / r.direction;
    var n = m * (r.origin - center);
    var k = abs(m) * rad;

    // Calculate intersection distances with all slabs
    var t1 = -n - k;  // Entry points
    var t2 = -n + k;  // Exit points

    // Find the latest entry and earliest exit points
    var tN = max(max(t1.x, t1.y), t1.z);  // Latest entry
    var tF = min(min(t2.x, t2.y), t2.z);  // Earliest exit

    // Check if ray misses box
    if (tN > tF || tF < 0.0)
    {
        record.hit_anything = false;
        return;
    }

    // Use entry point for intersection
    var t = tN;
    // Check if intersection is within valid distance range
    if (t < RAY_TMIN || t > t_max)
    {
        record.hit_anything = false;
        return;
    }

    // Record successful hit information
    record.t = t;
    record.p = ray_at(r, t);
    // Calculate normal based on which face was hit
    record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
    record.hit_anything = true;

    return;
}