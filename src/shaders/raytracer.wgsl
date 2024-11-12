const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record {
    // Get counts of different object types from uniforms
    var spheresCount = i32(uniforms[19]);
    var quadsCount = i32(uniforms[20]);
    var boxesCount = i32(uniforms[21]);
    var trianglesCount = i32(uniforms[22]);
    var meshCount = i32(uniforms[27]);

    // Initialize records for hit testing
    var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var closest = record;
    var closest_so_far = max;

    // Check all spheres
    for(var i = 0; i < spheresCount; i++) {
        var sphere = spheresb[i];
        var local_record = record;
        
        hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &local_record, closest_so_far);
        
        if(local_record.hit_anything) {
            closest_so_far = local_record.t;
            local_record.object_color = sphere.color;
            local_record.object_material = sphere.material;
            closest = local_record;
        }
    }

    // Check all quads
    for(var i = 0; i < quadsCount; i++) {
        var quad = quadsb[i];
        var local_record = record;
        
        hit_quad(r, quad.Q, quad.u, quad.v, &local_record, closest_so_far);
        
        if(local_record.hit_anything) {
            closest_so_far = local_record.t;
            local_record.object_color = quad.color;
            local_record.object_material = quad.material;
            closest = local_record;
        }
    }

    // Check all boxes
    for(var i = 0; i < boxesCount; i++) {
        var box = boxesb[i];
        var local_record = record;
        
        hit_box(r, box.center.xyz, box.radius.xyz, &local_record, closest_so_far);
        
        if(local_record.hit_anything) {
            closest_so_far = local_record.t;
            local_record.object_color = box.color;
            local_record.object_material = box.material;
            closest = local_record;
        }
    }

    // Check all triangles
    for(var i = 0; i < trianglesCount; i++) {
        var tri = trianglesb[i];
        var local_record = record;
        
        hit_triangle(r, tri.v0.xyz, tri.v1.xyz, tri.v2.xyz, &local_record, closest_so_far);
        
        if(local_record.hit_anything) {
            closest_so_far = local_record.t;
            // Note: Triangles might need color/material from a different source
            // like mesh or additional properties in triangle struct
            closest = local_record;
        }
    }

    // Check all meshes (collections of triangles)
    for(var i = 0; i < meshCount; i++) {
        var mesh = meshb[i];
        
        // Skip if mesh is hidden
        if(mesh.show_bb == 0.0) {
            continue;
        }
        
        // Check each triangle in the mesh
        for(var j = i32(mesh.start); j < i32(mesh.end); j++) {
            var tri = trianglesb[j];
            var local_record = record;
            
            hit_triangle(r, tri.v0.xyz, tri.v1.xyz, tri.v2.xyz, &local_record, closest_so_far);
            
            if(local_record.hit_anything) {
                closest_so_far = local_record.t;
                local_record.object_color = mesh.color;
                local_record.object_material = mesh.material;
                closest = local_record;
            }
        }
    }

    return closest;
}
// Simulates a diffuse (Lambertian) material reflection
// Parameters:
//   normal: Surface normal at hit point
//   absorption: Material roughness (affects scatter direction)
//   random_sphere: Random point in unit sphere for scatter direction
//   rng_state: Random number generator state for additional randomization
fn lambertian(normal: vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour {
    var scatter_direction: vec3f;
    
    // Perfect lambertian reflection would use normal + random unit vector
    // Modified by absorption (roughness) to allow for more control
    if (absorption <= 0.0) {
        // Pure diffuse reflection: scatter in random direction in hemisphere
        scatter_direction = normal + random_sphere;
    } else {
        // Mix between pure diffuse and more focused reflection
        // Higher absorption means more scattering closer to the normal
        scatter_direction = normalize(normal + absorption * random_sphere);
    }

    // Catch degenerate scatter direction
    // If scatter_direction is zero (or very close), use the normal
    if (length(scatter_direction) < 0.001) {
        scatter_direction = normal;
    }

    // Always scatter (true) with the calculated direction
    return material_behaviour(true, normalize(scatter_direction));
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflect = direction - 2 * dot(direction, normal) * normal;
  var scatter_direction = reflect + fuzz * random_sphere;

  return material_behaviour(true, normalize(reflect));
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{  
  return material_behaviour(false, vec3f(0.0));
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
  return material_behaviour(false, vec3f(0.0));
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f {
    var maxbounces = i32(uniforms[2]);
    var light = vec3f(0.0);     // Accumulated light/color
    var color = vec3f(1.0);     // Current ray color multiplier
    var r_ = r;                 // Current ray being traced
    
    // Get background colors from uniforms
    var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
    var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
    var behaviour = material_behaviour(true, vec3f(0.0));

    // Main ray bouncing loop
    for (var j = 0; j < maxbounces; j = j + 1) {
        // Check for intersection with any object in the scene
        var hit = check_ray_collision(r_, RAY_TMAX);
        
        // If we didn't hit anything, add background contribution and break
        if (!hit.hit_anything) {
            light += color * envoriment_color(r_.direction, backgroundcolor1, backgroundcolor2);
            break;
        }

        // Get material properties from hit record
        var mat_type = i32(hit.object_material.x);
        var roughness = hit.object_material.y;
        var ior = hit.object_material.z;     // Index of refraction
        var emission = hit.object_material.w; // Emission strength
        
        // Get a random point in unit sphere for material calculations
        var random_in_sphere = rng_next_vec3_in_unit_sphere(rng_state);

        // Handle different material types
        switch(mat_type) {
            // Diffuse
            case 0: {
                behaviour = lambertian(hit.normal, roughness, random_in_sphere, rng_state);
                color *= hit.object_color.rgb;
            }
            // Metallic
            case 1: {
                behaviour = metal(hit.normal, r_.direction, roughness, random_in_sphere);
                color *= hit.object_color.rgb;
            }
            // Glass/Dielectric
            case 2: {
                behaviour = dielectric(hit.normal, r_.direction, ior, hit.frontface, 
                                    random_in_sphere, roughness, rng_state);
                color *= vec3f(0.98);
            }
            // Emissive
            case 3: {
                behaviour = emmisive(hit.object_color.rgb, emission);
                light += color * hit.object_color.rgb * emission;
                break;
            }
            default: {
                behaviour = lambertian(hit.normal, roughness, random_in_sphere, rng_state);
                color *= hit.object_color.rgb;
            }
        }

        // If the material absorbed the ray (didn't scatter), break
        if (!behaviour.scatter) {
            break;
        }

        // Update ray for next bounce
        r_ = ray(hit.p, behaviour.direction);
    }

    return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera setup
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    
    // Get number of samples per pixel from uniforms
    var samples_per_pixel = i32(uniforms[4]);
    
    // Initialize accumulated color
    var color = vec3f(0.0);
    
    // 1. Loop for each sample per pixel
    for(var s = 0; s < samples_per_pixel; s++) {
        // 2. Get ray for this sample
        var r = get_ray(cam, uv, &rng_state);
        
        // 3. Call trace function to get color for this ray
        var sample_color = trace(r, &rng_state);
        
        // Accumulate the color
        color += sample_color;
    }
    
    // 4. Average the color by dividing by number of samples
    color = color / f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    
    var map_fb = mapfb(id.xy, rez);
    
    var should_accumulate = uniforms[3];
    var accumulated_color = rtfb[map_fb] * should_accumulate + color_out;

    rtfb[map_fb] = accumulated_color;
    fb[map_fb] = accumulated_color/accumulated_color.w;
}