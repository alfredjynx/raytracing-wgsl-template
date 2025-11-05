fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{

  var d = r.direction;
  
  var A = r.origin;

  var e = A - center;

  var a  = dot(d, d);
  
  var b = 2.0 * dot(d, e);

  var c = dot(e, e) - radius * radius;

  var delta = b * b - 4.0 * a * c;

  if (delta < 0.0)  {
    record.hit_anything = false;
    return;
  }

  var sqrt_delta = sqrt(delta);

  var t1 = (-b - sqrt_delta) / (2.0 * a);
  var t2 = (-b + sqrt_delta) / (2.0 * a);

  var t = t1;

  if (t2 < t1) {
    t = t2;
  } 

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;  
  record.p = ray_at(r, t);
  record.normal = (record.p - center) / radius;

  record.frontface = dot(record.normal, r.direction) < 0.0;
  record.hit_anything = true;
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}


// https://hugi.scene.org/online/hugi24/coding%20graphics%20chris%20dragan%20raytracing%20shapes.htm
fn hit_cylinder(r: ray,
                C: vec3f,        // base-center
                V: vec3f,        // axis (unit length)
                radius: f32,     // cylinder radius
                height: f32,     // finite: end cap at C + V * height
                record: ptr<function, hit_record>,
                t_max: f32)
{
  const EPS: f32 = 1e-6;

  // Ray and relative vectors
  let O = r.origin;
  let D = r.direction;
  let X = O - C;

  // Axis projections
  let Dv = dot(D, V);
  let Xv = dot(X, V);

  // ---- Side (exactly as in the article) ----
  // a   = D·D - (D·V)^2
  // b/2 = D·X - (D·V)*(X·V)
  // c   = X·X - (X·V)^2 - r^2
  let a = dot(D, D) - Dv * Dv;
  let half_b = dot(D, X) - Dv * Xv;
  let c = dot(X, X) - Xv * Xv - radius * radius;

  var have_side = false;
  var t_side    = 1e30;
  var p_side    = vec3f(0.0);
  var n_side    = vec3f(0.0);

  if (abs(a) > EPS) {
    let disc = half_b * half_b - a * c;
    if (disc >= 0.0) {
      let s  = sqrt(disc);
      let t0 = (-half_b - s) / a;
      let t1 = (-half_b + s) / a;

      // try nearer then farther, both inline, no helpers
      if (t0 >= RAY_TMIN && t0 <= t_max) {
        let m0 = Dv * t0 + Xv;               // axial coordinate at t0
        if (m0 >= 0.0 && m0 <= height) {
          have_side = true;
          t_side = t0;
          p_side = O + D * t0;
          n_side = normalize(p_side - (C + V * m0)); // side normal
        }
      }
      if (!have_side && t1 >= RAY_TMIN && t1 <= t_max) {
        let m1 = Dv * t1 + Xv;
        if (m1 >= 0.0 && m1 <= height) {
          have_side = true;
          t_side = t1;
          p_side = O + D * t1;
          n_side = normalize(p_side - (C + V * m1));
        }
      }
    }
  }

  // ---- Caps (planes + disk test, per article) ----
  // Start cap: plane (C, -V)
  var have_cap = false;
  var t_cap    = 1e30;
  var p_cap    = vec3f(0.0);
  var n_cap    = vec3f(0.0);

  let n0 = -V;
  let denom0 = dot(n0, D);
  if (abs(denom0) > EPS) {
    let t0p = dot(n0, (C - O)) / denom0;
    if (t0p >= RAY_TMIN && t0p <= t_max) {
      let P0 = O + D * t0p;
      // disk radius check
      if (length(P0 - C) <= radius + 1e-6) {
        have_cap = true;
        t_cap = t0p;
        p_cap = P0;
        n_cap = n0; // -V
      }
    }
  }

  // End cap: plane (C + V*height, +V)
  let Ce = C + V * height;
  let n1 = V;
  let denom1 = dot(n1, D);
  if (abs(denom1) > EPS) {
    let t1p = dot(n1, (Ce - O)) / denom1;
    if (t1p >= RAY_TMIN && t1p <= t_max) {
      let P1 = O + D * t1p;
      if (length(P1 - Ce) <= radius + 1e-6) {
        if (!have_cap || t1p < t_cap) {
          have_cap = true;
          t_cap = t1p;
          p_cap = P1;
          n_cap = n1; // +V
        }
      }
    }
  }

  // ---- Choose nearest valid intersection (side vs caps) ----
  var have_hit = false;
  var t_final  = 0.0;
  var p_final  = vec3f(0.0);
  var n_final  = vec3f(0.0);

  if (have_side && have_cap) {
    if (t_side < t_cap) {
      have_hit = true; t_final = t_side; p_final = p_side; n_final = n_side;
    } else {
      have_hit = true; t_final = t_cap;  p_final = p_cap;  n_final = n_cap;
    }
  } else if (have_side) {
    have_hit = true; t_final = t_side; p_final = p_side; n_final = n_side;
  } else if (have_cap) {
    have_hit = true; t_final = t_cap;  p_final = p_cap;  n_final = n_cap;
  }

  if (!have_hit) {
    record.hit_anything = false;
    return;
  }

  // Output (same style as your other hit_* functions)
  record.t = t_final;
  record.p = p_final;
  record.normal = n_final;
  record.hit_anything = true;
}
