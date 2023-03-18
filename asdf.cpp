#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;
  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading
    int num_intersects = 0;
    for (int i = 0; i < num_samples; i++) {
            Vector3D wj = hemisphereSampler->get_sample(); //point in obj frame
            Vector3D wj_world = (o2w * wj);
            double pdf = 1.0 / (2.0 * PI);
            Ray sample_ray(hit_p, wj_world); //rays should use world frame
            sample_ray.min_t = EPS_F;
            Intersection sample_isect;

            if (bvh->intersect(sample_ray, &sample_isect)) { //bsdf should use object frame
                Vector3D emission = sample_isect.bsdf->get_emission();
                Vector3D reflectance = isect.bsdf->f(w_out, wj);
                Vector3D L = (reflectance * emission * cos_theta(wj)) / pdf;
                L_out += L;
                num_intersects++;
            }
        }
  return L_out / (double) num_samples;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;
    int num_lights = scene->lights.size();
    for (int i = 0; i < num_lights; ++i) {

        Vector3D wi;
        double disToLight;
        double pdf;
        int num_samples = 1;

        if (!scene->lights[i]->is_delta_light()) {
            num_samples = ns_area_light;
        }

        Vector3D samp_L(0.0, 0.0, 0.0);

        for (int j = 0; j < num_samples; ++j) {
            Vector3D emitted_radiance = scene->lights[i]->sample_L(hit_p, &wi, &disToLight, &pdf);
            Vector3D wi_obj = (w2o * wi);

            Ray shadow_ray(hit_p, wi);
            shadow_ray.min_t = EPS_F;
            shadow_ray.max_t = disToLight - EPS_F;

            Intersection shadow_isect;
            if ((cos_theta(wi_obj) > 0.0) && (pdf != 0.0) && !bvh->intersect(shadow_ray, &shadow_isect)) {
                Vector3D reflectance = isect.bsdf->f(w_out, wi_obj);
                Vector3D L = (reflectance * emitted_radiance * cos_theta(wi_obj)) / pdf;
                samp_L += L;
            }
        }
        samp_L /= (double) num_samples;
        L_out += samp_L;
    }
    return L_out;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light
    return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
    if (direct_hemisphere_sample) {
        return estimate_direct_lighting_hemisphere(r, isect);
    } else {
        return estimate_direct_lighting_importance(r, isect);
    }
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.
    L_out += one_bounce_radiance(r, isect);

    if (r.depth <= 1) {
        return L_out;
    }

    double p = 0.65;
    if ((r.depth == max_ray_depth) || coin_flip(p)) {
        Vector3D wi;
        double pdf;
        Vector3D radiance = isect.bsdf->sample_f(w_out, &wi, &pdf);
        Vector3D wi_world = (o2w * wi);

        Ray bounce_ray(hit_p, wi_world);
        bounce_ray.min_t = EPS_F;
        bounce_ray.depth = r.depth - 1;
        Intersection bounce_isect;

        if ((pdf != 0.0) && bvh->intersect(bounce_ray, &bounce_isect)) {
            Vector3D bounced_light =  at_least_one_bounce_radiance(bounce_ray, bounce_isect);
            Vector3D L = ((bounced_light * radiance * cos_theta(wi)) / pdf);
            if (r.depth != max_ray_depth) {
                L /= p;
            }
            L_out += L;
        }

    }

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.

  // TODO (Part 3): Return the direct illumination.
    if (!bvh->intersect(r, &isect)) {
        return L_out;
    }

//    return normal_shading(isect.n); // Normal Shading

//    return zero_bounce_radiance(r, isect) + one_bounce_radiance(r, isect); // Direct Lighting
    return at_least_one_bounce_radiance(r, isect) - one_bounce_radiance(r, isect); // Indirect Lighting
//    return zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect); // Direct + Indirect Lighting

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;          // total samples to evaluate
    Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

    Vector3D pixel_color(0.0, 0.0, 0.0);
    double illum = 0.0;
    double illum2 = 0.0;
    int sample_cnt = 0;
    bool converged = false;
    for (int n = 0; !converged && (n < num_samples); n++) {
        Vector2D randSample = gridSampler->get_sample();
        double sampX = (origin.x + randSample.x) / (double)sampleBuffer.w;
        double sampY = (origin.y + randSample.y) / (double)sampleBuffer.h;
        Ray r = camera->generate_ray(sampX, sampY);
        r.depth = max_ray_depth;
        Vector3D L_out = est_radiance_global_illumination(r);
        pixel_color += L_out;
        sample_cnt++;
        illum += L_out.illum();
        illum2 += std::pow(L_out.illum(), 2);

        if ((sample_cnt % samplesPerBatch) == 0) {
            double mean = illum / sample_cnt;
            double std_dv2 = (1.0 / (sample_cnt - 1.0)) * (illum2 - (std::pow(illum, 2) / sample_cnt));

            double I = 1.96 * std::sqrt(std_dv2 / sample_cnt);
            if (I <= (maxTolerance * mean)) {
                converged = true;
            }
        }
    }

    pixel_color /= sample_cnt;

  sampleBuffer.update_pixel(pixel_color, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = sample_cnt;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
