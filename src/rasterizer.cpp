#include "rasterizer.h"
#include <CGL/vector3D.h>
using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    //加if循环
      double rate1d = sqrt(sample_rate);
      for (int i = 0; i < rate1d; i++) {  //limit the ditection area into boundery square of triangle
          for (int j = 0; j < rate1d; j++) {
              int buffer_index = ((y * rate1d + j) * width * rate1d) + (x * rate1d + i);
              sample_buffer[buffer_index] = c;
          }
      }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // helper function to determine point in triangle
  bool point_in_triangle(float x, float y,
      Vector2D A, Vector2D B, Vector2D C) {
      float valueAB = -(x - A.x) * (B.y - A.y) + (y - A.y) * (B.x - A.x);
      float valueBC = -(x - B.x) * (C.y - B.y) + (y - B.y) * (C.x - B.x);
      float valueCA = -(x - C.x) * (A.y - C.y) + (y - C.y) * (A.x - C.x);

      if ((valueAB > 0 && valueBC >= 0 && valueCA >= 0) || (valueAB <= 0 && valueBC <= 0 && valueCA <= 0)) {
          return true;
      }
      else {
          return false;
      }
      //return true;
  }



  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    float minX = min({ x0, min(x1, x2)});
    float maxX = max({ x0, max(x1, x2)});
    float minY = min({ y0, min(y1, y2)});
    float maxY = max({ y0, max(y1, y2)});

    Vector2D A(x0, y0);
    Vector2D B(x1, y1);
    Vector2D C(x2, y2);

    // TODO: Task 2: Update to implement super-sampled rasterization
    double rate1d = sqrt(sample_rate);

    for (int x = floor(minX); x < ceil(maxX); x++) {  //limit the ditection area into boundery square of triangle
        for (int y = floor(minY); y < ceil(maxY); y++) {
              //Image[x][y] = point_in_triangle(x + 0.5, y + 0.5);
              // calculate distance between side
              // super-sampled rasterization

            for (int i = 0; i < rate1d; i++) {
                for (int j = 0; j < rate1d; j++) {
                    float subpixel_x = (float)x+(float)i / rate1d;
                    float subpixel_y = (float)y+(float)j / rate1d;

                    if (point_in_triangle(subpixel_x, subpixel_y, A, B, C) == true) {
                        int buffer_index = ((y * rate1d + j) * width * rate1d) + (x * rate1d + i);
                        sample_buffer[buffer_index] = color;
                    }

                }
            }
        }
    }
  }

  //helper funxtion
  Vector3D computeBarycentric2D(float x, float y,
      Vector2D A, Vector2D B, Vector2D C) {
      float x0 = A.x, y0 = A.y;
      float x1 = B.x, y1 = B.y;
      float x2 = C.x, y2 = C.y;

      float denominator = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
      float lambda0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / denominator;
      float lambda1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / denominator;
      float lambda2 = 1.0f - lambda0 - lambda1;

      return Vector3D(lambda0, lambda1, lambda2);
  }

  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    // Step 1: Calculate bounding box
      float minX = min({ x0, min(x1, x2) });
      float maxX = max({ x0, max(x1, x2) });
      float minY = min({ y0, min(y1, y2) });
      float maxY = max({ y0, max(y1, y2) });
      Vector2D A(x0, y0);
      Vector2D B(x1, y1);
      Vector2D C(x2, y2);
      double rate1d = sqrt(sample_rate);
      // Step 2: Iterate over each pixel in the bounding box
      for (int x = floor(minX); x <= ceil(maxX); x++) {
          for (int y = floor(minY); y <= ceil(maxY); y++) {
              for (int i = 0; i < rate1d; i++) {
                  for (int j = 0; j < rate1d; j++) {
                      float subpixel_x = (float)x + (float)i / rate1d;
                      float subpixel_y = (float)y + (float)j / rate1d;
                      // Step 3: Compute barycentric coordinates
                      Vector3D barycentricCoords = computeBarycentric2D(x, y, Vector2D(x0, y0), Vector2D(x1, y1), Vector2D(x2, y2));

                      // Check if the pixel is inside the triangle
                      if (barycentricCoords.x >= 0 && barycentricCoords.y >= 0 && barycentricCoords.z >= 0) {

                          // Step 4: Use barycentric coordinates to interpolate the color
                          Color color = barycentricCoords.x * c0 + barycentricCoords.y * c1 + barycentricCoords.z * c2;

                          // Step 5: Set the pixel's color
                          //fill_pixel(x, y, color);
                          int buffer_index = ((y * rate1d + j) * width * rate1d) + (x * rate1d + i);
                          sample_buffer[buffer_index] = color;
                      }
                  }
              }
          }
      }




  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
      float minX = min({ x0, min(x1, x2) });
      float maxX = max({ x0, max(x1, x2) });
      float minY = min({ y0, min(y1, y2) });
      float maxY = max({ y0, max(y1, y2) });
      Vector2D A(x0, y0);
      Vector2D B(x1, y1);
      Vector2D C(x2, y2);
      double rate1d = sqrt(sample_rate);
      // Iterate over each pixel in the bounding box
      for (float x = floor(minX); x <= ceil(maxX); ++x) {
          for (float y = floor(minY); y <= ceil(maxY);++y) {
              for (int i = 0; i < rate1d; i++) {
                  for (int j = 0; j < rate1d; j++) {
                      float subpixel_x = (float)x + (float)i / rate1d;
                      float subpixel_y = (float)y + (float)j / rate1d;
                      // SampleParams structure
                      SampleParams sp;
                      // Check if the point is inside the triangle
                      if (point_in_triangle(subpixel_x, subpixel_y, A,B,C)) {

                          // Calculate barycentric coordinates for the current point
                          Vector3D bary = computeBarycentric2D(x, y, A,B,C);
                          Vector3D bary_dx = computeBarycentric2D(x + 1, y, A, B, C);
                          Vector3D bary_dy = computeBarycentric2D(x, y + 1, A, B, C);

                          // Interpolate UV coordinates using barycentric coordinates
                          sp.p_uv = Vector2D(bary.x * u0 + bary.y * u1 + bary.z * u2,
                              bary.x * v0 + bary.y * v1 + bary.z * v2);

                          // Calculat====e UV coordinates for neighboring pixels to get differentials
                          sp.p_dx_uv = Vector2D(bary_dx.x * u0 + bary_dx.y * u1 + bary_dx.z * u2,
                              bary_dx.x * v0 + bary_dx.y * v1 + bary_dx.z * v2);
                          sp.p_dy_uv = Vector2D(bary_dy.x * u0 + bary_dy.y * u1 + bary_dy.z * u2,
                              bary_dy.x * v0 + bary_dy.y * v1 + bary_dy.z * v2);
                          sp.psm = psm;
                          sp.lsm = lsm;
 

                          // Sample the texture color using the current texture sampling method
                          Color color = tex.sample(sp);

                          // Set the pixel color
                          int buffer_index = ((y * rate1d + j) * width * rate1d) + (x * rate1d + i);
                          sample_buffer[buffer_index] = color;
                      }
                  }
              }
          }
      }


  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
      int rate1d = sqrt(sample_rate);
      int samples_per_pixel = rate1d * rate1d;

      for (int x = 0; x < width; ++x) {
          for (int y = 0; y < height; ++y) {
              // Initialize the cumulative color value to 0
              Color pixel_color(0, 0, 0);

              // Iterate over all supersampled sub-pixels corresponding to each pixel
              for (int sx = 0; sx < rate1d; ++sx) {
                  for (int sy = 0; sy < rate1d; ++sy) {
                      // Calculate the index of a subpixel in the supersampled buffer
                      int sample_index = ((y * rate1d + sy) * width * rate1d) + (x * rate1d + sx);
                      // Accumulate the color values ​​of all subpixels
                      pixel_color += sample_buffer[sample_index];
                  }
              }

              // Calculate average color value
              pixel_color.r /= samples_per_pixel;
              pixel_color.g /= samples_per_pixel;
              pixel_color.b /= samples_per_pixel;

              // Convert the average color value to an integer in the range 0 to 255 and write it to the framebuffer
              int index = 3 * (y * width + x);
              rgb_framebuffer_target[index] = std::min(std::max(static_cast<int>(pixel_color.r * 255), 0), 255);
              rgb_framebuffer_target[index + 1] = std::min(std::max(static_cast<int>(pixel_color.g * 255), 0), 255);
              rgb_framebuffer_target[index + 2] = std::min(std::max(static_cast<int>(pixel_color.b * 255), 0), 255);
          }
      }
      
/*
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
    */
  }

  Rasterizer::~Rasterizer() { }


}// CGL
