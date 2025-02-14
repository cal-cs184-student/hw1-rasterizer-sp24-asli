#include "texture.h"
#include "CGL/color.h"
#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
      double level;
      int level0{};
      int level1{};
      int weight{};
      double level_f;

      // Determine the sampling method based on `lsm`
      if (sp.lsm == L_ZERO) {
        level = 0; // Use the highest resolution mipmap level
        level1 = 0;
        level0 = 0;
    } else if (sp.lsm == L_NEAREST) {
        level = get_level(sp); // Compute the nearest mipmap level
        level = round(level); // Clamp the level
        level1 = round(level);
        level0 = round(level);
        //printf("dd%f",level);
    } else if (sp.lsm == L_LINEAR) {
        // Interpolate between two mipmap levels
        level_f = get_level(sp);
        //printf("gg%f",level_f);
        level0 = std::floor(level_f);
        level1 = std::ceil(level_f);
        //level0 = std::max(0, std::min(level0, (int)mipmap.size() - 1));
        //level1 = std::max(0, std::min(level1, (int)mipmap.size() - 1));

        //Color color0 = sample_nearest(sp.p_uv, level0);
        //Color color1 = sample_nearest(sp.p_uv, level1);
        weight = level_f - level0;
        
      }

      // Sample from the determined mipmap level
      if (sp.psm == P_NEAREST) {
          Color color0 = sample_nearest(sp.p_uv, level0);
          Color color1 = sample_nearest(sp.p_uv, level1);
          return color0 * weight+color1 * (1-weight);
      }
      else if (sp.psm == P_LINEAR) {
          Color color0 = sample_bilinear(sp.p_uv, level0);
          Color color1 = sample_bilinear(sp.p_uv, level1);
          return color0 * weight + color1 * (1-weight);
      }
      //else {
          //return Color(1, 0, 1); // Magenta for invalid sampling method
      //}
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    // Step 1: Compute Derivatives
      Vector2D duv_dx = sp.p_dx_uv - sp.p_uv;
      Vector2D duv_dy = sp.p_dy_uv - sp.p_uv;

      // Step 2: Calculate Mipmap Level
      // Scale the derivatives by the texture size
      duv_dx.x *= width;
      duv_dx.y *= width;
      duv_dy.x *= height;
      duv_dy.y *= height;

      // Choose the largest magnitude between the two derivatives
      float L = std::max(std::sqrt(duv_dx.x * duv_dx.x + duv_dx.y * duv_dx.y),
          std::sqrt(duv_dy.x * duv_dy.x + duv_dy.y * duv_dy.y));

      // Calculate level as the log2 of the largest derivative
      return L > 0 ? std::min(static_cast<double>(std::log2(L)), (double)mipmap.size()) : 0;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }


  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
       // Step 1: Check for a valid mipmap level
      if (level < 0 || level >= mipmap.size()) {
          return Color(1, 0, 1); // magenta for an invalid level
      }

      // Step 2: Retrieve the appropriate mipmap level (assuming mipmap is a std::vector<MipLevel>)
      auto& mip = mipmap[level];

      // Step 3: Scale the UV coordinates
      float u = uv.x * (mip.width - 1);
      float v = uv.y * (mip.height - 1);

      // Step 4: Find the nearest texel coordinates
      int x = static_cast<int>(std::round(u));
      int y = static_cast<int>(std::round(v));

      // Step 5: Retrieve the texel color
      Color color = mip.get_texel(x, y);

      // Step 6: Return the texel color
      return color;
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
      // Step 1: Validate Mipmap Level
      if (level < 0 || level >= mipmap.size()) {
          return Color(1, 0, 1); // Magenta for an invalid level
      }

      // Step 2: Retrieve Mipmap Level
      auto& mip = mipmap[level];

      // Step 3: Scale UV Coordinates
      float u = uv.x * (mip.width - 1);
      float v = uv.y * (mip.height - 1);

      // Step 4: Calculate Texel Indices
      int x0 = static_cast<int>(floor(u));
      int y0 = static_cast<int>(floor(v));
      int x1 = static_cast<int>(ceil(u));
      int y1 = static_cast<int>(ceil(v));

      // Step 5: Retrieve Texel Colors
      Color c00 = mip.get_texel(x0, y0);
      Color c10 = mip.get_texel(x1, y0);
      Color c01 = mip.get_texel(x0, y1);
      Color c11 = mip.get_texel(x1, y1);

      // Step 6: Interpolate Horizontally
      float s = u - x0;
      Color top = c00 * (1 - s) + c10 * s;
      Color bottom = c01 * (1 - s) + c11 * s;

      // Step 7: Interpolate Vertically
      float t = v - y0;
      Color finalColor = top * (1 - t) + bottom * t;

      // Step 8: Return Final Color
      return finalColor;
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
