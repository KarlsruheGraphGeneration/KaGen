#pragma once

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

class CairoPainter {
 public:
  CairoPainter(float x1 = 0, float y1 = 0, float x2 = 1, float y2 = 1) {
    area_bounds_[0][0] = x1;
    area_bounds_[0][1] = y1;

    area_bounds_[1][0] = x2;
    area_bounds_[1][1] = y2;

    for (uint d = 0; d < 2; ++d) {
      area_offset_[d] = area_bounds_[1][d] -
                        area_bounds_[0][d];  // offset bounds in middle of image

      img_bounds_[0][d] = 0;  // image starts at 0
      img_bounds_[1][d] =
          img_bounds_[0][d] + 3 * (area_bounds_[1][d] - area_bounds_[0][d]) *
                                  RESOLUTION;  // 9 quadrants
    }

    cs_ = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
                                      static_cast<int>(imgDim(0)),
                                      static_cast<int>(imgDim(1)));
    cr_ = Cairo::Context::create(cs_);

    // draw background white
    setColor(1, 1, 1);
    cr_->paint();

    cr_->set_line_width(1.0);
    setColor(0, 0, 0);

    // set font options
    cr_->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                          Cairo::FONT_WEIGHT_NORMAL);
    cr_->set_font_size(FONT_SIZE);

    // draw bounds
    drawRectangle(area_bounds_[0][0], area_bounds_[0][1], area_bounds_[1][0],
                  area_bounds_[1][1]);
    stroke();
  }

  CairoPainter(const CairoPainter &a) {
    area_bounds_[0][0] = a.area_bounds_[0][0];
    area_bounds_[0][1] = a.area_bounds_[0][1];
    area_bounds_[1][0] = a.area_bounds_[1][0];
    area_bounds_[1][1] = a.area_bounds_[1][1];

    img_bounds_[0][0] = a.img_bounds_[0][0];
    img_bounds_[0][1] = a.img_bounds_[0][1];
    img_bounds_[1][0] = a.img_bounds_[1][0];
    img_bounds_[1][1] = a.img_bounds_[1][1];

    area_offset_[0] = a.area_offset_[0];
    area_offset_[1] = a.area_offset_[1];

    cs_ = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
                                      static_cast<int>(imgDim(0)),
                                      static_cast<int>(imgDim(1)));
    cr_ = Cairo::Context::create(cs_);

    // set font options
    cr_->select_font_face("serif", Cairo::FONT_SLANT_NORMAL,
                          Cairo::FONT_WEIGHT_NORMAL);
    cr_->set_font_size(FONT_SIZE);

    // copy a's surface
    cr_->save();
    cr_->set_source(a.cs_, 0, 0);
    cr_->paint();

    cr_->restore();
  }

  void setColor(float r, float g, float b, float alpha = 1.0) {
    cr_->set_source_rgba(r, g, b, alpha);
  }

  // draw objects
  void drawPoint(float x, float y, uint id = ~((uint)0)) {
    cr_->arc(translatePoint(x, 0), translatePoint(y, 1), 5, 0, 2 * M_PI);
    cr_->fill();

    if (id != ~((uint)0)) {
      cr_->move_to(translatePoint(x, 0) + 7, translatePoint(y, 1) + 7);
      cr_->set_font_size(10);
      cr_->show_text(std::to_string(id));
      cr_->set_font_size(FONT_SIZE);
    }
  }

  void drawLine(float x1, float y1, float x2, float y2, bool dashed = false) {
    cr_->save();

    if (dashed) cr_->set_dash(std::vector<double>({2, 2}), 0);

    cr_->move_to(translatePoint(x1, 0), translatePoint(y1, 1));
    cr_->line_to(translatePoint(x2, 0), translatePoint(y2, 1));
    stroke();

    cr_->unset_dash();
    cr_->restore();
  }

  void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3,
                    bool dashed = false) {
    drawLine(x1, y1, x2, y2, dashed);
    drawLine(x2, y2, x3, y3, dashed);
    drawLine(x3, y3, x1, y1, dashed);

    stroke();
  }

  void drawRectangle(float x1, float y1, float x2, float y2,
                     bool dashed = false) {
    cr_->rectangle(translatePoint(x1, 0), translatePoint(y1, 1),
                   translateLength(x2 - x1, 0), translateLength(y2 - y1, 1));

    stroke();
  }

  void drawCircle(float x, float y, float r, bool dashed = false) {
    cr_->save();

    cr_->arc(translatePoint(x, 0), translatePoint(y, 1), 5, 0, 2 * M_PI);
    cr_->fill();

    if (dashed) cr_->set_dash(std::vector<double>({2, 2}), 0);

    cr_->arc(translatePoint(x, 0), translatePoint(y, 1), translateLength(r, 0),
             0, 2 * M_PI);
    stroke();

    cr_->unset_dash();
    cr_->restore();
  }

  void save(const std::string &file) const { cs_->write_to_png(file + ".png"); }

 private:
  float imgDim(const uint dim) {
    return img_bounds_[1][dim] - img_bounds_[0][dim];
  }

  float translatePoint(float in, uint dim) {
    return ((area_offset_[dim] + in - area_bounds_[0][dim]) /
            (3 * (area_bounds_[1][dim] - area_bounds_[0][dim]))) *
           imgDim(dim);
  }

  float translateLength(float in, uint dim) {
    return (in / (3 * (area_bounds_[1][dim] - area_bounds_[0][dim]))) *
           imgDim(dim);
  }

  void stroke() { cr_->stroke(); }

 private:
  float area_bounds_[2][2];
  float img_bounds_[2][2];
  float area_offset_[2];

  const uint RESOLUTION = 1000;
  const uint PADDING = 10;
  const uint FONT_SIZE = 15;

  Cairo::RefPtr<Cairo::ImageSurface> cs_;
  Cairo::RefPtr<Cairo::Context> cr_;
};