// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef CAMERA_MANIP_ROT_BPA_H
#define CAMERA_MANIP_ROT_BPA_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera_manip.h"

class Camera_manip_rot_bpa : public Camera_manip {
public:
  Camera_manip_rot_bpa(Camera& camera);

protected:
  virtual void mouse_press_event(QMouseEvent* e) override;
  virtual void mouse_move_event(QMouseEvent* e) override;
  virtual void mouse_release_event(QMouseEvent* e) override;
  virtual void resize(int w, int h) override;

private:
  int m_vp_width, m_vp_height;
};

#endif
