// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <stdexcept>
#include <omp.h>
#include <array>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <memory>
#include <cmath>
#include <limits>
#include <execution>
#include <optional>
#include <conio.h>
#include <chrono>
#include <fstream>
#include <string_view>
#include <filesystem>

#include "boost/container/static_vector.hpp"
#include "tinyxml2.h"
#include "cxxopts.hpp"
#include "GL\freeglut.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPLYReader.h"
#include "vtkSelectEnclosedPoints.h"
#include "vtkSmartPointer.h"
#include "vtkSurfaceReconstructionFilter.h"
#include "vtkContourFilter.h"
#include "vtkReverseSense.h"
#include "vtkMarchingCubes.h"
#include "vtkVoxelModeller.h"
#include "vtkImageData.h"
 
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2);