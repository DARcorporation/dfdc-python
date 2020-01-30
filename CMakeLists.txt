cmake_minimum_required(VERSION 3.12)
project(dfdc Fortran)

enable_language(Fortran)
add_library(dfdc SHARED
        src/adjpanl.f
        src/aero.f
        src/aic.f
        src/airio.f
        src/atmo.f
        src/dfdcsubs.f
        src/forces.f
        src/gauss.f
        src/geom.f
        src/geutil.f
        src/grdutils.f
        src/inigrd.f
        src/lamp.f
        src/pnsubs.f
        src/qaic.f
        src/rotoper.f
        src/sgutil2.f
        src/solve.f
        src/spline.f
        src/system.f
        src/userio.f
        src/vels.f
        src/viscvel.f
        src/wakesubs.f
        src/xutils.f)