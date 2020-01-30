C=========================================================================
C DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
C and analysis tool for aircraft with propulsors in ducted fan
C configurations.
C 
C This software was developed under the auspices and sponsorship of the
C Tactical Technology Office (TTO) of the Defense Advanced Research
C Projects Agency (DARPA).
C 
C Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
C
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version.
C 
C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C 
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
C
C Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
C Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
C
C=========================================================================
C
C Version 070-ES2
C PORTFAC variable added. Scales portrait oriented windows vs landscape
C Version 070-ES3
C PAR & LSLOPE variables added for compatibility with xrotor subs
C
C=========================================================================

      LOGICAL LPLOT, LCURS, LLAND, LEGEND, LGRID, LSLOPE
C
      COMMON /PR_PLT/ SIZE, SCRNFR, 
     &                CHGT, SHGT, CSIZE, PAR,
     &                XPAGE,YPAGE,
     &                XWIND,YWIND,
     &                XMARG,YMARG,
     &                XMIN,YMIN,
     &                XMAX,YMAX,
     &                XYOFF(2),XYFAC(2),
     &                PORTFAC
C
      COMMON /PI_PLT/ IDEV, IDEVPS, IPSLU, NCOLOR
      COMMON /PL_PLT/ LPLOT,LCURS,LLAND,LEGEND,LGRID,LSLOPE
C
C-------------------------------------------------------------------------
C

