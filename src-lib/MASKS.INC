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
C---- dot-pattern masks for use with PLGRID, NEWPAT, etc.
C
C     mask0:  _________________________   (solid)
C         1:  .........................
C         2:  . . . . . . . . . . . . .
C         3:  .   .   .   .   .   .   .
C         4:  .       .       .       .
C
      DATA LMASK0, LMASK1, LMASK2, LMASK3, LMASK4
     &   / -1    , -21846, -30584, -32640, -32768 /
C
C      DATA LMASK0,     LMASK1,     LMASK2,     LMASK3,     LMASK4
C     &   /Z'FFFFFFFF',Z'AAAAAAAA',Z'CCCCCCCC',Z'88888888',Z'21084210'/



