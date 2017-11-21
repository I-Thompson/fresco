REM ***********************************************************************
REM      
REM     Copyright (c) 2012, Lawrence Livermore National Security, LLC.
REM                         Produced at the Lawrence Livermore National
REM                         Laboratory.
REM                         Written by Ian Thompson, I-Thompson@llnl.gov
REM      
REM     LLNL-CODE-XXXXX All rights reserved.
REM 
REM     Copyright 2012, I.J. Thompson
REM      
REM     This file is part of FRESCO.
REM 
REM     FRESCO is free software: you can redistribute it and/or modify it
REM     under the terms of the GNU General Public License as published by
REM     the Free Software Foundation, either version 3 of the License, or
REM     (at your option) any later version.
REM      
REM     FRESCO is distributed in the hope that it will be useful, but
REM     WITHOUT ANY WARRANTY; without even the implied warranty of
REM     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
REM     GNU General Public License for more details.
REM      
REM     You should have received a copy of the GNU General Public License
REM     along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
REM      
REM     OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
REM     LICENSE
REM      
REM     Our Preamble Notice
REM 
REM       A. This notice is required to be provided under our contract
REM          with the U.S. Department of Energy (DOE). This work was
REM          produced at the Lawrence Livermore National Laboratory under
REM          Contract No. DE-AC52-07NA27344 with the DOE.
REM       B. Neither the United States Government nor Lawrence Livermore
REM          National Security, LLC nor any of their employees, makes any
REM          warranty, express or implied, or assumes any liability or
REM          responsibility for the accuracy, completeness, or usefulness
REM          of any information, apparatus, product, or process disclosed,
REM          or represents that its use would not infringe privately-owned
REM          rights.
REM       C. Also, reference herein to any specific commercial products,
REM          process, or services by trade name, trademark, manufacturer
REM          or otherwise does not necessarily constitute or imply its
REM          endorsement, recommendation, or favoring by the United States
REM          Government or Lawrence Livermore National Security, LLC. The
REM          views and opinions of authors expressed herein do not
REM          necessarily state or reflect those of the United States
REM          Government or Lawrence Livermore National Security, LLC, and
REM          shall not be used for advertising or product endorsement
REM          purposes.
REM         
REM     The precise terms and conditions for copying, distribution and
REM     modification are contained in the file COPYING.
REM         
REM ***********************************************************************
ftn95 fresco_std.f /optimize /fpp
ftn95 sfresco_std.f /optimize /fpp
pause
