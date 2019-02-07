/*  This file is part of the PopART IBM.

    The PopART IBM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The PopART IBM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the PopART IBM.  If not, see <http://www.gnu.org/licenses/>.
 */

/* compat.h handles issues related to compatibility between different compilers. */

#ifndef COMPAT_H
#define COMPAT_H

#ifdef _WIN32
   #define THISOS 0
   #define POPEN _popen
   #define PCLOSE _pclose
   #define _CRT_SECURE_NO_WARNINGS 1
   #define _CRT_SECURE_NO_DEPRECATION 1
#else
   #define THISOS 1
   #define POPEN popen
   #define PCLOSE pclose
#endif

#endif
