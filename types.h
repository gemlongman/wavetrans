/*****************************************************************************
 * wavetrans - a simple wave transformation tool                             *
 * Copyright (C) 2005 Radu Rendec <rrendec@yahoo.com>                        *
 *                                                                           *
 * This program is free software; you can redistribute it and/or modify      *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation; either version 2 of the License, or         *
 * (at your option) any later version.                                       *
 *                                                                           *
 * Foobar is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 * GNU General Public License for more details.                              *
 *                                                                           *
 * You should have received a copy of the GNU General Public License         *
 * along with Foobar; if not, write to the Free Software                     *
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
 *                                                                           *
 *****************************************************************************/

/**********************************************************************/
/*                                                                    */
/* Tipuri de date intregi. Aceste tipuri sunt dependente de           */
/* platforma si ar trebui configurate automat (de exemplu cu          */
/* autoconf + automake).                                              */
/*                                                                    */
/* Definitiile de aici au fost scrise pentru GNU Linux / i386         */
/*                                                                    */
/**********************************************************************/

#ifndef _TYPES_H
#define _TYPES_H

typedef short t_s16;
typedef long long t_s64;

typedef unsigned int t_u32;
typedef int t_s32;

#endif
