/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef SCHED_H
#define SCHED_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern void free_comm_schedule(int **sched, int npes);
extern void empty_comm_schedule(int **sched, int npes);
extern int **make_comm_schedule(int npes);
extern void fill_comm_schedule(int **sched, int npes);
extern int check_comm_schedule(int **sched, int npes);
extern void invert_comm_schedule(int **sched, int npes);
extern void sort_comm_schedule(int **sched, int npes, int sort_pe);
extern void print_comm_schedule(int **sched, int npes);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* SCHED_H */
