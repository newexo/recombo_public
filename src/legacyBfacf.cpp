#include "legacyBfacf.h"
#include "bfacfMove.h"
#include "edgePool.h"
#include "clkRecombination.h"
#include "clkValidation.h"
#include "clkLattice.h"
#include "clkTopology.h"

#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);


bool add_edge_to_knot(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end)
{
   EdgePtr ep = (EdgePtr) calloc(1, sizeof (Edge)); // create a new edge

   sub_ivector(ep->increment, end, start); // increment is (end - start)
   ep->dir = clk_direction(ep->increment); // direction associated with increment
   copy_ivector(ep->start, start); // starting location of this edge
   ep->frozen = false; // edge is free to move
   ep->comp = (void *) comp;

   if (comp->first_edge)
   { // some edges have already been added
      comp->last_edge->next = ep;
      ep->prev = comp->last_edge;
      comp->last_edge = ep;
   }
   else
   { // this is the first edge
      comp->first_edge = comp->last_edge = ep;
   }

   comp->first_edge->prev = comp->last_edge;
   comp->last_edge->next = comp->first_edge;
   comp->nedges++;
   clkp->nedges_total++;

   return clk_check_increment(ep->increment);
}

void bfacf_init_pool_and_lattice(CubicLatticeKnotPtr clkp, bool ignore_mid, int poolsize, ivector min, ivector max)
{
   ComponentCLKPtr comp = clkp->fcomp;

   clkp->auto_recentre = true;

   // create a `swimming space' (a pool) for our edges
   // this is to avoid frequent memory allocation / deallocation

   clkp->poolsize = poolsize;
   clkp->edgepool = (EdgePtr *) calloc(clkp->poolsize, sizeof (EdgePtr));
   int index = 0;
   int ID = 0;

   comp = clkp->fcomp;
   while (comp)
   {
      comp->ID = ID++;
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (index > poolsize)
         {
            printf("yikes! %d %d\n", index, poolsize);
            exit(11);
         }
         clkp->edgepool [index++] = ep;
         ep = ep->next;
      }
      comp = comp->next;
   }

   while (index < poolsize)
      clkp->edgepool [index++] = (EdgePtr) calloc(1, sizeof (Edge)); // new, unused edge

   for (index = 0; index < clkp->poolsize; index++)
   {
      clkp->edgepool [index]->locpool = index;
      clkp->edgepool [index]->ID = index;
   }

   clkp->lattice = (char *) calloc(LATTICE_TOTAL_SIZE, sizeof (char));
   ivector mid;
   midpoint(mid, min, max);

   if (ignore_mid)
      set_ivector(mid, 0, 0, 0);

   init_lattice(clkp, mid);


   comp = clkp->fcomp;
   while (comp)
   {
      ivector closing_incr;
      sub_ivector(closing_incr, comp->first_edge->start, comp->last_edge->start);
      // if the edge that would close the knot is not valid,
      // then we have an open-ended string.  
      if (!clk_check_increment(closing_incr))
      {
         comp->flags |= COMPONENT_CLK_FLAG_OPEN;
         delete_Edge(clkp, comp, comp->last_edge);
         freezeEdge(clkp, comp->first_edge);
         freezeEdge(clkp, comp->last_edge);
         comp->first_edge->frozen = comp->last_edge->frozen = true;
      }

      comp = comp->next;
   }
}


EdgePtr clk_get_edge(CubicLatticeKnotPtr clkp, ivector start)
{
   ComponentCLKPtr comp = clkp->fcomp;
   add_ivector(start, start, clkp->loffset);
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(ep->start, start))
         {
            return ep;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL;
}


// from clk_util.cpp

int clk_count_edges(EdgePtr ep)
{
   int kount = 1;
   EdgePtr start = ep;
   ep = ep->next;
   while (ep != start)
   {
      kount++;
      ep = ep->next;
   }
   return kount;
}




int clk_arc_length_distance(EdgePtr ep1, EdgePtr ep2)
{
   ComponentCLKPtr comp1 = (ComponentCLKPtr) ep1->comp;
   ComponentCLKPtr comp2 = (ComponentCLKPtr) ep2->comp;
   if (comp1 != comp2)
      return ARC_LENGTH_DISTANCE_INFINITE;

   EdgePtr epn = ep1->next;
   EdgePtr epp = ep1->prev;
   int distance = 1;
   while (epn != ep2 && epp != ep2)
   {
      epn = epn->next;
      epp = epp->prev;
      ++distance;
   }

   return distance;
}

int clk_arc_length_distance(CubicLatticeKnotPtr clkp, ivector s1, ivector s2)
{
   EdgePtr ep1 = clk_get_edge(clkp, s1);
   if (!ep1) return ARC_LENGTH_DISTANCE_NO_EDGE1;
   EdgePtr ep2 = clk_get_edge(clkp, s2);
   if (!ep2) return ARC_LENGTH_DISTANCE_NO_EDGE2;
   return clk_arc_length_distance(ep1, ep2);
}

void clk_get_extent(ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot)
{
   if (!knot) return;
   copy_ivector(minbb, knot->fcomp->first_edge->start);
   copy_ivector(maxbb, minbb);
   ComponentCLKPtr comp = knot->fcomp;

   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         min_ivector(minbb, minbb, ep->start);
         max_ivector(maxbb, maxbb, ep->start);
         ep = ep->next;
      }
      comp = comp->next;
   }
   sub_ivector(minbb, minbb, knot->loffset);
   sub_ivector(maxbb, maxbb, knot->loffset);
}

void clk_bounding_box(int &X, int &Y, int &Z, CubicLatticeKnotPtr knot)
{
   ivector minbb, maxbb;
   clk_get_extent(minbb, maxbb, knot);
   X = maxbb [0] - minbb [0];
   Y = maxbb [1] - minbb [1];
   Z = maxbb [2] - minbb [2];
}




#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

int clk_direction(ivector incr)
{
   // return the direction associated with an increment
   switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
   {
      case NORTH:
         return MOVE_NORTH;
         break;
      case EAST:
         return MOVE_EAST;
         break;
      case WEST:
         return MOVE_WEST;
         break;
      case SOUTH:
         return MOVE_SOUTH;
         break;
      case UP:
         return MOVE_UP;
         break;
      case DOWN:
         return MOVE_DOWN;
         break;
   }
   return MOVE_INVALID;
}

int clk_direction(ivector end, ivector start)
{
   ivector incr;
   sub_ivector(incr, end, start);
   return clk_direction(incr);
}


void clk_get_statistics(CubicLatticeKnotPtr clkp, int &n1, int &n2, int &n3, int &n4)
{
   // gather statistics on local edge configurations with an eye to possibly
   // reordering the cases
   // currently, this function is used only by the KnotPlot driver
   if (!clkp) return;
   n1 = n2 = n3 = n4 = 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         // Case 1: chosen edge is parallel to both adjacent edges,
         // all four moves are +2 moves.

         // directions of adjacent edges same as chosen edge

         if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
         {
            ++n1;
         }

            // Case 2: chosen edge is perpendicular to both adjacent edges
            // and those adjacent edges are anti-parallel to each other. 
            // One move is a -2 move, other three are +2 moves.

         else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
         {
            ++n2;
         }

            // Case 3: chosen edge is perpendicular to one adjacent edge
            // and parallel to the other adjacent edge.
            // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
            // because the case of both being parallel has been ruled out above.
            // One move is a 0 move, other three are +2 moves.

            // NOTE: test this before Case 1 above

         else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
         {
            ++n3;
         }

            // Case 4: Only possibility remaining: chosen edge is perpendicular
            // to both adjacent edges, and the adjacent edges are not anti-parallel.
            // Two moves are 0 moves and two moves are +2 moves.

         else
         {
            ++n4;
         }

         ep = ep->next;

      }
      comp = comp->next;
   }
}

bool clkp_always_turns(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (ep->dir == ep->next->dir) return false;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return true;
}


// currently, set_lattice(), fill_lattice(), invert_lattice(), dilate_lattice(), and clear_lattice()
// are used only by the KnotPlot driver to implement topological obstructions and cavities









#define ABS(X)      ((X) < 0 ? -(X) : (X))

void bfacf_set_probabilities(ComponentCLKPtr comp, double pm2, double p0, double pp2)
{
   comp->p_minus2 = pm2;
   comp->p_0 = p0;
   comp->p_plus2 = pp2;

   // probably doesn't make much of a difference,
   // but precompute the following anyway:

   comp->p_4p2 = 4.0 * comp->p_plus2;
   comp->p_03p2 = comp->p_0 + 3.0 * comp->p_plus2;
   comp->p_m23p2 = comp->p_minus2 + 3.0 * comp->p_plus2;
   comp->p_2p0 = 2.0 * comp->p_0;
   comp->p_2p02p2 = 2.0 * comp->p_0 + 2.0 * comp->p_plus2;
}

void bfacf_set_probabilities(ComponentCLKPtr comp, double z)
{
   double p_plus2 = (z * z) / (1.0 + 3.0 * z * z);
   double p_0 = (1.0 + z * z) / (2.0 * (1.0 + 3.0 * z * z));
   double p_minus2 = 1.0 / (1.0 + 3.0 * z * z);
   comp->z = z;
   bfacf_set_probabilities(comp, p_minus2, p_0, p_plus2);
}

void bfacf_set_probabilities(CubicLatticeKnotPtr knot, double z)
{
   ComponentCLKPtr comp = knot->fcomp;
   while (comp)
   {
      bfacf_set_probabilities(comp, z);
      comp = comp->next;
   }
}

void bfacf_lattice_info(CubicLatticeKnotPtr knot)
{
   int min, max;
   int minloc, maxloc;
   minloc = maxloc = 0;
   min = max = knot->alt_lattice [0];
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (knot->alt_lattice [i] > max)
      {
         max = knot->alt_lattice [i];
         maxloc = i;
      }
      if (knot->alt_lattice [i] < min)
      {
         min = knot->alt_lattice [i];
         minloc = i;
      }
   }
   printf("%d, min value of %d at %d, max value of %d at %d\n", knot->nedges_total, min, minloc, max, maxloc);
}

void pexit(char *s)
{
   fprintf(stderr, "%s\n", s);
   fflush(stderr);
   exit(102);
}
