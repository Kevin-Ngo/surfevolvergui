/* vor2fe.c
   Generates large number of random voronoi cells in unit torus.
   Adapted from vtor2d.c to get initial few cells, 
   then does Bowyer algorithm to add rest, from vor2plane.c.
   Also uses high-quality kb_drand() random number generator.
   Uses random selection of Delaunay triangles to insert incremental
     random seed.  Selection done by generating future event time for
     each triangle when it is generated, and keeping sorted heap for
     finding next time.

 Usage: vor2fe [-sS] [-nN] [-f] > datafile.fe
    Options:  -sS   to set random number generator seed to integer S
              -nN   do N cells, default 1000.
              -f    do foam model with fixed volumes.
              -m r c   multiprocessor model, rows and cols.

*/

/*  
   Programmer: Ken Brakke, Susquehanna University, brakke@susqu.edu
   Begun: April 24, 1999
   Last modified: May 16, 1999
*/
#define XVORDEBUG

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#ifndef M_PI
#define M_PI pi
#endif
double pi;  // stupid Microsoft math.h doesn't have M_PI
int seed = 1;

// random number functions from kb_drand.c
extern double kb_drand(void);   
extern int kb_initr(int);
#define frand kb_drand
#define srand kb_initr
void check_orders(void);
void do_datafile(void);
void redo_edges(void);

int phases;  /* number of phases, if used */
int phaseflag; /* whether phases used */
int foam_flag; /* whether to do foam model */

#define MAXSIDES 25
struct cell { double x[2];  /* seed coords */
              int v[MAXSIDES];    /* vertex list */
              int vwraps[MAXSIDES][2]; /* wraps to vertices */
              struct vertex *vv[MAXSIDES];  /* link into vertex net */
              int e[MAXSIDES];    /* edge list */
              int count;   /* number of vertices */
              int ccount;  /* number of neighbor cells */
              int proc;    /* processor number */
            } *cells;
int ccount; /* number of cells done */ 

struct vertex { double x[2]; /* vertex coords */
                double rr;   /* square radius of void */
                double area; /* area of Delaunay triangle */
                double eventtime; /* when next seed appears in triangle */
                int cell[3]; /* nbr cells */
                int cwraps[3][2];  /* wraps to nbr cells */
                struct vertex *nbrv[3]; /* nbr vertices, in counterclockwise order */
                                          /* cell[i] before nbrv[]i */
                int nwrap[3][2];   /* wraps to neighbors */
                struct edge *nbre[3];  /* edges to neighbors */
                struct vertex *nextfree;  /* for freelist */
                int flag;  /* VFREE or VUSED */
                int num;   /* number in datafile */
                int proc;    /* processor number */
              } *verts;
#define VFREE 0
#define VUSED 1
int vcount;  /* vertices done so far */
struct vertex **tree[30];  /* for next event heap */
struct vertex *newvlist[30];
struct vertex *oldvlist[30]; // to be deleted
struct vertex *freehead;  /* vertex freelist */
int newvcount; /* number of new vertices */
int oldvcount; /* number deleted */

struct edge { int v[2];   /* endpoints */
              int wrap[2]; /* wraps for each coord */
              int proc;    /* processor number */
            } *edges;
int ecount;  /* edges done so far */
int edgealloc; /* structures allocated */

int N = 1000; /* seeds */
int maxlevel;
int multirows = 1;
int multicols = 1;
int multiflag = 0;

int stamp = 1;

#define BIGTIME 1e30
struct vertex bigtime;  /* placeholder in heap */
double now;

int calc_voronoi(int);
int assemble_net(void);
void do_datafile(void);

void setup(void);
void findtree(struct vertex *v, struct vertex *parent,struct cell *c);
void addpoint(void);
void batch(void);

/* assign processor numbers */
void assign_procs()
{ int k;
  struct vertex *v;
  struct edge *e;
  struct cell *c;

  /* vertices */
  for ( k = 0, v = verts ; k < 2*N+100 ; k++,v++ )
  { if ( v->flag & VUSED )
    { int row,col;
 
      row = (int)floor(v->x[0]*multirows);
      col = (int)floor(v->x[1]*multicols + ((row & 1) ? 0.5 : 0.0));
      v->proc = 1 + row*multicols + col;
    }
  }
  /* edges */
  for ( k = 0, e = edges ; k < ecount ; k++,e++ )
  { int row,col;
    double x = (verts[e->v[0]].x[0] + verts[e->v[1]].x[0])/2;
    double y = (verts[e->v[0]].x[1] + verts[e->v[1]].x[1])/2;

    x += 0.5*e->wrap[0];
    y += 0.5*e->wrap[1];

    row = ((int)floor(x*multirows) + multirows) % multirows;
    col = ((int)floor(y*multicols + ((row & 1) ? 0.5 : 0.0)) + multicols) % multicols;
    e->proc = 1 + row*multicols + col;
  }
  /* faces */
  for ( k = 0, c = cells ; k < N ; k++,c++ )
  {
    c->proc = edges[abs(c->e[0])-1].proc;
  }
}

/* add triangle to heap */
void add_triangle(struct vertex *v)
{ int m,n;
  double t;

  t = v->eventtime = now + (-log(frand())/v->area);
  for ( m = maxlevel, n = v-verts ; m > 0 ; m-- )
  { tree[m][n] = v;
    if ( t >= tree[m][n^1]->eventtime ) return;
    n>>=1;
  }
  tree[0][0] = v;
}

/* delete triangle from heap */
void del_triangle(struct vertex *v)
{ int m,n;

  v->eventtime = BIGTIME;
  for ( m = maxlevel, n = v-verts ; m > 0 ; m--, n >>= 1 )
   if ( tree[m-1][n>>1] == v ) 
     { if ( tree[m][n]->eventtime < tree[m][n^1]->eventtime )
          tree[m-1][n>>1] = tree[m][n];
       else
          tree[m-1][n>>1] = tree[m][n^1];
     }
   else return;
}

/* finding which triangle has next event */
struct vertex *find_next(void)
{ struct vertex *v;
  
   v = tree[0][0];
   now = v->eventtime;

   del_triangle(v);

   return v;
}

/* sees if vertex already exists */
/* called only from calc_voronoi() in initial setup */
int swaptemp;
#define swap(a,b) {swaptemp=a;a=b;b=swaptemp;}
int find_vertex(a,b,c,newv)
int a,b,c;  /* cell numbers of neighbors */
struct vertex *newv; /* location */
{ int k,j;
  struct vertex *v;
  
  /* sort */
  if ( a > b ) { swap(a,b); swap(newv->cell[0],newv->cell[1]);
                 for(j=0;j<2;j++) swap(newv->cwraps[0][j],newv->cwraps[1][j]);
               }
  if ( b > c ) { swap(b,c); swap(newv->cell[1],newv->cell[2]);
                 for(j=0;j<2;j++) swap(newv->cwraps[1][j],newv->cwraps[2][j]);
               }
  if ( a > b ) { swap(a,b); swap(newv->cell[0],newv->cell[1]);
                 for(j=0;j<2;j++) swap(newv->cwraps[0][j],newv->cwraps[1][j]);
               }

  /* linear search */
  for ( k = 0, v = verts ; k < vcount ; k++, v++ )
   if ( (a==v->cell[0]) && (b==v->cell[1]) && (c==v->cell[2])
       && (fabs(v->x[0]-newv->x[0])+fabs(v->x[1]-newv->x[1]) < 0.00000001))
       return k;

  /* here is where algorithm failure shows up */
  if ( vcount >= 2*N ) 
  {fprintf(stderr,"."); return -1; }

  /* not found, so add */
  freehead = v->nextfree; v->nextfree = NULL; v->flag = VUSED;
  /* coordinates and cells and wraps were already set up in calc_voronoi() */
  vcount++;
  return k;
}

void print_help()
{
  fprintf(stderr," vor2fe: Creates random Voronoi tessellation in 2d unit flat torus.\n");
  fprintf(stderr," Usage: vor2fe [-sS] [-nN] [-f] > datafile.fe \n");
  fprintf(stderr,"    Options:  -sS   to set random number generator seed to integer S\n");
  fprintf(stderr,"              -nN   do N cells, default 1000.\n");
  fprintf(stderr,"              -f    do foam model with fixed volumes.\n");
  fprintf(stderr,"              -m r c   multiprocessor model, rows and cols.\n");
}


main(argc,argv)
int argc;
char *argv[];
{
  srand(1);

  if ( argc < 2 ) 
  { fputs("ERROR: Need number of cells on command line.\n",stderr); 
    print_help();
    exit(1); 
  }
  for ( ; argv[1] && argv[1][0] == '-'; argv++, argc-- )
  { switch ( argv[1][1] )
    {
      case 's': srand(atoi(argv[1]+2)); break;
      case 'n': N = atoi(argv[1]+2); break;
      case 'f': foam_flag = 1; break;
      case 'm': multiflag = 1; 
                multirows = atoi(argv[2]);
                multicols = atoi(argv[3]);
                if ( multirows<1 || multicols<1 )
                { print_help(); exit(1); }
                break;

      default : print_help(); exit(1);
    }
  } 

  if ( N < 2 ) 
  { fputs("ERROR: Too few cells.  Need at least 2.\n",stderr);
    exit(2); 
  }
 
  cells = (struct cell*)calloc(N+1+100,sizeof(struct cell));
  edgealloc = 10*N;
  edges = (struct edge*)calloc(edgealloc,sizeof(struct edge));
  verts = (struct vertex*)calloc(2*N+100,sizeof(struct vertex));

  batch();
  if ( multiflag )
    assign_procs();
  do_datafile();

  return 0;
}

void batch()
{ int level,bins,i;

  memset(cells,0,(N+1)*sizeof(struct cell));
  memset(edges,0,edgealloc*sizeof(struct edge));
  memset(verts,0,(2*N+100)*sizeof(struct vertex));

  /* vertex freelist */
  freehead = verts;
  for ( i = 0 ; i < 2*N+99 ; i++ ) verts[i].nextfree = verts+i+1;
  ccount = vcount = ecount = 0;

  /* next event heap */
  bigtime.eventtime = BIGTIME;
  for ( level = 0, bins = 1  ; ; level++,bins*=2 )
  { tree[level] = (struct vertex **)calloc(bins,sizeof(struct vertex*));
    for ( i = 0 ; i < bins ; i++ ) tree[level][i] = &bigtime;
    if ( bins > 2*N+30 ) break; /* have enough */
  }
  maxlevel = level; 

  setup();         
                                                                  

#ifdef VORDEBUG
  check_orders();
#endif

  for ( ; ccount < N ; )
  {  addpoint();           
#ifdef VORDEBUG
     check_orders();
     getchar();
#endif
  }
  if ( N > 99 )
    redo_edges();
 
  for ( i = 0 ; i < 30 ; i++ ) 
    if ( tree[i] ) free((char*)tree[i]);

}

void setup(void)
{ int i;
  int n;
  
  for ( i = 0; i < 100 ; i++ )
  { cells[i].x[0] = frand();
    cells[i].x[1] = frand();
    cells[i].count = 0;
  }
  ecount = 0;
  vcount = 0;
  n = N > 100 ? 100 : N;
  calc_voronoi(n);
  assemble_net();

  for ( i = 0 ; i < ccount ; i++ )
  { for ( n = 0; n < cells[i].count ; n++ )
      cells[i].vv[n] = verts + cells[i].v[n];
  }

}

/* calc_voronoi() finds vertices around each seed. Sets cells and
   coordinates for each vertex, but not neighbor vertices.
   Sets cells->v[] in order around cells.
*/
int calc_voronoi(int ncells)
{ int i,j,k;
  struct vertex *tempnext;

  /* find voronoi cell of each seed */
  for ( i = 0 ; i < ncells ; i++ )
  { struct cell *ci = cells + i;
    double rr,minrr = 1e12;
    int thisj,nextj,minj = -1;
    double dx[2];
    double xt[2];
    int wx,wy,minwx,minwy;
    int thiswrap[2],nextwrap[2];

    /* first, find nearest neighbor */
    for ( j = 0 ; j < ncells ; j++ )
    { //if ( j == i ) continue;
      for ( wx = -1 ; wx <= 1 ; wx++ )
        for ( wy = -1 ; wy <= 1 ; wy++ )
        { // testing all possible wraps
          if ( (j == i) && (wx == 0) && (wy == 0) )
            continue; // don't do self
          dx[0] = ci->x[0] - (cells[j].x[0]+wx);
          dx[1] = ci->x[1] - (cells[j].x[1]+wy);
         
          rr = dx[0]*dx[0] + dx[1]*dx[1];
          if ( rr < minrr )
          { minrr = rr;
            minj = j;
            minwx = wx;
            minwy = wy;
          }
        }
    }

   /* now go around clockwise */
   thisj = minj; thiswrap[0] = minwx; thiswrap[1] = minwy;
   do
   { 
     double cx[2];  /* midpoint of segment */
     double u[2];   /* direction of perp bisector */
     double t,mint = 1e12;
     double denom;
     struct cell *ct;
     int newv;

     ct = cells + thisj;
     for ( k = 0 ; k < 2 ; k++ )  /* unwrap second vertex */
       xt[k] = ct->x[k] + thiswrap[k]; 

     for ( k = 0 ; k < 2 ; k++ )
       cx[k] = (ci->x[k] + xt[k])/2;  // midpoint between seeds
     u[0] = xt[1] - ci->x[1];  // perp bisector direction
     u[1] = -xt[0] + ci->x[0];
     for ( j = 0 ; j < ncells ; j++ )
     { double xj[2];
       for ( wx = -1 ; wx <= 1 ; wx++ )  // try all possible wraps
         for ( wy = -1 ; wy <= 1 ; wy++ )
         { if ( (j == i) && (wx == 0) && (wy == 0) ) continue;
           if ( (j == thisj) && (wx == thiswrap[0]) && (wy == thiswrap[1]) )
             continue;
           /* unwrap trial vertex */
           xj[0] = cells[j].x[0] + wx;
           xj[1] = cells[j].x[1] + wy;
      
           for ( k = 0, denom = 0.0 ; k < 2 ; k++ )
             denom += u[k]*(xj[k] - ci->x[k]);
           if ( denom <= 0.0 ) continue;
           for ( k = 0, t = 0.0 ; k < 2 ; k++ )
             t += (xj[k]*xj[k] - ci->x[k]*ci->x[k])/2 - (xj[k] - ci->x[k])*cx[k];
           t /= denom;

           if ( t < mint )
           { nextj = j;
             mint = t;
             nextwrap[0] = wx;
             nextwrap[1] = wy;
           }
        }
     }

     /* now have nextj as third cell */
       newv = vcount;  // tentative new vertex
       tempnext = verts[newv].nextfree;
       memset(&verts[newv],0,sizeof(struct vertex));
       verts[newv].nextfree = tempnext;
       for ( k = 0 ; k < 2 ; k++ )
         verts[newv].x[k] = cx[k] + mint*u[k];
       verts[newv].rr = 
           (verts[newv].x[0] - ci->x[0]) * (verts[newv].x[0] - ci->x[0])
         + (verts[newv].x[1] - ci->x[1]) * (verts[newv].x[1] - ci->x[1]);
                     // get into fundamental domain
       verts[newv].cell[0] = i;
       verts[newv].cell[1] = thisj;
       verts[newv].cell[2] = nextj;
 
       if ( verts[newv].x[0] < 0.0 ) 
       { verts[newv].x[0]++; ci->vwraps[ci->count][0] = -1; 
         verts[newv].cwraps[0][0] = 1;
       }
       else if ( verts[newv].x[0] > 1.0 ) 
       { verts[newv].x[0]--; ci->vwraps[ci->count][0] = 1; 
         verts[newv].cwraps[0][0] = -1;
       }
       if ( verts[newv].x[1] < 0.0 ) 
       { verts[newv].x[1]++; ci->vwraps[ci->count][1] = -1; 
         verts[newv].cwraps[0][1] = 1; 
       }
       else if ( verts[newv].x[1] > 1.0 ) 
       { verts[newv].x[1]--; ci->vwraps[ci->count][1] = 1; 
         verts[newv].cwraps[0][1] = -1;
       }
       verts[newv].cwraps[1][0] = verts[newv].cwraps[0][0]-thiswrap[0];
       verts[newv].cwraps[1][1] = verts[newv].cwraps[0][1]-thiswrap[1];
       verts[newv].cwraps[1][0] = verts[newv].cwraps[0][0]-nextwrap[0];
       verts[newv].cwraps[1][1] = verts[newv].cwraps[0][1]-nextwrap[1];
         
     // see if already exists
     newv = find_vertex(i,thisj,nextj,verts+newv);
     if ( newv < 0 ) 
       return 0;  /* failure */
     ci->v[ci->count++] = newv;

     thisj = nextj; thiswrap[0] = nextwrap[0]; thiswrap[1] = nextwrap[1];
    }
    while ( (nextj != minj) || (nextwrap[0] != minwx) || (nextwrap[1] != minwy) );  /* end single cell loop */
    ccount++;
  } /* end main loop */
  return 1; /* success */
}

/* returns oriented edge number, based at 1, not 0 */
int old_find_edge(a,b,wrapa,wrapb)
int a,b;  /* vertex numbers of neigbors */
int *wrapa,*wrapb;  /* wraps from seed, to figure edge wrap */
{ int k;
  struct edge *e;

  /* linear search; only have to test opposite direction
     of existing edge */
  for ( k = 0, e = edges ; k < ecount ; k++, e++ )
   if ( (a==e->v[1]) && (b==e->v[0]) && (wrapa[0]-wrapb[0]==e->wrap[0])
         && (wrapa[1]-wrapb[1]==e->wrap[1]) ) return -(k+1);

  /* not found, so add */
  e->v[0] = a; e->v[1] = b;
  /* set wrap */
  for ( k = 0 ; k < 2  ; k++ )  /* unwrap second vertex */
   e->wrap[k] = wrapb[k]-wrapa[k];
 

  ecount++;
  if ( ecount > 3*N ) {fputs(".",stderr); return 0; }
  return ecount;
}

/* sees if edge already exists, adds if not */
/* returns oriented edge number, based at 1, not 0 */
int find_edge(a,b)
int a,b;  /* vertex numbers of neigbors */
{ int k;
  struct edge *e;
  struct vertex *v, *vv;

  /* search around first vertex for existing edge */
  v = verts + a;
  vv = verts + b;
  for ( k = 0 ; k < 3 ; k++ )
    if ( v->nbrv[k] == vv )
      break;
  if ( k == 3 )
  { fprintf(stderr,"Cannot find neighbor vertex in find_edge().\n");
    return 0;
  }
  if ( v->nbre[k] )
  { /* already existing, so check orientation */
    if ( v->nbre[k]->v[0] == a ) return v->nbre[k] - edges + 1;
    else return -(v->nbre[k] - edges + 1);
  }

  /* not found, so add */
  e = edges + ecount;
  e->v[0] = a; e->v[1] = b;
  v->nbre[k] = e;
  for ( k = 0 ; k < 3 ; k++ )
    if ( vv->nbrv[k] == v )
    { vv->nbre[k] = e; break; }

  /* set wrap */
  for ( k = 0 ; k < 2  ; k++ )  /* unwrap second vertex */
    { if (verts[a].x[k] -verts[b].x[k] > 0.5 ) e->wrap[k] = 1;
      else if (verts[a].x[k] -verts[b].x[k] < -0.5 ) e->wrap[k] = -1;
      else e->wrap[k] = 0;
    }

  ecount++;
  if ( ecount > edgealloc ) {fputs(".",stderr); return 0; }
  return ecount;
}

/* assemble_net() takes vertex-cell data from calc_voronoi() and
   creates edges and sets up vertex-vertex neighbor
   relations.
*/
int assemble_net()
{ int i,j,k;
  struct cell *c;
 
  /* do edges */
  for ( i = 0, c = cells ; i < ccount ; i++,c++ )
   { c->v[c->count] = c->v[0];  /* wrap around list */
     for ( j = 0 ; j < 2 ; j++ )
       c->vwraps[c->count][j] = c->vwraps[0][j];
     for ( k = 0 ; k < c->count ; k++ )
     { int e = old_find_edge(c->v[k],c->v[k+1],c->vwraps[k],c->vwraps[k+1]);
       c->e[k] = e;
     }
   }

  /* assign nbr vertices to vertices */
  
  for ( i = 0 ; i < vcount ; i++ )
    for ( j = 0 ; j < 3 ; j++ )
      verts[i].nbrv[j] = NULL;
      
  for ( i = 0 ; i < ecount ; i++ )
  { struct vertex *v;
  
    v = verts + edges[i].v[0];
    for ( j = 0 ; j < 3 ; j++ )
      if ( v->nbrv[j] == NULL ) 
        { v->nbrv[j] = verts + edges[i].v[1]; 
          v->nwrap[j][0] = edges[i].wrap[0]; 
          v->nwrap[j][1] = edges[i].wrap[1];
          break; 
        }
    if ( j == 3 ) 
      printf("too many neighbors.\n");
    
    v = verts + edges[i].v[1];
    for ( j = 0 ; j < 3 ; j++ )
      if ( v->nbrv[j] == NULL ) 
        { v->nbrv[j] = verts + edges[i].v[0]; 
          v->nwrap[j][0] = -edges[i].wrap[0]; 
          v->nwrap[j][1] = -edges[i].wrap[1];
          break; 
        }
    if ( j == 3 ) 
      printf("too many neighbors.\n");
  }
  /* get in counterclockwise order */
  for ( i = 0 ; i < vcount ; i++ )
  { struct vertex *v  = verts+i;
    struct vertex *v0 = v->nbrv[0];
    struct vertex *v1 = v->nbrv[1];
    double x0 = v0->x[0]-v->x[0]+v->nwrap[0][0];
    double y0 = v0->x[1]-v->x[1]+v->nwrap[0][1];
    double x1 = v1->x[0]-v->x[0]+v->nwrap[1][0];
    double y1 = v1->x[1]-v->x[1]+v->nwrap[1][1];
    
    if ( x0*y1 - y0*x1 < 0.0 )
    { struct vertex *temp = v->nbrv[0];
      v->nbrv[0] = v->nbrv[1];
      v->nbrv[1] = temp;
      swap(v->nwrap[0][0],v->nwrap[1][0]);
      swap(v->nwrap[0][1],v->nwrap[1][1]);
    }
  }

  /* get cell links in counterclockwise order, cell[i] before nbrv[i] */
  for ( i = 0 ; i < vcount ; i++ )
  { struct vertex *v = verts+i;
    int n;
    for ( n = 0 ; n < 3; n++ )
    {
      struct vertex *va = v->nbrv[n];
      struct vertex *vb = v->nbrv[(n+1)%3];
      for ( j = 0 ; j < 3 ; j++ )
        for ( k = 0 ; k < 3 ; k++ )
          if ( va->cell[j] == vb->cell[k] )
          { int m;
            for ( m = 0 ; m < 3 ; m++ )
            { if ( v->cell[m] == va->cell[j] )
              { int tmp = v->cell[(n+1)%3];
                v->cell[(n+1)%3] = v->cell[m];
                v->cell[m] = tmp;
                goto ccexit;
              }
            }
          }
ccexit: ;
        }
  }

  /* calculate triangle areas and initialize area tree */
  for ( i = 0 ; i < vcount ; i++ )
  { struct vertex *v = verts + i;
    double x0,x1,y0,y1; /* unwrapped triangle sides */
    struct cell *c0 = cells + v->cell[0];
    struct cell *c1 = cells + v->cell[1];
    struct cell *c2 = cells + v->cell[2];
    x0 = (c1->x[0]+v->cwraps[1][0]) - (c0->x[0]+v->cwraps[0][0]);
    x1 = (c2->x[0]+v->cwraps[2][0]) - (c0->x[0]+v->cwraps[0][0]);
    y0 = (c1->x[1]+v->cwraps[1][1]) - (c0->x[1]+v->cwraps[0][1]);
    y1 = (c2->x[1]+v->cwraps[2][1]) - (c0->x[1]+v->cwraps[0][1]);
    
    v->area = (x0*y1 - x1*y0)/2;
    add_triangle(v);
  }

  return 1;
}

void check_orders()
{ int i,j,k;
  /* get cell links in counterclockwise order, cell[i] before nbrv[i] */
  for ( i = 0 ; i < 2*N+10 ; i++ )
  { struct vertex *v = verts+i;
    int n;
    if (v->flag == VFREE ) continue;
    for ( n = 0 ; n < 3; n++ )
    {
      struct vertex *va = v->nbrv[n];
      struct vertex *vb = v->nbrv[(n+1)%3];
      for ( j = 0 ; j < 3 ; j++ )
      { for ( k = 0 ; k < 3 ; k++ )
          if ( va->cell[j] == vb->cell[k] )
          { if ( v->cell[(n+1)%3] == va->cell[j] ) goto ccc;                  
          }
      }
                   
      if (  j == 3 ) 
      { printf("Can't find adjacent cell.\n");
      }
ccc: ;
    }
  }
}

/* Distance squared between points, accounting for wraps */
double dist(double *a, double *b)
{ double p = a[0] - b[0];
  double q = a[1] - b[1];
  double dd = p*p + q*q;
  if ( dd < 0.2 ) return dd;
  if ( p > .5 ) p -= 1.0;
  else if ( p < -.5 ) p += 1.0;
  if ( q > .5 ) q -= 1.0;
  else if ( q < -.5 ) q += 1.0;
  return p*p + q*q;
}


/* re-make edge list in its entirety */
/* to be called after inserting all cells */
void redo_edges()
{ int i,k;
  struct cell *c;

  /* do edges again */
  ecount = 0;
  for ( i = 0, c = cells ; i < ccount ; i++,c++ )
   { c->v[c->count] = c->v[0];  /* wrap around list */
     for ( k = 0 ; k < c->count ; k++ )
     { int e = find_edge(c->v[k],c->v[k+1]);
       c->e[k] = e;
     }
   }
}

/* Add a vertex */
void addpoint(void)
{ 
  struct cell *c = cells+ccount;
  double x = c->x[0],y = c->x[1];
  struct vertex *minv;
  int i;

  double a,b; /* barycentric coordinates */

  minv = find_next();
  if ( minv->flag == VFREE )
    printf ("Bad minv.\n");

  /* get random barycentric coordinates */
  a = frand();
  b = frand();
  if ( a + b > 1.0 ) { a = 1.0 - a; b = 1.0 - b; }

  { double x0,x1,y0,y1; /* unwrapped triangle sides */
    struct cell *c0 = cells + minv->cell[0];
    struct cell *c1 = cells + minv->cell[1];
    struct cell *c2 = cells + minv->cell[2];
    x0 = c1->x[0] - c0->x[0];
    x1 = c2->x[0] - c0->x[0];
    y0 = c1->x[1] - c0->x[1];
    y1 = c2->x[1] - c0->x[1];
    if ( x0 < -.5 ) x0 += 1.0; else if ( x0 > .5 ) x0 -= 1.0;
    if ( y0 < -.5 ) y0 += 1.0; else if ( y0 > .5 ) y0 -= 1.0;
    if ( x1 < -.5 ) x1 += 1.0; else if ( x1 > .5 ) x1 -= 1.0;
    if ( y1 < -.5 ) y1 += 1.0; else if ( y1 > .5 ) y1 -= 1.0;
    x = c0->x[0] + a*x0 + b*x1;
    y = c0->x[1] + a*y0 + b*y1;
    if ( x < 0.0 ) x += 1.0; else if ( x > 1.0 ) x -= 1.0;
    if ( y < 0.0 ) y += 1.0; else if ( y > 1.0 ) y -= 1.0;
    c->x[0] = x; c->x[1] = y;
#ifdef VORDEBUG
    printf("New seed %d at %f %f\n",ccount+1,c->x[0],c->x[1]);
#endif
  }


  /* now find tree of deletable vertices */

  newvcount = oldvcount = 0;
  findtree(minv,NULL,c);
  ccount++;
  c->count = newvcount;
  for ( i = 0 ; i < newvcount ; i++ )
  { c->v[newvcount-i-1] = newvlist[i] - verts;
    c->vv[newvcount-i-1] = newvlist[i];
  }

  /* set neighbors on new vertices */
  for ( i = 0 ; i < newvcount ; i++ )
  { newvlist[i]->nbrv[1] = i == newvcount-1 ? newvlist[0] : newvlist[i+1];
    newvlist[i]->nbrv[2] = i == 0 ? newvlist[newvcount-1] : newvlist[i-1];
  }

  /* remove old vertices from cell lists */
  for ( i = 0 ; i < oldvcount ; i++ )
  { int n,k;
    for ( n = 0 ; n < 3 ; n++ )
    { struct cell *cc = cells+oldvlist[i]->cell[n];
      for ( k = 0 ; k < cc->count ; k++ )
        if ( cc->vv[k] == oldvlist[i] )
          break;
      for ( ; k < cc->count-1 ; k++ )
        { cc->vv[k] = cc->vv[k+1];
          cc->v[k] = cc->v[k+1];
        }
      cc->count--;
    }
    oldvlist[i]->area = 0.0;
    del_triangle(oldvlist[i]);
    oldvlist[i]->nextfree = freehead;
    freehead = oldvlist[i];
    freehead->flag &= ~VUSED;
  }
}

/* One recursive stage in finding deletable tree */
/* New vertices put in cell list in order, but old ones not deleted yet */
void findtree(struct vertex *v, struct vertex *parent,struct cell *c)
{ struct vertex *vv;
  int k,j,i,m,mm;
#ifdef VORDEBUG
printf("Deletable vertex %d at %f %f\n",v-verts,v->x[0],v->x[1]);
#endif
   v->flag = VFREE;     
   oldvlist[oldvcount++] = v;

  /* go counterclockwise from parent */
  for ( k = 0 ; k < 3 ; k++ )
   if ( v->nbrv[k] == parent ) break;

  for ( j = 0, k++ ; j < 3 ; j++, k++ )  // keep on going round
  { vv = v->nbrv[k%3];
    if ( vv == parent ) break;
        
    if ( (vv->flag == VUSED) && (dist(c->x,vv->x) < vv->rr) )
    {
      /* recurse */
      findtree(vv,v,c);
    }
    else // have edge with new vertex on it
    { struct cell *c1,*c2;
      struct vertex *newv;
      double x,y,lambda,dx,dy;
      double c1x,c1y,c2x,c2y;
#ifdef VORDEBUG
printf("Nbr vertex of %d at %f %f\n",v-verts,vv->x[0],vv->x[1]);
#endif
      /* find coordinates of new vertex */
      c1 = cells+v->cell[k%3];
      c2 = cells+v->cell[(k+1)%3];

      c1x = c1->x[0] - c->x[0]; 
      if ( c1x > .5 ) c1x--;
      else if ( c1x < -.5 ) c1x++;

      c1y = c1->x[1] - c->x[1]; 
      if ( c1y > .5 ) c1y--;
      else if ( c1y< -.5 ) c1y++;

      c2x = c2->x[0] - c->x[0]; if ( c2x > .5 ) c2x--;
      else if ( c2x < -.5 ) c2x++;

      c2y = c2->x[1] - c->x[1]; if ( c2y > .5 ) c2y--;
      else if ( c2y < -.5 ) c2y++;

      lambda = ((c2x - c1x)*( - c2x) - (c2y - c1y)*(c2y ))/
                            (( c1y)*( - c2x) - (-c1x )*(c2y ));
      newv = freehead; 
      freehead = newv->nextfree; 
      newv->nextfree = NULL;
      newv->flag = VUSED;

      newvlist[newvcount++] = newv;

      dx = (c1x + lambda*(c1y))/2;
      dy = (c1y + lambda*(-c1x))/2;
      newv->rr = dx*dx + dy*dy;
      x = c->x[0] + dx;
      y = c->x[1] + dy;
      if ( x < 0.0 ) x += 1.0; else if ( x > 1.0 ) x -= 1.0;
      if ( y < 0.0 ) y += 1.0; else if ( y > 1.0 ) y -= 1.0;
      newv->x[0] = x; newv->x[1] = y;
      newv->nbrv[0] = vv;
      newv->area = (c1x*c2y - c1y*c2x)/2;
      add_triangle(newv);

      newv->cell[0] = v->cell[k%3];
      newv->cell[1] = v->cell[(k+1)%3];
      newv->cell[2] = ccount;
          
#ifdef VORDEBUG
      if ( newv->cell[0] == newv->cell[1] ||
        newv->cell[1] == newv->cell[2] ||
        newv->cell[2] == newv->cell[1] )
        printf("common cell!!\n");
#endif

      /* put new vertex in cells' vertex lists */
      for ( m = 0 ; m < c1->count-1 ; m++ )
      { if ( ((c1->vv[m] == v) && (c1->vv[m+1] == vv) )
           || ((c1->vv[m] == vv) && (c1->vv[m+1] == v) ) )
        break;
      }
      for ( mm = c1->count-1 ; mm > m ; mm-- )
      { c1->vv[mm+1] = c1->vv[mm];
        c1->v[mm+1] = c1->v[mm];
      }
      c1->vv[m+1] = newv;
      c1->v[m+1] = newv-verts;
      c1->count++;

      for ( m = 0 ; m < c2->count-1 ; m++ )
      { if ( ((c2->vv[m] == v) && (c2->vv[m+1] == vv) )
           || ((c2->vv[m] == vv) && (c2->vv[m+1] == v) ) )
        break;
      }      
      for ( mm = c2->count-1 ; mm > m ; mm-- )
      { c2->vv[mm+1] = c2->vv[mm];
        c2->v[mm+1] = c2->v[mm];
      }
      c2->vv[m+1] = newv;
      c2->v[m+1] = newv-verts;
      c2->count++;

      // get vertex into fundamental domain
      if ( newv->x[0] < 0.0 ) newv->x[0]++;
      else if ( newv->x[0] > 1.0 ) newv->x[0]--;
      if ( newv->x[1] < 0.0 ) newv->x[1]++;
      else if ( newv->x[1] > 1.0 ) newv->x[1]--;
      // reset neighbor vertex neighbor to this
      for ( i = 0 ; i < 3 ; i++ )
        if ( vv->nbrv[i] == v )
           { vv->nbrv[i] = newv; break; }
    }
  }

}


void do_datafile()
{ int k,j;
  struct vertex *v;
  struct edge *e;
  struct cell *c;
  int vcounter;

  printf("// 2D Voronoi diagram of %d cells in unit torus.\n\n",N);
  puts("STRING\nSPACE_DIMENSION 2\n");
  puts("TORUS_FILLED\nPERIODS\n1.0 0.0\n0.0 1.0\n");
  if ( phaseflag )
    puts("PHASEFILE \"vorphase.dat\"\n");

  if ( foam_flag )
  { puts("\nlength_method_name \"circular_arc_length\"\n");
    puts("\narea_method_name \"circular_arc_area\"\n");
  }
  else
  { /* for automatic evolution */
    puts("autopop");
    printf("autochop %f\n",0.15/sqrt((double)N));
    puts("area_normalization");
/*  puts("effective_area"); */
    printf("scale %g fixed\n",5e-4/N);
    puts("");
  }
 
  if ( multiflag )
  { printf("define vertex attribute vpart real\n");
    printf("define edge   attribute epart real\n");
    printf("define facet  attribute fpart real\n\n");
  }
  printf("\ndefine facet attribute center real[2] // seed coordinates \n");

  puts("VERTICES");
  for ( k = 0, v = verts, vcounter = 0 ; k < 2*N+100 ; k++,v++ )
  { if ( v->flag & VUSED )
    { vcounter++;
      if ( multiflag )
        printf("%d@%d   %f  %f vpart %d\n",vcounter,v->proc,v->x[0],v->x[1],
             v->proc);
      else
        printf("%d   %f  %f\n",vcounter,v->x[0],v->x[1]);
      v->num = vcounter;
    }
  }
  puts("\nEDGES");
  for ( k = 0, e = edges ; k < ecount ; k++,e++ )
  { char wx = (e->wrap[0] == 0) ? '*' : (e->wrap[0]==1 ? '+' : '-');
    char wy = (e->wrap[1] == 0) ? '*' : (e->wrap[1]==1 ? '+' : '-');
    if ( multiflag )
      printf("%d@%d   %d@%d  %d@%d   %c %c epart %d\n",k+1,e->proc,
       verts[e->v[0]].num,verts[e->v[0]].proc,
       verts[e->v[1]].num,verts[e->v[1]].proc, wx,wy,
       e->proc);
    else
      printf("%d   %d  %d   %c %c\n",k+1,verts[e->v[0]].num,
       verts[e->v[1]].num, wx,wy);
  }
  puts("\nFACES");
  for ( k = 0, c = cells ; k < N ; k++,c++ )
    { if ( multiflag )
      { printf("%d@%d  ",k+1,c->proc);
        for ( j = c->count - 1 ; j >= 0;  j-- )
          printf("%d@%d ",-c->e[j],edges[abs(c->e[j])-1].proc); 
        printf(" fpart %d\n",c->proc);
      }
      else
      { printf("%d  ",k+1);
        for ( j = c->count - 1 ; j >= 0;  j-- )
          printf("%d ",-c->e[j]);  /* getting positive area */
      }
      printf(" center {%18.15f,%18.15f} ",cells[k].x[0],cells[k].x[1]);
      if ( phaseflag )
          printf(" phase %d ",((rand()+rand()/1000)%phases)+1);
      printf("\n");
    }
  puts("\nBODIES");
  for ( k = 0 ; k < N ; k++ )
  { if ( multiflag )
      printf("%d@%d  %d@%d\n",k+1,cells[k].proc,k+1,cells[k].proc);
    else
      printf("%d  %d\n",k+1,k+1);
  }

  puts("\nread");
  puts("");
  puts("clipped"); /* for clipped torus mode display */
  /*puts("// set volconsts so body areas correct");
  puts("set body volconst -floor(volume+0.5)");
  if ( foam_flag )
  {
    puts("set body target volume");
    puts("");
    puts("set_pressures := {");
    puts("    ph := 1;");
    puts("    while ( ph <= 5 ) do");
    puts("    { surften := random;");
    puts("      set body bb pressure -surften where bb.facet[1].phase == ph;");
    puts("      ph += 1;");
    puts("    }");
    puts("}");
  }
  else
  { puts("// re-define 'g' so autochop and scale factor change with number of cells"); 
    printf("g :::= { autochop 0.15/sqrt(facet_count); \n");
    printf("         scale := 5e-4/facet_count;\n");
    printf("         'g'; \n");
    printf("       }\n");
  }*/
}




/**************************************************************************

  Random number generation functions adapted from comprand.c 
  ("MathLink Program for High-Quality Random Number" by
    A. Compagner, Delft U. of Technology and 
    A. S. Berdnikov, S. B. Turtia, Institute dor Analytical
    Instrumentation, St. Petersburg, Russia

  Modifications: File variables declared static; 
                 Removed timer init since nonportable.
                 Some functions renamed, others removed.
                 kb_drand() extended to full random bits

****************************************************************************/

/*
   The following is the list of primitive trinomials
   for Mersenne primes. The list is composed from 
   data taken from:
        N.Zierler (up to p=9689), 
               in: Inform.Control 15(1969)67 
        J.R.Heringa, H.W.J.Bloete, A.Compagner,
               in: Int.J. of Mod.Phys. C3(1992)561 

   The first number in the list is the number of the
   Mersenne prime, the second number is the Mersenne 
   prime p itself, and the last number is the degree 
   of the third term of the characteristic trinomial, q. 
   Sometimes there is no such trinomial (this is listed 
   as NONE), sometimes there are different possibilities 
   for q (they will be seperated by comma's). 

       #           p             q

       1           2             1 
       2           3             1
       3           5             2
       4           7             1, 3
       5          13             NONE
       6          17             3, 5, 6
       7          19             NONE
       8          31             3, 6, 7, 13
       9          61             NONE
      10          89             38
      11         107             NONE
      12         127             1, 7, 15, 30, 63
      13         521             32, 48, 158, 168
      14         607             105, 147, 273
      15        1279             216, 418
      16        2203             NONE
      17        2281             715, 915, 1029 
      18        3217             67, 576
      19        4253             NONE
      20        4423             271, 369, 370, 649, 1393, 1419, 2098
      21        9689             84, 471, 1836, 2444, 4187
      22        9941             NONE
      23       11213             NONE
      24       19937             881, 7083, 9842
      25       21701             NONE
      26       23209             1530, 6619, 9739
      27       44497             8575, 21034
      28       86243             NONE
      29      110503             25230, 53719
      30      132049             7000, 33912, 41469, 52549, 54454
      31      216091             NONE

   Of course, also the complements p-q give rise to
   primitive trinomials, when they are used instead of q. 

   Note that it does not help to multiply m trinomials with 
   different p if 3**m is larger than one half of the sum of
   the values of p of the factors. 

   For efficiency, combining 4 trinomials is very likely 
   to give reliable Monte Carlo results up to the end of 
   time, and probably 3 would alreday be good enough for 
   that. 

   To specify constants below the recepies
           #10, #12, #21 and #27
   are selected which gives the period of 
   pseudorandom sequence equal to: 
        (2^89-1)*(2^127-1)*(2^9689-1)*(2^21034-1)
   which is enough for any reasonable application.
 
   In principle the user can vary the number of recepies 
   and the values (p,q) used in intiation of global arrays.
   All "p" used in recepies are different, and it is better
   (but not necessary) if q ~ p/2. 
   Note that since the length of the word is equal to 32,
   it is essential that q > 32 for all recepies!
*/

/*
   Procedures use explicitly 
   that data types:
           unsigned long int, 
           long int, 
           int
   contain 32 bits
 */

#define wordsize 32
#define intmask  0x7FFFFFFF
#define IREGmask 0xFFFFFFFF

/* 
   define constants for shift registers taking
   into account that unsigned long int contains 32 bits

      numreg  - the total number of shift registers
      ndim.i = pval[i]/32 + 1
      pval[i], qval[i] - primitive trinomial constants
*/

#define numreg 4

#define dim0    3 
#define dim1    4
#define dim2  303
#define dim3 1391 

static int Pval[numreg] = {89, 127, 9689, 44497};
static int Qval[numreg] = {38,  63, 4187, 21034};

/* global internal variables */ 

static double rfactor = 0.0;
static int    notinit = 1;

typedef unsigned long int IREG;

static IREG Reg0[dim0], Reg1[dim1], Reg2[dim2], Reg3[dim3];

static IREG * Regs[numreg] = {Reg0, Reg1, Reg2, Reg3};

static int Index[numreg],
    NshR [numreg],  NshL[numreg],
    KshR [numreg],  KshL[numreg],
    Nbase[numreg], Kbase[numreg];

/* Initialization procedure */
int kb_initr(ind)
    int ind;
{
    int i, k;
    unsigned long int IJKL;
    unsigned long int m1, m2, m3, m4;
    IREG mm;

    /* 
       calculate scaling factor to convert 
       unsigned long int to double 
    */
    rfactor=1.0/(IREGmask+1.0);

    /*
        calculate index constants
    */
    for (k=0; k < numreg; k++)
    {
       Nbase[k]=Pval[k]/wordsize;
       NshL[k]= Pval[k] - Nbase[k]*wordsize;
       NshR[k]=wordsize - NshL[k];       

       Kbase[k]=Qval[k]/wordsize;
       KshL[k]= Qval[k]  - Kbase[k]*wordsize;
       KshR[k]= wordsize - KshL[k]; 

       Index[k]=0;      
    }

    /*
       initiate shift registers
    */
    IJKL=ind;
    if (! IJKL) IJKL=1;
    for (k=0; k < numreg; k++)
    {
        for (i=0; i <= Nbase[k]; i++)
        {
           IJKL=(1726*IJKL) % 8191; m1=IJKL;
           IJKL=(1726*IJKL) % 8191; m2=IJKL;
           IJKL=(1726*IJKL) % 8191; m3=IJKL;
           IJKL=(1726*IJKL) % 8191; m4=IJKL;
           mm=(m1 >> 5) | ((m2 >> 5) << 8)
                        | ((m3 >> 5) << 16)
                        | ((m4 >> 5) << 24);
           if (! mm) mm=1; 
           (Regs[k])[i]=mm; 
        }
    }

    notinit=0;

    IJKL=0;
    for (k=0; k < numreg; k++)
    {
        IJKL=IJKL + Pval[k];
    }

    IJKL=IJKL/wordsize;
    return IJKL;
}

static IREG IREGrand (void);
static IREG IREGrand()
{
   int k, ind, i1, i2, j1, j2;
   IREG base, baseH;

   if (notinit) kb_initr(1);

   base=0;
   for(k=0; k < numreg; k++)
   {
      ind=Index[k];
      if (ind) {j2=ind--;} else {j2=0; ind=Nbase[k];};
      Index[k]=ind; j1=ind;

      i1=ind + Kbase[k];
      if (i1 > Nbase[k]) i1=i1-Nbase[k];
      i2=i1+1;
      if (i2 > Nbase[k]) i2=0;

      baseH=IREGmask & (
            (((Regs[k])[i1] << KshL[k]) | ((Regs[k])[i2] >> KshR[k])) 
            ^ 
            (((Regs[k])[j1] << NshL[k]) | ((Regs[k])[j2] >> NshR[k]))
            );

      (Regs[k])[ind]=baseH;
       base =base  ^ baseH; 
   } 

   return base;
}


double kb_drand()
{
   return (IREGrand()*rfactor + IREGrand())*rfactor;
}


