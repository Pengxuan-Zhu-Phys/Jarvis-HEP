/*
	parallel.c
		parallel execution of SquaredMEHel
		this file is part of FormCalc
		last modified 6 Sep 12 th
*/


#define HAVE_FORK
//#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <sys/select.h>

#if NOUNDERSCORE
#define sqmeprep_ sqmeprep
#define sqmeexec_ sqmeexec
#define sqmesync_ sqmesync
#define sqmewait_ sqmewait
#define restorecache_ restorecache
#define clearcache_ clearcache
#define cubasetexit_ cubasetexit
#endif

#define RealType double
#define NCOMP 2

typedef void (*subroutine)(RealType *, const int *);

#ifdef HAVE_FORK

typedef unsigned long long int seq_t;

typedef struct {
  int fd;
  seq_t seq;
} childinfo;

#define MINCORES 1
#define SEQ_S(seq) ((seq) & 0xffffffff)
#define SEQ_ANGLE(seq) ((seq) & 0xffffffff00000000)

#ifdef DEBUG
#define TERM_RED "\e[31m"
#define TERM_BLUE "\e[34m"
#define TERM_RESET "\e[0m\n"
#define MASTER(s, ...) \
fprintf(stderr, TERM_RED "SquaredME master: " s TERM_RESET, __VA_ARGS__)
#define WORKER(core, s, ...) \
fprintf(stderr, TERM_BLUE "SquaredME core %d: " s TERM_RESET, core, __VA_ARGS__)
#else
#define MASTER(s, ...)
#define WORKER(core, s, ...)
#endif

/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))

/* (a < 0) ? 0 : a */
#define IDim(a) ((a) & NegQ(-(a)))

/* (a < b) ? a : b */
#define IMin(a, b) ((a) - IDim((a) - (b)))

/* (a > b) ? a : b */
#define IMax(a, b) ((b) + IDim((a) - (b)))

static int ncores = -1, nlaunched, nrunning;
static int fdmax;
static childinfo *child;
static fd_set children;
static struct {
  char *h, *he;
  char *v, *ve;
  char *a, *ae;
  char *s, *se;
} mem;

#define mem_seq *(seq_t *)mem.h

#ifdef HAVE_GETLOADAVG
static double loadavg = -1;
#endif

/*********************************************************************/

#ifndef MSG_WAITALL
/* Windows */
#define MSG_WAITALL 0
#endif

static inline int readsock(int fd, void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = recv(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

static inline int writesock(int fd, const void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = send(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

/*********************************************************************/

static inline void newcore(subroutine foo, const int flags)
{
  int fd[2];
  pid_t pid;
  int core = nlaunched++;

  assert(
    socketpair(AF_LOCAL, SOCK_STREAM, 0, fd) != -1 &&
    (pid = fork()) != -1 );

  if( pid ) {
    MASTER("forked core %d pid %d pipe %d(master) -> %d(worker) seq %llx",
      core, pid, fd[0], fd[1], mem_seq);
    close(fd[1]);
    child[core].fd = fd[0];
    child[core].seq = mem_seq;
    FD_SET(fd[0], &children);
    fdmax = IMax(fd[0], fdmax);
    return;
  }

  close(fd[0]);

  for( ; ; ) {
    RealType res[NCOMP];
    seq_t seq = mem_seq;

    res[0] = res[1] = 0;
    foo(res, &flags);
    WORKER(core, "writing result(%ld)", sizeof res);
    writesock(fd[1], res, sizeof res);

    WORKER(core, "reading mem_hel(%p#%ld)", mem.h, mem.he - mem.h);
    if( !readsock(fd[1], mem.h, mem.he - mem.h) ) exit(0);
    WORKER(core, "seq %llx  new %llx", seq, mem_seq);
    seq ^= mem_seq;
    if( SEQ_ANGLE(seq) ) {
      WORKER(core, "reading mem_angle(%p#%ld+%p#%ld)",
        mem.v, mem.ve - mem.v, mem.a, mem.ae - mem.a);
      readsock(fd[1], mem.v, mem.ve - mem.v);
      readsock(fd[1], mem.a, mem.ae - mem.a);
      restorecache_();
    }
    if( SEQ_S(seq) ) {
      WORKER(core, "reading mem_s(%p#%ld)", mem.s, mem.se - mem.s);
      readsock(fd[1], mem.s, mem.se - mem.s);
      clearcache_();
    }
    seq = mem_seq;
  }
}

/*********************************************************************/

static inline int readycore(RealType *result)
{
  static int nwaiting, core;
  static fd_set ready;
  RealType res[NCOMP];
  int c;

  if( nwaiting == 0 ) {
    memcpy(&ready, &children, sizeof ready);
    nwaiting = select(fdmax + 1, &ready, NULL, NULL, NULL);
    core = 0;
  }
  --nwaiting;

  for( ; core < ncores; ++core )
    if( FD_ISSET(child[core].fd, &ready) ) break;

  MASTER("reading result(%ld) from core %d", sizeof res, core);

  readsock(child[core].fd, res, sizeof res);
  for( c = 0; c < NCOMP; ++c ) result[c] += res[c];

  --nrunning;
  return core++;
}

/*********************************************************************/

static inline void oldcore(const int core)
{
  int fd = child[core].fd;
  seq_t seq = child[core].seq;
  child[core].seq = mem_seq;

  MASTER("sending mem_hel(%p#%ld) to core %d",
    mem.h, mem.he - mem.h, core);
  writesock(fd, mem.h, mem.he - mem.h);
  seq ^= mem_seq;
  if( SEQ_ANGLE(seq) ) {
    MASTER("sending mem_angle(%p#%ld+%p#%ld) to core %d",
      mem.v, mem.ve - mem.v, mem.a, mem.ae - mem.a, core);
    writesock(fd, mem.v, mem.ve - mem.v);
    writesock(fd, mem.a, mem.ae - mem.a);
  }
  if( SEQ_S(seq) ) {
    MASTER("sending mem_s(%p#%ld) to core %d",
      mem.s, mem.se - mem.s, core);
    writesock(fd, mem.s, mem.se - mem.s);
  }
}

#endif

/*********************************************************************/

void sqmeprep_(char *h, char *he, char *v, char *ve,
  char *a, char *ae, char *s, char *se)
{
#ifdef HAVE_FORK
  void sqmewait_(void);
  void cubasetexit_(void (*)(void), void *);

  mem.h = h;
  mem.he = he;
  mem.v = v;
  mem.ve = ve;
  mem.a = a;
  mem.ae = ae;
  mem.s = s;
  mem.se = se;

  if( ncores < 0 ) {
    const char *env = getenv("FCCORES");
    ncores = env ? atoi(env) : sysconf(_SC_NPROCESSORS_ONLN);
#ifdef HAVE_GETLOADAVG
    if( env == NULL || ncores < 0 ) {
      if( loadavg < 0 ) getloadavg(&loadavg, 1);
      ncores = abs(ncores) - floor(loadavg);
    }
#endif
    printf("SquaredME: using %d cores\n", ncores);
    fflush(NULL);
    assert( child = malloc(ncores*sizeof *child) );
    cubasetexit_(sqmewait_, NULL);
  }
#endif
}

/*********************************************************************/

void sqmeexec_(subroutine foo, RealType *result, const int *flags)
{
#ifdef HAVE_FORK
  if( ncores >= MINCORES ) {
    if( nlaunched < ncores ) newcore(foo, *flags);
    else oldcore((nrunning < ncores) ? nrunning : readycore(result));
    ++nrunning;
  }
  else
#endif
  foo(result, flags);
}

/*********************************************************************/

void sqmesync_(RealType *result)
{
#ifdef HAVE_FORK
  while( nrunning ) readycore(result);
#endif
}

/*********************************************************************/

void sqmewait_()
{
#ifdef HAVE_FORK
  if( nlaunched ) {
    int core;
    pid_t pid;
    for( core = 0; core < nlaunched; ++core ) {
      MASTER("closing core %d", core);
      close(child[core].fd);
    }
    for( core = 0; core < ncores; ++core ) {
      MASTER("waiting for core %d", core);
      wait(&pid);
      MASTER("core %d pid %d terminated", core, pid);
    }
    free(child);
    fdmax = 0;
    nlaunched = 0;
  }
#endif
}

