/*JS*********************************************************************
*
*    File    : GENSORT.H
*    Language: ANSI-C
*    Author  : Joerg Schoen
*    Purpose : To generate automatically code for a sorting routine.
*
*              $Id: gensort.h,v 1.1 2004/07/28 20:17:57 fregeau Exp $
*
*  You have to define the following preprocessor symbols:
*  GENSORT_NAME       will be the name of the generated routine.
*  GENSORT_TYPE       is the type of the array to be sorted.
*  GENSORT_KEYTYPE    is the type of the key of the array elements.
*  GENSORT_GETKEY(a)  is a macro to get the key element of the array
*                     element a.
*  GENSORT_COMPAREKEYS(k1,k2) compares the two key elements k1 and k2 with
*                     condition "lower than" and yields TRUE / FALSE values.
*                     If the two keys match, it must return FALSE!
*
*  With the line
*     #include "gensort.h"
*  this file will automatically generate a sort routine.
*
*  The generated routine will have the prototype
*     void GENSORT_NAME(GENSORT_TYPE *array,long n)
*   and sorts the n array entries array[0] to array[n - 1].
*
*  At the end of this file the symbols GENSORT_* are undefined for reuse.
*
*  Flags:
*   If the flag "GENSORT_USEPOINTERS" is defined, the routine uses pointer
*    gymnastics instead of indices, which is faster on some machines.
*
*   If the flag "GENSORT_NISINT" is defined, the size argument n to the
*    generated routine is an integer instead of the default long.
*
*   If the flag "GENSORT_NOSMALL" is defined, the quicksort routines do
*    not use a straight insertion algorithm for small n. The threshold
*    for n can be changed by defining "GENSORT_LIMIT" as the threshold,
*    the default is given below. Straight insertion is not used if the
*    macro "GENSORT_SWAP" is defined.
*
*   If the flag "GENSORT_HEAPSORT" is defined, the sort routine uses heap-
*    sort instead of the default quicksort algorithm.
*
*   If the symbol "GENSORT_ARGSPROTO" and "GENSORT_ARGS" is defined, the
*    function has additional arguments.
*
*   If the symbol "GENSORT_SWAP" is defined, exchange of two elements is
*    done with the macro GENSORT_SWAP instead with a temporary variable.
*    The macro is called via "GENSORT_SWAP(a,i,j);" enclosed in
*    curly brackets and is expected to exchange the elements with index i
*    and j in the array a. Usage of "GENSORT_SWAP" together with
*    "GENSORT_HEAPSORT" is not possible!
*
*  The quicksort is a recursive routine which uses at most
*      log_2(n) * (3 * sizeof(GENSORT_TYPE*) + sizeof(GENSORT_KEYTYPE)
*                  sizeof(GENSORT_TYPE))
*   bytes of stack (may be smaller if optimizer is good).
*
*  Example 1: To generate a quicksort routine which sorts an integer array
*   in increasing order, use the following lines:
*      | Left margin
*      V
*      #define GENSORT_NAME                 intquicksort
*      #define GENSORT_TYPE                 int
*      #define GENSORT_KEYTYPE              int
*      #define GENSORT_GETKEY(a)            a
*      #define GENSORT_COMPAREKEYS(k1,k2)   k1 < k2
*
*      #include "gensort.h"
*
*  Example 2: To generate a pointer quicksort routine that sorts a structure
*   according to a string in the key field (<string.h> is included to get
*   the prototype for the strcmp function):
*      | Left margin
*      V
*      #include <string.h>
*
*      struct Test {
*        char *KeyString;
*        int   OtherFields;
*         ....
*      };
*
*      #define GENSORT_NAME                 namequicksort
*      #define GENSORT_TYPE                 struct Test
*      #define GENSORT_KEYTYPE              char *
*      #define GENSORT_GETKEY(a)            a.KeyString
*      #define GENSORT_COMPAREKEYS(k1,k2)   strcmp(k1,k2) < 0
*
*      #define GENSORT_USEPOINTERS
*
*      #include "gensort.h"
*
*  Example 3: Generate a routine to sort an indirection table 'ind' for
*   a string array 'table', i. e. 'table' stays unsorted, but 'table[ind[i]]'
*   will access the original array in sorted order
*      | Left margin
*      V
*      #include <string.h>
*
*      ...
*
*      #define GENSORT_NAME                 indquicksort
*      #define GENSORT_ARGS                 ,strings
*      #define GENSORT_ARGSPROTO            ,char **strings
*      #define GENSORT_TYPE                 int
*      #define GENSORT_KEYTYPE              char *
*      #define GENSORT_GETKEY(a)            strings[a]
*      #define GENSORT_COMPAREKEYS(k1,k2)   strcmp(k1,k2) < 0
*
*      #define GENSORT_NISINT
*
*      #include "gensort.h"
*
*      ...
*
*      {
*        char **strings;
*        int i,n,*ind;
*
*        ...
*
*        for(i = 0 ; i < n ; i++) ind[i] = i;
*
*        indquicksort(ind,n,strings);
*
*        for(i = 0 ; i < n ; i++)
*          printf("%s\n",strings[ind[i]]);
*
*      }
*
*************************************************************************/

/*  Check if all symbols we need are defined, so that user gets a
 *   message from the preprocessor in case he forgot something.
 */
#ifndef GENSORT_NAME
# error "gensort.h": Symbol "GENSORT_NAME" undefined!
#endif

#ifndef GENSORT_TYPE
# error "gensort.h": Symbol "GENSORT_TYPE" undefined!
#endif

#ifndef GENSORT_KEYTYPE
# error "gensort.h": Symbol "GENSORT_KEYTYPE" undefined!
#endif

#ifndef GENSORT_GETKEY
# error "gensort.h": Macro "GENSORT_GETKEY"  undefined!
#endif

#ifndef GENSORT_COMPAREKEYS
# error "gensort.h": Macro "GENSORT_COMPAREKEYS"  undefined!
#endif

#ifndef GENSORT_ARGSPROTO
 /*  If not given, set them to nothing  */
# define GENSORT_ARGSPROTO
# define GENSORT_ARGS
#endif

#ifndef GENSORT_LIMIT
/*  If array size is GENSORT_LIMIT or below, use straight insertion  */
# define GENSORT_LIMIT   14  /*  value is hand optimized  */
#endif

#ifndef GENSORT_INTTYPE
# ifdef GENSORT_NISINT
#  define GENSORT_INTTYPE     int
# else
#  define GENSORT_INTTYPE     long
# endif
#endif

#ifndef GENSORT_HEAPSORT
# ifndef GENSORT_USEPOINTERS
/* ******   Quicksort without pointers     *********** */
void GENSORT_NAME(GENSORT_TYPE *__array,GENSORT_INTTYPE __n GENSORT_ARGSPROTO)
{
  register GENSORT_INTTYPE __i,__j;
  GENSORT_KEYTYPE __x;

  /*   We only need n - 1 in most cases    */
  __n--;

  /*    Only sort array of at least two elements, all other are sorted  */
  while(__n > 0) {
#  if !defined(GENSORT_NOSMALL) && !defined(GENSORT_SWAP)
    /*  If __n is small enough, use straight insertion  */
    if(__n < GENSORT_LIMIT) {
      for(__j = 1 ; __j <= __n ; __j++) {
	GENSORT_TYPE __save;

#   ifdef GENSORT_GETELEMENT
	__save = GENSORT_GETELEMENT(__array,__j);
#   else
	__save = __array[__j];
#   endif

	for(__i = __j - 1 ; __i >= 0 &&
#   ifdef GENSORT_GETELEMENT
	    (GENSORT_COMPAREKEYS((GENSORT_GETKEY(__save)),
				 (GENSORT_GETKEY((GENSORT_GETELEMENT(__array,
								     __i))))))
#   else
	    (GENSORT_COMPAREKEYS((GENSORT_GETKEY(__save)),
				 (GENSORT_GETKEY((__array[__i])))))
#   endif
	      ; __i--) {
#   ifdef GENSORT_GETELEMENT
	  GENSORT_GETELEMENT(__array,(__i + 1)) =
	    GENSORT_GETELEMENT(__array,__i);
#   else
	  __array[__i + 1] = __array[__i];
#   endif
	}

#   ifdef GENSORT_GETELEMENT
	GENSORT_GETELEMENT(__array,(__i + 1)) = __save;
#   else
	__array[__i + 1] = __save;
#   endif
      }

      return;
    }
#  endif

    /*   Preset pointers   */
    __i = 0;    /*  Index to first element  */
    __j = __n;  /*  Index to last element   */

#  ifdef GENSORT_GETELEMENT
    /*      Get middle element      */
    __x = GENSORT_GETKEY(GENSORT_GETELEMENT(__array,((__n + 1) >> 1)));
#  else
    /*      Get middle element      */
    __x = GENSORT_GETKEY((__array[ (__n + 1) >> 1]));
#  endif

    do {
#  ifdef GENSORT_GETELEMENT
      while(GENSORT_COMPAREKEYS((GENSORT_GETKEY(
                (GENSORT_GETELEMENT(__array,__i)))),__x))
	__i++;

      while(GENSORT_COMPAREKEYS(__x,(GENSORT_GETKEY(
		(GENSORT_GETELEMENT(__array,__j))))))
	__j--;
#  else
      while(GENSORT_COMPAREKEYS((GENSORT_GETKEY((__array[__i]))),__x))
	__i++;

      while(GENSORT_COMPAREKEYS(__x,(GENSORT_GETKEY((__array[__j])))))
	__j--;
#  endif

      if(__i > __j) break;

#  ifdef GENSORT_SWAP
      {
        /*    Use user supplied macro   */
        GENSORT_SWAP(__array,__i,__j);
      }
      __i++; __j--;
#  else
      {
	GENSORT_TYPE __temp; /*   Only needed temporary  */

	/*     Exchange elements   */
        __temp = __array[__i];  __array[__i++] = __array[__j];
	__array[__j--] = __temp;
      }
#  endif
    } while(__i <= __j);

    if(__j < (__n - __i)) {
      /*  Call recursive with left part  */
      if(0 < __j)
	GENSORT_NAME(__array,(__j + 1) GENSORT_ARGS);

      /*  Loop again with right part   */
#  ifdef GENSORT_GETELEMENT
      __array = &GENSORT_GETELEMENT(__array,__i);
#  else
      __array += __i;
#  endif
      __n -= __i;
    } else {
      /*  Call recursive with right part  */
      if(__i < __n)
#  ifdef GENSORT_GETELEMENT
	GENSORT_NAME(&GENSORT_GETELEMENT(__array,__i),
		     (__n + 1 - __i) GENSORT_ARGS);
#  else
	GENSORT_NAME(__array + __i,__n + 1 - __i GENSORT_ARGS);
#  endif
      /*  Loop again with left part   */
      __n = __j;
    }
  }

  return;
}

# else

#  ifdef GENSORT_GETELEMENT
#   error GENSORT.H: Cannot use GENSORT_GETELEMENT with pointers!
#  endif

/* ******   Quicksort with pointers     *********** */
void GENSORT_NAME(GENSORT_TYPE *__array,GENSORT_INTTYPE __n GENSORT_ARGSPROTO)
{
  register GENSORT_TYPE *__arrayi;
  register GENSORT_TYPE *__arrayj;
  GENSORT_KEYTYPE __x;

  /*   We only need __n - 1 in most cases    */
  __n--;

  /*    Only sort array of at least two elements, all other are sorted  */
  while(__n > 0) {
#  if !defined(GENSORT_NOSMALL) && !defined(GENSORT_SWAP)
    /*  If __n is small enough, use straight insertion  */
    if(__n < GENSORT_LIMIT) {
      for(__arrayj = &__array[1] ; __arrayj <= &__array[__n] ; __arrayj++) {
	GENSORT_TYPE __save;

	__save = *__arrayj;

	for(__arrayi = __arrayj - 1 ; __arrayi >= __array &&
	    (GENSORT_COMPAREKEYS((GENSORT_GETKEY(__save)),
				 (GENSORT_GETKEY((*__arrayi)))))
	      ; __arrayi--) {
	  __arrayi[1] = *__arrayi;
	}

	__arrayi[1] = __save;
      }

      return;
    }
#  endif

    /*   Preset pointers   */
    __arrayi = __array;        /*  Pointer to first element  */
    __arrayj = __array + __n;  /*  Pointer to last element   */

    /*      Get middle element      */
    __x = GENSORT_GETKEY( (__array[ (__n + 1) >> 1]) );

    do {
      while(GENSORT_COMPAREKEYS((GENSORT_GETKEY((*__arrayi))),__x))
	__arrayi++;

      while(GENSORT_COMPAREKEYS(__x,(GENSORT_GETKEY((*__arrayj)))))
	__arrayj--;

      if(__arrayi > __arrayj)
	break;

#  ifdef GENSORT_SWAP
      {
        /*    Use user supplied macro   */
        GENSORT_SWAP(__array,(__arrayi - __array),(__arrayj - __array));
      }
      __arrayi++; __arrayj--;
#  else
      {
	GENSORT_TYPE __temp; /*   Only needed temporary  */

	/*     Exchange elements   */
	__temp = *__arrayi;  *__arrayi++ = *__arrayj;  *__arrayj-- = __temp;
      }
#  endif
    } while(__arrayi <= __arrayj);

    if( (__arrayj - __array) < (__array - __arrayi + __n)) {
      /*  Call recursive with left part  */
      if(__array < __arrayj)
	GENSORT_NAME(__array,__arrayj - __array + 1 GENSORT_ARGS);

      /*  Loop again with right part   */
      __n -= __arrayi - __array;
      __array = __arrayi;
    } else {
      /*  Call recursive with right part  */
      if(__arrayi < (__array + __n))
	GENSORT_NAME(__arrayi,__array + __n + 1 - __arrayi GENSORT_ARGS);

      /*  Loop again with left part   */
      __n = __arrayj - __array;
    }
  }

  return;
}
# endif

#else
/* ******   Heapsort routine              *********** */

#ifdef GENSORT_SWAP
# error "gensort.h": Cannot use GENSORT_SWAP and GENSORT_HEAPSORT!
#endif

/*   Preprocessor macro to generate second routine    */
/*    with temporary name (sift routine).             */
#define GENSORT_TEMPNAME(name)   GENSORT_TEMPNAME2(name)
#define GENSORT_TEMPNAME2(name)  name##__sift

/* ****    Now the main heapsort routine  **** */
void GENSORT_NAME(GENSORT_TYPE *__array,GENSORT_INTTYPE __n GENSORT_ARGS)
{
  register GENSORT_INTTYPE __i;

  /*   Prototype for sift routine (Compiler does not like "static" here!  */
  void GENSORT_TEMPNAME(GENSORT_NAME)
    (GENSORT_TYPE *,GENSORT_INTTYPE,GENSORT_INTTYPE GENSORT_ARGS);

  /*  Build up heap    */
  for(__i = (__n >> 1) ; __i > 0 ; )
    GENSORT_TEMPNAME(GENSORT_NAME) (__array,--__i,__n GENSORT_ARGS);

  for(__i = __n - 1 ; __i > 0 ; ) {
    {
      GENSORT_TYPE __temp;

#  ifdef GENSORT_GETELEMENT
      __temp = GENSORT_GETELEMENT(__array,0);
      GENSORT_GETELEMENT(__array,0) = GENSORT_GETELEMENT(__array,__i);
      GENSORT_GETELEMENT(__array,i) = __temp;
#  else
      /*     Exchange elements   */
      __temp = *__array;  *__array = __array[__i];  __array[__i] = __temp;
#  endif
    }

    /*   Call sift routine    */
    GENSORT_TEMPNAME(GENSORT_NAME) (__array,0,__i-- GENSORT_ARGS);
  }

  return;
}

/* ****   Sift routine is global! Compiler does not like static here  **** */
void GENSORT_TEMPNAME(GENSORT_NAME) (GENSORT_TYPE *__array,GENSORT_INTTYPE __i,
				     GENSORT_INTTYPE __n GENSORT_ARGS)
{
  register GENSORT_INTTYPE __j;
  GENSORT_TYPE __temp;

  /*   Get first element to sink in heap   */
# ifdef GENSORT_GETELEMENT
  __temp = GENSORT_GETELEMENT(__array,__i);
# else
  __temp = __array[__i];
# endif

  while( (__j = ((__i << 1) + 1)) < __n) {
# ifdef GENSORT_GETELEMENT
    if( (__j < (__n - 1)) && (GENSORT_COMPAREKEYS((GENSORT_GETKEY(
		GENSORT_GETELEMENT(__array,__j))),(GENSORT_GETKEY(
		GENSORT_GETELEMENT(__array,(__j + 1)))))) )
      ++__j;

    if( !(GENSORT_COMPAREKEYS((GENSORT_GETKEY(__temp)),(GENSORT_GETKEY(
		GENSORT_GETELEMENT(__array,__j))))) )
      break;

    /*   Copy lower element to top    */
    GENSORT_GETELEMENT(__array,__i) = GENSORT_GETELEMENT(__array,__j);
    __i = __j;
# else
    if( (__j < (__n - 1)) && (GENSORT_COMPAREKEYS((GENSORT_GETKEY((__array[__j]))),
            (GENSORT_GETKEY((__array[__j + 1]))))) )
      ++__j;

    if( !(GENSORT_COMPAREKEYS((GENSORT_GETKEY(__temp)),
	   (GENSORT_GETKEY((__array[__j]))))) )
      break;

    /*   Copy lower element to top   */
    __array[__i] = __array[__j];
    __i = __j;
# endif
  }

# ifdef GENSORT_GETELEMENT
  /*   Now the first element has reached its place   */
  GENSORT_GETELEMENTS(__array,__i) = __temp;
# else
  /*   Now the first element has reached its place   */
  __array[__i] = __temp;
# endif

  return;
}

# undef GENSORT_TEMPNAME
# undef GENSORT_TEMPNAME2
#endif

/* *******  Undefine all symbols for reuse    ****** */
#undef GENSORT_NAME
#undef GENSORT_TYPE
#undef GENSORT_KEYTYPE
#undef GENSORT_GETKEY
#undef GENSORT_COMPAREKEYS
#undef GENSORT_HEAPSORT
#undef GENSORT_ARGSPROTO
#undef GENSORT_ARGS
#undef GENSORT_LIMIT
#undef GENSORT_NOSMALL
#undef GENSORT_SWAP
#undef GENSORT_GETELEMENT
#undef GENSORT_NISINT
#undef GENSORT_INTTYPE
