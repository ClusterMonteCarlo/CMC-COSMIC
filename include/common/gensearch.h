/* gensearch.h: binary search routine generator
 * 
 * inspired by Joerg Schoen's gensort.h
 *
 * you need to define the following before #include'ing this header file:
 * 
 * GENSEARCH_NAME        : name of the routine generated
 * GENSEARCH_TYPE        : type of the array being searched
 * GENSEARCH_KEYTYPE     : type of the key being used for searching
 * GENSEARCH_GETKEY      : a method for obtaining the key from an element 
 *                        of the array
 * GENSEARCH_COMPAREKEYS : a method for comparing two keys
 *
 * example 1, a binary search routine for a double array:
 *
 * #define GENSEARCH_NAME 			bsearchex1
 * #define GENSEARCH_TYPE 			double
 * #define GENSEARCH_KEYTYPE			double
 * #define GENSEARCH_GETKEY(a)		a
 * #define GENSEARCH_COMPAREKEYS(k1, k2)	k1 < k2
 *
 * example 2, a binary search routine for a struct, sorted with a long int key:
 * 
 * struct {
 * 	double	strvalue;
 * 	long int	strkey;
 * } mystr;
 * 
 * #define GENSEARCH_NAME 			bsearchex2
 * #define GENSEARCH_TYPE 			struct mystr
 * #define GENSEARCH_KEYTYPE			long int
 * #define GENSEARCH_GETKEY(a)		a.strkey
 * #define GENSEARCH_COMPAREKEYS(k1, k2)	k1 < k2

 */
size_t GENSEARCH_NAME(GENSEARCH_TYPE *__x, size_t __min, size_t __max, GENSEARCH_KEYTYPE __y) {
	size_t __try;

	if ( GENSEARCH_COMPAREKEYS(GENSEARCH_GETKEY(__x[__min]), 
				         GENSEARCH_GETKEY(__x[__max])) ) {
		do {
			__try = (__min+__max+1)/2;
			if ( GENSEARCH_COMPAREKEYS(GENSEARCH_GETKEY(__x[__try]), __y) ) {
				__min = __try;
			} else {
				__max = __try-1;
			} 
		} while (__max != __min);
	} else {
		do {
			__try = (__min+__max+1)/2;
			if ( GENSEARCH_COMPAREKEYS(__y, GENSEARCH_GETKEY(__x[__try])) ) {
				__min = __try;
			} else {
				__max = __try-1;
			} 
		} while (__max != __min);
	}

	return __min;
}
