#ifndef _CUPMALLOC_H_
#define _CUPMALLOC_H_

#ifdef __cplusplus
extern "C" {
#endif

void pmalloc(void ** /*h_ptr*/, size_t /*nbytes*/);
void pmalloc_wc(void **h_ptr, size_t nbytes);
void pfree(void * /*h_ptr*/); 

#ifdef __cplusplus
}
#endif
#endif  /* __CUPMALLOC_H_ */


