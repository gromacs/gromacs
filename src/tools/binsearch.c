#include <types/simple.h>

/*Make range-array (Permutation identity) for sorting */
void rangeArray(int *ar,int size)
{
int i;	
	for (i=0;i<size;i++){
		ar[i]=i;
	}
}

void pswap(int *v1, int *v2)
{
	int temp;
	temp=*v1;
	*v1=*v2;
	*v2=temp;
}


void Swap (real *v1, real *v2)
{
real temp;
temp = *v1;
*v1 = *v2;
*v2 = temp;
}



void insertionSort(real *arr, int *perm, int startndx, int endndx, int direction) {
      int i, j;

	if(direction>=0){
      		for (i = startndx; i <=endndx; i++) {
            		j = i;

            		while (j > startndx && arr[j - 1] > arr[j]) {
                  		Swap(&arr[j],&arr[j-1]);
				pswap(&perm[j],&perm[j-1]);
                  		j--;
            		}

      		}

	}

	if(direction<0){
		 for (i = startndx; i <=endndx; i++) {
                        j = i;

                        while (j > startndx && arr[j - 1] < arr[j]) {
                                Swap(&arr[j],&arr[j-1]);
				pswap(&perm[j],&perm[j-1]);
                                j--;
                        }

                }
	}
}


int BinarySearch (real *array, int low, int high, real key,int direction)
{
int mid, max, min;
max=high+2;
min=low+1;

//Iterative implementation

if (direction>=0){
while (max-min>1){
	mid=(min+max)>>1;
	if(key<array[mid-1]) max=mid;
	else min=mid;
}
return min;
}

else if (direction<0){
while(max-min>1){
	mid=(min+max)>>1;
	if(key>array[mid-1]) max=mid;
	else min=mid;
}
return min-1;

    }//end -ifelse direction
}


int start_binsearch(real *array, int *perm, int low, int high, real key, int direction)
{
	insertionSort(array,perm,low,high,direction);
	return BinarySearch(array,low,high,key,direction);
}
