/* Bi-directional linear search - iterative, returns the array dim of the closest element, assume array sorted small->large - also bidirectionally*/

#include <stdio.h>


int LinearSearch (double *array,int startindx, int stopindx, double key,int *count, int direction){
	//Iterative implementation - assume elements sorted
int i;
int keyindex;

if(direction>=0){
	keyindex=startindx;
	for (i=startindx;i<=stopindx;i++){
		(*count)++;
		if(array[i]>key) {	
			keyindex=i-1;
			return keyindex;
		}
	}
}
else if(direction<0 ){
	keyindex=stopindx;
	for(i=stopindx;i>=startindx;i--){
		(*count)++;
		if (array[i]>key){
			 keyindex=i+1;
			 return keyindex;
		}
	}
}

else printf("Error: startindex=stopindex");

}
			

