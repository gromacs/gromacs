

/* generate a unique array of size size, based on rank, N and testnr. */
void unique_array(int rank, int N, int testnr, int size, int *array);

/* check two arrays */
int check_arrays(int *a, int *b, int size);

/* print and exit */
void print_match_err(int myrank, int *buf_recv, int *buf_expt, int size);



void send_recv_data_tester(void);
