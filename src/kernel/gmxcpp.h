/* The possible return codes for these functions */
enum { eCPP_OK, eCPP_FILE_NOT_FOUND, eCPP_EOF, eCPP_SYNTAX, eCPP_INTERRUPT,
       eCPP_INVALID_HANDLE,
       eCPP_FILE_NOT_OPEN, eCPP_UNKNOWN, eCPP_NR };

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. Return integer status */
extern int cpp_open_file(char *filenm,void **handlep);

/* Turn debugging (printing to stderr) on or off. */
extern void cpp_debug_on(void);
extern void cpp_debug_off(void);
   
/* Return one whole line from the file into buf which holds at most n
   characters, for subsequent processing. Returns integer status */
extern int cpp_read_line(void **handlep,int n,char buf[]);

/* Close the file! Return integer status. */
extern int cpp_close_file(void **handlep);

/* Return a string containing the error message coresponding to status
   variable */
extern char *cpp_error(void **handlep,int status);
