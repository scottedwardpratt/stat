/*********************************************************************
process_pipe.h
Copyright 2012, The University of North Carolina at Chapel Hill.

This software was written in 2012 by Hal Canary <cs.unc.edu/~hal>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/

#ifndef __process_pipe_h_
#define __process_pipe_h_

#ifdef __cplusplus
#include <cstdio>
namespace madai {
	extern "C" {
		typedef struct process_pipe {
			std::FILE * question;
			std::FILE * answer;
			long int pid; /* in case other signals need to be sent. */
		} process_pipe;

#else /* NOT __cplusplus */

#include <stdio.h>
typedef struct process_pipe {
	FILE * question;
	FILE * answer;
	long int pid; /* in case other signals need to be sent. */
} process_pipe;

#endif /* __cplusplus */


/** returns EXIT_FAILURE on error, EXIT_SUCCESS otherwise.  Note that
 argv[0] doesn't seem to respect your PATH evironment variable, so you
 may need to pass an absolute path. argv should be NULL terminated. */
		int create_process_pipe(process_pipe * pp, char * const * argv);

#ifdef __cplusplus
	}
}
#endif /* __cplusplus */

#endif	/* __process_pipe_h_ */
