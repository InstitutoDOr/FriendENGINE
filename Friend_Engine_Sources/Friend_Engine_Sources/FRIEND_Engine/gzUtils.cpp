#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gzUtils.h"

char prog[] = "gzUtils";
void error(const char *msg)
{
	fprintf(stderr, "%s : %s\n", prog, msg);
//	exit(1);
}

void gz_compress(FILE *in, gzFile out)
{
	char buf[BUFLEN];
	int len;
	int err;

#ifdef USE_MMAP
	/* Try first compressing with mmap. If mmap fails (minigzip used in a
	* pipe), use the normal fread loop.
	*/
	if (gz_compress_mmap(in, out) == Z_OK) return;
#endif
	for (;;) {
		len = (int)fread(buf, 1, sizeof(buf), in);
		if (ferror(in)) {
			perror("fread");
			exit(1);
		}
		if (len == 0) break;

		if (gzwrite(out, buf, (unsigned)len) != len) error(gzerror(out, &err));
	}
	fclose(in);
	if (gzclose(out) != Z_OK) error("failed gzclose");
}

#ifdef USE_MMAP /* MMAP version, Miguel Albrecht <malbrech@eso.org> */

/* Try compressing the input file at once using mmap. Return Z_OK if
* if success, Z_ERRNO otherwise.
*/
int gz_compress_mmap(in, out)
FILE   *in;
gzFile out;
{
	int len;
	int err;
	int ifd = fileno(in);
	caddr_t buf;    /* mmap'ed buffer for the entire input file */
	off_t buf_len;  /* length of the input file */
	struct stat sb;

	/* Determine the size of the file, needed for mmap: */
	if (fstat(ifd, &sb) < 0) return Z_ERRNO;
	buf_len = sb.st_size;
	if (buf_len <= 0) return Z_ERRNO;

	/* Now do the actual mmap: */
	buf = mmap((caddr_t)0, buf_len, PROT_READ, MAP_SHARED, ifd, (off_t)0);
	if (buf == (caddr_t)(-1)) return Z_ERRNO;

	/* Compress the whole file at once: */
	len = gzwrite(out, (char *)buf, (unsigned)buf_len);

	if (len != (int)buf_len) error(gzerror(out, &err));

	munmap(buf, buf_len);
	fclose(in);
	if (gzclose(out) != Z_OK) error("failed gzclose");
	return Z_OK;
}
#endif /* USE_MMAP */

/* ===========================================================================
* Uncompress input to output then close both files.
*/
void gz_uncompress(gzFile in, FILE *out)
{
	local char buf[BUFLEN];
	int len;
	int err;

	for (;;) {
		len = gzread(in, buf, sizeof(buf));
		if (len < 0) error(gzerror(in, &err));
		if (len == 0) break;

		if ((int)fwrite(buf, 1, (unsigned)len, out) != len) {
			error("failed fwrite");
		}
	}
	if (fclose(out)) error("failed fclose");

	if (gzclose(in) != Z_OK) error("failed gzclose");
}


/* ===========================================================================
* Compress the given file: create a corresponding .gz file and remove the
* original.
*/
void file_compress(char  *file, char  *mode)
{
	local char outfile[MAX_NAME_LEN];
	FILE  *in;
	gzFile out;

	if (strlen(file) + strlen(GZ_SUFFIX) >= sizeof(outfile)) {
		fprintf(stderr, "%s: filename too long\n", prog);
		exit(1);
	}

	strcpy(outfile, file);
	strcat(outfile, GZ_SUFFIX);

	in = fopen(file, "rb");
	if (in == NULL) {
		perror(file);
		exit(1);
	}
	out = gzopen(outfile, mode);
	if (out == NULL) {
		fprintf(stderr, "%s: can't gzopen %s\n", prog, outfile);
		exit(1);
	}
	gz_compress(in, out);

	unlink(file);
}


/* ===========================================================================
* Uncompress the given file and remove the original.
*/
int file_uncompress(char  *file, char *outputFile)
{
	local char buf[MAX_NAME_LEN];
	char *infile, *outfile;
	FILE  *out;
	gzFile in;
	size_t len = strlen(file);

	if (len + strlen(GZ_SUFFIX) >= sizeof(buf)) {
		fprintf(stderr, "%s: filename too long\n", prog);
		return(1);
	}

	strcpy(buf, file);

	if (len > SUFFIX_LEN && strcmp(file + len - SUFFIX_LEN, GZ_SUFFIX) == 0) {
		infile = file;
		outfile = buf;
		outfile[len - 3] = '\0';
	}
	else {
		outfile = file;
		infile = buf;
		strcat(infile, GZ_SUFFIX);
	}
	in = gzopen(infile, "rb");
	if (in == NULL) {
		fprintf(stderr, "%s: can't gzopen %s\n", prog, infile);
		return(1);
	}
	out = fopen(outfile, "wb");
	if (out == NULL) {
		perror(file);
		return(1);
	}

	gz_uncompress(in, out);
	strcpy(outputFile, outfile);
	//unlink(infile);
	return 0;
}
