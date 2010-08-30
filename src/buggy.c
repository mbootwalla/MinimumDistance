#include <Rinternals.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

/* XXX: this is a version of the alpha_freq code with bugs introduced so
   that we can demo the use of gdb
 */

static int update_freq(const char *file, int *freq_data)
{
    FILE *fh = fopen(file, "r");
    int c, res;
    while (EOF != (c = getc(fh))) {
        int i = tolower(c) - 'a';
        if (i > -1 && i < 26) freq_data[i] += 1;
    }
    res = ferror(fh);
    fclose(fh);
    return res;
}

SEXP sapply_mean(SEXP fname)
{
    SEXP ans = Rf_allocVector(INTSXP, 26);
    const char *file;
    int *freq_data;

    PROTECT(ans);

    freq_data = INTEGER(ans);
    memset(freq_data, 0, sizeof(int) * 26);
    file  = CHAR(VECTOR_ELT(fname, 0));
    update_freq(NULL, freq_data);

    UNPROTECT(1);
    return ans;
}
