/* freadmat.h May 28 2001 H.Mine */

double **fread_mat_lls(FILE *fp, int *mp, int *np);
double **fread_mat_lfh(FILE *fp, int *mp, int *np);
double **fread_mat_paup(FILE *fp, int *mp, int *np);
double **fread_mat_puzzle(FILE *fp, int *mp, int *np);
double **fread_mat_phyml(FILE *fp, int *mp, int *np);

enum seqfile {SEQ_MT, SEQ_MOLPHY, SEQ_PAML, SEQ_PAUP, SEQ_PUZZLE, SEQ_PHYML};

int seqmode=SEQ_MT;

char *fext_mt=".mt";
char *fext_molphy=".lls";
char *fext_paml=".lnf";
char *fext_paup=".txt";
char *fext_puzzle=".sitelh";
char *fext_phyml=".txt";
