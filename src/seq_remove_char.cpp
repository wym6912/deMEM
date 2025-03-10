#include "lib/seqio.hpp"

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Error: this program has only three parameters: %s in out\n", argv[0]);
		return 1;
	}
	FILE* in = fopen(argv[1], "r"), *out = fopen(argv[2], "wb");
	if (in == NULL || out == NULL)
	{
		fprintf(stderr, "Error: can not open file %s or %s.\n Please check these files and try again. Program will exit.\n", argv[1], argv[2]);
		return 2;
	}
	kseq_t* stdin_kseq = kseq_init(fileno(in));
	char* seq_temp = NULL, *seqbk = NULL, *raw = NULL;
	while (kseq_read(stdin_kseq) >= 0)
	{
		seq_temp = (char*)malloc(sizeof(char) * (stdin_kseq->seq.l + 2));
		if (seq_temp == NULL) { fprintf(stderr, "Error: can not detect chars in sequence file. Program will exit.\n"); return 1; }
		seqbk = seq_temp;
		raw = stdin_kseq->seq.s;
		for(int i = 0; i < stdin_kseq->seq.l; ++  i, ++ raw)
			if (isalpha(*raw) || *raw == '-' || *raw == '*' || *raw == '.')
			{
				*seq_temp++ = *raw;
			}
		*seq_temp = 0;
		fprintf(out, ">%s%s\n%s\n", stdin_kseq->name.s, stdin_kseq->comment.s ? stdin_kseq->comment.s : "", seqbk);
		free(seqbk);
	}
	fcloseall();
	return 0;
}